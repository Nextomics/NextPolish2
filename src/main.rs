use crossbeam_channel::bounded;
use crossbeam_utils::thread;
use fxhash::FxHashMap as HashMap;
use kseq::parse_path;
use rust_htslib::bam::{
    ext::BamRecordExtensions, record::Cigar, record::CigarStringView, IndexedReader, Read, Record,
};
use std::cell::Cell;
use std::cmp::{max, Reverse};
use std::fmt;
use std::iter::zip;

mod utils;
use utils::{
    kmer::{iter2kmer, KmerInfo, SEQ_NUM},
    louvain::{assign_data, insert_data, new_data, phase_communities},
    option::Option as Opt,
    resource::resource_str,
    secondary::{retrieve_secondary_seq_from_bam, reverse_complement_seq_u8},
};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

const LQSEQ_MAX_CAN_COUNT: usize = 60;
const INVALID_KMER: u64 = u64::MAX;

#[derive(Default, Debug, Clone, Copy, PartialEq)]
struct AlignBase {
    q_base: u8,
    delta: u16,
    t_pos: u32,
} // an aligned base

impl AlignBase {
    fn head(t_pos: u32, delta: u16) -> Self {
        Self {
            q_base: 0b1111, //Max to 4 bits, should be different from valid bases
            delta,
            t_pos,
        }
    }

    fn is_head(&self) -> bool {
        self.q_base == 0b1111
    }
}

type Kctype = u32;
#[derive(Clone, Debug, Default)]
struct Kmer {
    // 3-mer
    delta: u16, // the delta position of 1-base

    // bases: u16 = u2-u2-u4-u4-u4
    //    u2         u2       u4     u4     u4
    //    --         --      ----   ----   ----
    //  2-base     3-base   1-base 2-base 3-base
    //delta-flag delta-flag
    bases: u16,
    count: Kctype,       //count of this 3-mer
    besti: Cell<Kctype>, //index of the previous best 3-mer
    score: Cell<i64>,    //score of this 3-mer
}

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}{}{}",
            SEQ_NUM[(self.bases >> 8 & 0b1111) as usize] as char,
            SEQ_NUM[(self.bases >> 4 & 0b1111) as usize] as char,
            SEQ_NUM[(self.bases & 0b1111) as usize] as char
        )
    }
}

impl Kmer {
    fn new(base1: AlignBase, base2: AlignBase, base3: AlignBase) -> Self {
        let mut bases = 0u16;
        if base2.t_pos == base1.t_pos {
            bases |= 0b0100;
        }

        if base2.t_pos == base3.t_pos {
            bases |= 0b0001;
        }

        Self {
            delta: base1.delta,
            bases: ((bases << 4 | base1.q_base as u16) << 4 | base2.q_base as u16) << 4
                | base3.q_base as u16,
            count: 1,
            besti: Cell::new(0),
            score: Cell::new(0),
        }
    }

    // p is the position of 3-base
    fn bases(&self, p: u32) -> (AlignBase, AlignBase, AlignBase) {
        //the 1-3 bases
        if self.bases & 0b0101_0000_0000_0000 == 0b0101_0000_0000_0000 {
            //A--
            (
                AlignBase {
                    t_pos: p,
                    delta: self.delta,
                    q_base: (self.bases >> 8 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p,
                    delta: self.delta + 1,
                    q_base: (self.bases >> 4 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p,
                    delta: self.delta + 2,
                    q_base: (self.bases & 0b1111) as u8,
                },
            )
        } else if self.bases & 0b0001_0000_0000_0000 != 0 {
            //AA-
            (
                AlignBase {
                    t_pos: p - 1,
                    delta: self.delta,
                    q_base: (self.bases >> 8 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p,
                    delta: 0,
                    q_base: (self.bases >> 4 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p,
                    delta: 1,
                    q_base: (self.bases & 0b1111) as u8,
                },
            )
        } else if self.bases & 0b0100_0000_0000_0000 != 0 {
            //A-A
            (
                AlignBase {
                    t_pos: p - 1,
                    delta: self.delta,
                    q_base: (self.bases >> 8 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p - 1,
                    delta: self.delta + 1,
                    q_base: (self.bases >> 4 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p,
                    delta: 0,
                    q_base: (self.bases & 0b1111) as u8,
                },
            )
        } else {
            //AAA
            (
                AlignBase {
                    t_pos: p - 2,
                    delta: self.delta,
                    q_base: (self.bases >> 8 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p - 1,
                    delta: 0,
                    q_base: (self.bases >> 4 & 0b1111) as u8,
                },
                AlignBase {
                    t_pos: p,
                    delta: 0,
                    q_base: (self.bases & 0b1111) as u8,
                },
            )
        }
    }
}

#[derive(Clone, Default)]
struct Msa {
    kmers: Vec<Kmer>,
}

impl Msa {
    fn push(&mut self, value: Kmer) {
        let mut is_exist = false;
        for kmer in self.kmers.iter_mut() {
            if kmer.bases == value.bases && kmer.delta == value.delta {
                kmer.count += 1;
                assert!(kmer.count < Kctype::MAX, "kmer count overflow!");
                is_exist = true;
                break;
            }
        }

        if !is_exist {
            self.kmers.push(value);
        }
    }

    fn get(&self, base2: AlignBase, base3: AlignBase) -> impl Iterator<Item = (Kctype, &Kmer)> {
        let base23 = base2.q_base << 4 | base3.q_base;
        let delta23 = u16::from(base2.t_pos == base3.t_pos);
        self.kmers.iter().enumerate().filter_map(move |(i, v)| {
            if v.bases as u8 == base23 && v.bases >> 12 & 1 == delta23 {
                let bases = v.bases(base3.t_pos);
                if bases.1 == base2 && bases.2 == base3 {
                    return Some((
                        i.try_into()
                            .expect("The length of kmers in Msa > Kctype::MAX."),
                        v,
                    ));
                }
            }
            None
        })
    }

    fn sort(&mut self) {
        self.kmers.sort_by_cached_key(|v| v.bases(0).2.delta);
    }

    //sort() must be called before calling this function
    fn coverage(&self) -> i64 {
        let mut c = 0;
        for kmer in &self.kmers {
            if kmer.bases(0).2.delta != 0 {
                break;
            }
            c += kmer.count as i64;
        }
        c
    }

    fn shrink(&mut self) {
        self.kmers.shrink_to_fit();
    }

    fn clear(&mut self) {
        self.kmers.clear();
    }
}

fn clear_msas(msas: &mut [Msa]) {
    for msa in msas {
        msa.clear();
    }
}

fn shrink_msas(msas: &mut [Msa]) {
    for msa in msas {
        msa.shrink();
    }
}

fn sort_msas(msas: &mut [Msa]) {
    for msa in msas {
        msa.sort();
    }
}

// the field aln_t_s always < 1<<30, so its highest bit can be used as a lable
const ALN_T_S_LABLE: u32 = 1 << 31;
struct AlignSeq {
    aln_t_s: u32,
    aln_t_e: u32,
    align_bases: Vec<u8>,
}

impl AlignSeq {
    fn new(aln: &Alignment) -> Self {
        let len = (aln.aln_len() + 1) >> 1;
        let mut alignseq = AlignSeq {
            aln_t_s: aln.aln_t_s,
            aln_t_e: aln.aln_t_s,
            align_bases: vec![0; len + 1],
        };

        let mut b;
        let mut i = 0;
        for (tb, qb) in
            zip(aln.t_aln_str.as_bytes(), aln.q_aln_str.as_bytes()).skip(aln.shift as usize)
        {
            b = SEQ_NUM[*qb as usize];
            if *tb == b'-' {
                b |= 8;
            } else if i != 0 {
                alignseq.aln_t_e += 1;
            }

            if i & 1 == 0 {
                b <<= 4;
            }
            alignseq.align_bases[i >> 1] |= b;
            i += 1;
        }

        if i & 1 == 0 {
            alignseq.align_bases[i >> 1] |= 255;
        } else {
            alignseq.align_bases[i >> 1] |= 15;
        }
        alignseq
    }

    fn get_align_tag(&self, p: &mut usize, alignbase: &mut AlignBase) -> bool {
        let mut t: u8 = self.align_bases[*p >> 1];
        if *p & 1 == 0 {
            t >>= 4;
        }

        if (t & 15) == 15 {
            return false;
        }

        alignbase.q_base = t & 7;
        if *p != 0 {
            if t & 8 != 0 {
                alignbase.delta += 1;
            } else {
                alignbase.delta = 0;
                alignbase.t_pos += 1;
            }
        } else {
            alignbase.t_pos = self.aln_t_s;
            alignbase.delta = 0;
        }
        *p += 1;
        true
    }

    fn set_lable(&mut self) {
        self.aln_t_s |= ALN_T_S_LABLE;
    }

    fn has_lable(&self) -> bool {
        self.aln_t_s & ALN_T_S_LABLE != 0
    }

    fn unset_lable(&mut self) {
        self.aln_t_s ^= ALN_T_S_LABLE;
    }
}

#[derive(Default, Debug)]
struct Alignment {
    shift: u32,
    aln_t_s: u32,      // ref.
    aln_t_e: u32,      //not include, ref.
    aln_q_s: u32,      // read
    aln_q_e: u32,      //not include, read
    q_aln_str: String, // read
    t_aln_str: String, // ref.
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "shift: {:>6}\ntseq: {:>7}{:>7} {}\nqseq: {:>7}{:>7} {}\n",
            self.shift,
            self.aln_t_s,
            self.aln_t_e,
            &self.t_aln_str[self.shift as usize..],
            self.aln_q_s,
            self.aln_q_e,
            &self.q_aln_str[self.shift as usize..]
        )
    }
}

impl Alignment {
    fn new() -> Self {
        Default::default()
    }

    // self.aln_t_s must be set
    fn fill_with_cigar<T>(&mut self, cigars: &CigarStringView, tseq: &str, qseq: &T)
    where
        T: std::ops::Index<usize, Output = u8>,
    {
        let mut qs = 0; // read
        let mut ts = 0; // ref.
        let mut is_first = true;
        for c in cigars {
            match c {
                Cigar::SoftClip(l) => {
                    qs += l;
                    if is_first {
                        self.aln_q_s = qs;
                    } else {
                        self.aln_q_e = qs - l;
                    }
                }
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                    (0..*l).for_each(|_| {
                        self.q_aln_str.push(qseq[qs as usize] as char);
                        qs += 1;
                    });
                    self.t_aln_str
                        .push_str(&tseq[ts as usize..(ts + l) as usize]);
                    ts += l;
                }
                Cigar::Ins(l) => {
                    (0..*l).for_each(|_| {
                        self.q_aln_str.push(qseq[qs as usize] as char);
                        qs += 1;
                    });
                    std::iter::repeat('-')
                        .take(*l as usize)
                        .for_each(|v| self.t_aln_str.push(v));
                }
                Cigar::Del(l) => {
                    std::iter::repeat('-')
                        .take(*l as usize)
                        .for_each(|v| self.q_aln_str.push(v));
                    self.t_aln_str
                        .push_str(&tseq[ts as usize..(ts + l) as usize]);
                    ts += l;
                }
                Cigar::HardClip(_) => (),
                _ => {
                    panic!("Unknown cigar: {:?}", c);
                }
            }
            is_first = false;
        }
        if self.aln_q_e == 0 {
            self.aln_q_e = qs;
        }
        self.aln_t_e = self.aln_t_s + ts;
    }

    fn aln_len(&self) -> usize {
        self.t_aln_str.len() - self.shift as usize
    }

    //ensure alignment start/end with a matched len-mer
    fn trim(&mut self, len: u32) {
        let mut j = 0;
        let t_aln_str = self.t_aln_str.as_bytes();
        let q_aln_str = self.q_aln_str.as_bytes();

        // start
        for i in 0..self.t_aln_str.len() {
            if t_aln_str[i] == q_aln_str[i] {
                j += 1;
                self.aln_t_s += 1;
                self.aln_q_s += 1;
            } else {
                if t_aln_str[i] != b'-' {
                    self.aln_t_s += 1;
                }

                if q_aln_str[i] != b'-' {
                    self.aln_q_s += 1;
                }

                j = 0;
            }

            if j == len {
                self.aln_t_s -= len;
                self.aln_q_s -= len;
                self.shift = i as u32 + 1 - len;
                break;
            }
        }

        // end
        if j == len {
            j = 0;
            for i in (0..self.t_aln_str.len()).rev() {
                if t_aln_str[i] == q_aln_str[i] {
                    j += 1;
                    self.aln_t_e -= 1;
                    self.aln_q_e -= 1;
                } else {
                    if t_aln_str[i] != b'-' {
                        self.aln_t_e -= 1;
                    }

                    if q_aln_str[i] != b'-' {
                        self.aln_q_e -= 1;
                    }

                    j = 0;
                }

                if j == len {
                    self.aln_t_e += len;
                    self.aln_q_e += len;

                    let new_len = i + len as usize;
                    if new_len < self.t_aln_str.len() {
                        self.t_aln_str.truncate(new_len);
                        self.q_aln_str.truncate(new_len);
                    }
                    break;
                }
            }
        } else {
            self.shift = self.t_aln_str.len() as u32;
        }
    }

    fn clear(&mut self) {
        self.shift = 0;
        self.aln_t_s = 0;
        self.aln_t_e = 0;
        self.aln_q_s = 0;
        self.aln_q_e = 0;
        self.q_aln_str.clear(); // read
        self.t_aln_str.clear(); // ref.
    }

    fn shrink_to(&mut self, min_capacity: usize) {
        self.t_aln_str.shrink_to(min_capacity);
        self.q_aln_str.shrink_to(min_capacity);
    }
}

fn filter_alignseqs_by_clip(alignseqs: &mut [AlignSeq]) {
    fn in_ranges(aln_ranges: &[(u32, u32)], start: u32, end: u32) -> bool {
        for (s, e) in aln_ranges {
            if *s <= start && end <= *e {
                return true;
            } else if end < *s {
                break;
            }
        }
        false
    }

    let offset = 50;
    let mut aln_ranges = Vec::with_capacity(1024);
    let (mut s, mut e) = (0, 0);
    for (aln_t_s, aln_t_e) in alignseqs
        .iter()
        .filter(|x| !x.has_lable())
        .map(|x| (x.aln_t_s + offset, x.aln_t_e - offset))
    {
        if s == e {
            s = aln_t_s;
            e = aln_t_e;
        } else if aln_t_s > e {
            aln_ranges.push((s, e));
            s = aln_t_s;
            e = aln_t_e;
        } else if e < aln_t_e {
            e = aln_t_e;
        }
    }

    if s != e {
        aln_ranges.push((s, e));
    }

    for alignseq in alignseqs.iter_mut().filter(|x| x.has_lable()) {
        alignseq.unset_lable();
        if in_ranges(&aln_ranges, alignseq.aln_t_s, alignseq.aln_t_e) {
            // println!("{:?} {:?}", alignseq.aln_t_s, alignseq.aln_t_e);
            alignseq.align_bases = Vec::new();
        }
    }
}

fn update_msas(msas: &mut [Msa], alignseqs: &[AlignSeq]) {
    for alignseq in alignseqs.iter().filter(|x| !x.align_bases.is_empty()) {
        let mut p = 0;
        let mut b1 = AlignBase::head(alignseq.aln_t_s - 1, 0);
        let mut b2 = AlignBase::head(alignseq.aln_t_s - 1, 1);
        let mut b3: AlignBase = Default::default();
        while alignseq.get_align_tag(&mut p, &mut b3) {
            let msa = &mut msas[b3.t_pos as usize];
            msa.push(Kmer::new(b1, b2, b3));
            b1 = b2;
            b2 = b3;
        }
    }
}

#[derive(Clone, Copy)]
struct ConsensusBase {
    pos: u32,
    // qv: char,
    base: char,
}

fn display_consensusbase_vec(
    tid: &str,
    consensusbases: &[ConsensusBase],
    uppercase: bool,
    out_pos: bool,
) {
    if out_pos {
        for base in consensusbases {
            println!(
                "{}\t{}\t{}",
                tid,
                if uppercase {
                    base.base.to_ascii_uppercase()
                } else {
                    base.base
                },
                base.pos
            );
        }
    } else {
        println!(
            ">{} start:{} end:{}",
            tid,
            consensusbases.first().unwrap().pos,
            consensusbases.last().unwrap().pos
        );
        for base in consensusbases.iter().map(|x| {
            if uppercase {
                x.base.to_ascii_uppercase()
            } else {
                x.base
            }
        }) {
            print!("{}", base);
        }

        println!();
    }
}

#[derive(Default, Debug)]
struct LqSeq {
    order: u32,
    kscore: u16,
    kmer: u64,   // 0 means not a valid kmer
    seq: String,
}

const LQSEQS_LABLE_TEMP: u8 = 0b0000_0001; //a temporary flag
const LQSEQS_LABLE_SUCC: u8 = 0b1000_0000; //LqSeqs are successfully corrected, sudoseed is set.
const LQSEQS_LABLE_HETE: u8 = 0b0100_0000; //heterozygous marker
const LQSEQS_LABLE_RECH: u8 = 0b0010_0000; //2 or more dominant genotypes, need to do recheck

#[derive(Default, Debug)]
struct LqSeqs {
    lable: u8, // a lable
    start: u32,
    end: u32,
    sudoseed: String,
    seqs: Vec<LqSeq>,
}

impl fmt::Display for LqSeqs {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "start:{} end:{} SUCC:{}|HETE:{}|RECH:{} sudoseed: {} {}",
            self.start,
            self.end,
            self.has_lable(LQSEQS_LABLE_SUCC),
            self.has_lable(LQSEQS_LABLE_HETE),
            self.has_lable(LQSEQS_LABLE_RECH),
            self.sudoseed.len(),
            self.sudoseed
        )?;
        for seq in &self.seqs {
            writeln!(
                f,
                "order:{} kscore:{} seqlen:{} seq:{}",
                seq.order,
                seq.kscore,
                seq.seq.len(),
                seq.seq
            )?;
        }
        Ok(())
    }
}

impl LqSeqs {
    fn set_lable(&mut self, lable: u8) {
        self.lable |= lable;
    }

    #[allow(dead_code)]
    fn unset_lable(&mut self, lable: u8) {
        self.lable ^= lable;
    }

    fn has_lable(&self, lable: u8) -> bool {
        self.lable & lable != 0
    }

    fn clean_seqs(&mut self) {
        self.seqs = Vec::new();
    }

    fn retain_sort_seqs(&mut self, stat: &HashMap<u32, usize>, min_c: usize) {
        self.seqs
            .sort_by_key(|v| Reverse(stat.get(&v.order).unwrap_or(&0))); // this sort maintains the order of self.seqs.order
        let mut c = 0;
        for order in self.seqs.iter().map(|v| v.order) {
            if *stat.get(&order).unwrap_or(&0) < min_c {
                break;
            }
            c += 1;
        }
        self.seqs.truncate(c);
        self.seqs.shrink_to_fit();
    }
}

#[allow(dead_code)]
fn display_lqseqs_vec(lqseqs: &[LqSeqs]) {
    println!("#count of lqseqs:{:?}", lqseqs.len());
    for (p, lqseq) in lqseqs.iter().enumerate() {
        if !lqseq.has_lable(LQSEQS_LABLE_HETE) {
            // continue;
        }
        println!("#{p} {lqseq}\n");
    }
}

fn retrieve_kmer_count(lqseqs: &mut [LqSeqs], kmer_info: &KmerInfo, min_kmer_count: u16) {
    let ksize = kmer_info.ksize as usize;
    let mut kmer_hash: HashMap<u64, u16> = HashMap::default();
    for lqseq in lqseqs.iter() {
        for seq in &lqseq.seqs {
            // for lqseqs with length > ksize
            if seq.seq.len() > ksize {
                for kmer in iter2kmer(seq.seq.chars().map(|x| x as u8), ksize) {
                    kmer_hash.insert(kmer_info.to_dump_key(kmer), 0);
                }
            } else if seq.kmer != INVALID_KMER {
                kmer_hash.insert(seq.kmer, 0);
            }
        }
    }

    kmer_info.retrieve_kmers_count_from_dump(&mut kmer_hash);

    for lqseq in lqseqs.iter_mut() {
        for seq in &mut lqseq.seqs {
            if seq.seq.len() > ksize {
                seq.kscore = iter2kmer(seq.seq.chars().map(|x| x as u8), ksize)
                    .map(|x| {
                        let kmer = kmer_info.to_dump_key(x);
                        kmer_hash
                            .get(&kmer)
                            .map_or(0, |&v| if v <= min_kmer_count { 0 } else { v })
                    })
                    .min()
                    .unwrap_or(0);
            } else if seq.kmer != INVALID_KMER {
                seq.kscore =
                    kmer_hash
                        .get(&seq.kmer)
                        .map_or(0, |&v| if v <= min_kmer_count { 0 } else { v });
            }
        }
    }
}

fn is_valid_snp(seq1: &[u8], seq2: &[u8]) -> bool {
    let mut i = 0;
    let mut j = 0;
    //check whether the SSR compressed lqseqs are the same, becasue HiFi data are prone to errors in the SSR region.
    while i < seq1.len() && j < seq2.len() {
        if seq1[i] != seq2[j] {
            return true;
        }

        while i + 1 < seq1.len() && seq1[i] == seq1[i + 1] {
            i += 1;
        }
        while j + 1 < seq2.len() && seq2[j] == seq2[j + 1] {
            j += 1;
        }
        i += 1;
        j += 1;
    }

    // such as GAGCTCT vs GAGCTCTCT, so here we always return false
    false
}

fn get_min_count(c: usize)-> usize{

    if c >= 9 {
        3
    } else if c >= 6 {
        2
    } else {
        1
    }
}

fn fill_order_stat(lqseq: &LqSeqs, stats: &mut[usize], order_stat: &mut HashMap<u32, usize>) -> (usize, usize, usize, usize){
    let (mut max1_c, mut max1_p, mut max2_c, mut max2_p) = (0, 0, 0, 0);
    stats.fill(0);
    order_stat.clear();

    for (p1, seq) in lqseq.seqs.iter().enumerate().filter(|(_, v)| v.kscore > 0) {
        if stats[p1] > 0 {
            continue;
        }

        let c = lqseq.seqs[p1..].iter().filter(|x| x.seq == seq.seq).count();
        order_stat.insert(lqseq.seqs[p1].order, c);

        for (p2, _) in lqseq.seqs[p1..]
            .iter()
            .enumerate()
            .filter(|(_, v)| v.seq == seq.seq)
        {
            stats[p1 + p2] = c;
        }

        if c > max1_c || (c == max1_c && seq.order == 0) {
            max2_c = max1_c;
            max2_p = max1_p;
            max1_c = c;
            max1_p = p1;
        } else if max1_p == max2_p || c > max2_c {
            max2_c = c;
            max2_p = p1;
        }
    }
    (max1_c, max1_p, max2_c, max2_p)
}

fn fill_seed_lqseqs(lqseqs: &mut [LqSeqs], select_by_count: bool) {
    let mut stats = [0; LQSEQ_MAX_CAN_COUNT];
    let mut order_stat = HashMap::default();

    for lqseq in lqseqs.iter_mut() {
        let (max1_c, max1_p, _, _) = fill_order_stat(lqseq, &mut stats, &mut order_stat);
        lqseq.sudoseed = lqseq.seqs[max1_p].seq.to_owned();
        lqseq.set_lable(LQSEQS_LABLE_SUCC);

        lqseq.set_lable(LQSEQS_LABLE_RECH); //this lqseq need do be re-check
        let min_c = get_min_count(lqseq.seqs.len());
        
        if !select_by_count  {
            assert_eq!(lqseq.seqs[0].order, 0, "the first lqseq is not ref.");
            
            // make sure the lqseq from ref will be saved, here *v/c > 1 is to avoid switch err.
            if let Some(v) = order_stat.get_mut(&0) {
                if *v > 1 && *v < min_c{
                    *v = min_c;
                }
            }else {
                let c = lqseq.seqs.iter().filter(|x| x.seq == lqseq.seqs[0].seq).count();
                if c > 1 {
                    order_stat.insert(0, min_c);
                }
            }

            if max1_p != 0 && max1_c < min_c { // ln case max1_p is the only correct kmer.
                *order_stat.get_mut(&lqseq.seqs[max1_p].order).unwrap() = min_c;
            }
        };

        lqseq.retain_sort_seqs(&order_stat, min_c);
        
        if lqseq.seqs.len() <= 1 {
            if !lqseq.seqs.is_empty() {//IN case the sudoseed/lqseq is not from ref
                lqseq.sudoseed.clear();
                lqseq.sudoseed.push_str(&lqseq.seqs[0].seq);
            }
            lqseq.unset_lable(LQSEQS_LABLE_RECH);
            lqseq.clean_seqs();
        }  
    }
}

fn mark_hete_lqseqs(lqseqs: &mut [LqSeqs]) {
    let mut stats = [0; LQSEQ_MAX_CAN_COUNT];
    let mut order_stat = HashMap::default();

    for lqseq in lqseqs.iter_mut() {
        let (_, max1_p, max2_c, max2_p) = fill_order_stat(lqseq, &mut stats, &mut order_stat);
        let min_c = get_min_count(lqseq.seqs.len());

       if max2_c >= min_c && (lqseq.seqs.len() >= 6 || lqseq.seqs[max1_p].seq.len() == lqseq.seqs[max2_p].seq.len()) 
            && is_valid_snp(lqseq.seqs[max1_p].seq.as_bytes(), lqseq.seqs[max2_p].seq.as_bytes()) {
            lqseq.set_lable(LQSEQS_LABLE_HETE);
            
            for (p, seq) in lqseq.seqs.iter_mut().enumerate().filter(|(_, v)| v.kscore > 0){
                if stats[p] < min_c {
                    seq.kscore = 0;
                }
            }
        }
    }
}

fn phase_reads_by_lqseqs(lqseqs: &[LqSeqs], asref: bool) -> Vec<u32> {
    let mut data = new_data();
    let mut dif = new_data();
    let mut ref_data = new_data();
    for lqseq in lqseqs.iter().filter(|x| x.has_lable(LQSEQS_LABLE_HETE)) {
        for i in 0..lqseq.seqs.len() {
            let seq1 = &lqseq.seqs[i];
            if seq1.kscore == 0 {
                continue;
            }
            for j in i + 1..lqseq.seqs.len() {
                let seq2 = &lqseq.seqs[j];
                if seq2.kscore == 0 {
                    continue;
                }
                let w = if seq1.seq == seq2.seq {
                    // here, may be need to consider SNP quality, not always 1/-1
                    1.
                } else {
                    -1.
                };

                //ref. should not be used for louvain algorithm.
                if seq1.order == 0 {
                    if asref {
                        insert_data(&mut ref_data, seq1.order, seq2.order, w);
                    }
                    continue;
                }
                assert!(seq2.order != 0, "seq2 order is equal to 0");

                if w == -1. {
                    insert_data(&mut dif, seq1.order, seq2.order, -1.);
                    insert_data(&mut dif, seq2.order, seq1.order, -1.);
                }

                insert_data(&mut data, seq1.order, seq2.order, w);
                insert_data(&mut data, seq2.order, seq1.order, w);
            }
        }
    }
    
    // Here we may need to re-consider?
    // If A and B contain >= 3 different SNPs, then A and B are from different phases.
    for (n1, n1_v) in dif.into_iter() {
        for (n2, w) in n1_v.into_iter() {
            if w <= -3. {
                assign_data(&mut data, n1, n2, w);
            }
        }
    }

    //some valid ids are not in valid_ids, so here we should to return invalid_ids
    phase_communities(data, ref_data.into_values().next())
}

fn get_lqseqs_next_idx_by_lable(lqseqs: &[LqSeqs], mut lqseqs_i: usize, lable: u8) -> usize {
    //here use lqseqs_i < lqseqs.len() can make sure lqseqs_i >= 0
    lqseqs_i -= 1;
    while lqseqs_i < lqseqs.len() && !lqseqs[lqseqs_i].has_lable(lable) {
        lqseqs_i -= 1;
    }

    lqseqs_i
}

fn update_consensus_with_lqseqs(
    lqseqs: &[LqSeqs],
    consensus: Vec<ConsensusBase>,
    lable: u8,
) -> Vec<ConsensusBase> {
    let mut consensusbases: Vec<ConsensusBase> = Vec::with_capacity(consensus.len());

    let mut i = 0;
    let mut lqseqs_i = get_lqseqs_next_idx_by_lable(lqseqs, lqseqs.len(), lable);
    while i < consensus.len() {
        let p = consensus[i].pos;
        if lqseqs_i < lqseqs.len() && p == lqseqs[lqseqs_i].start {
            for base in lqseqs[lqseqs_i].sudoseed.chars() {
                consensusbases.push(ConsensusBase {
                    pos: p,
                    // qv: '!',
                    base,
                });
            }

            while i < consensus.len() && consensus[i].pos <= lqseqs[lqseqs_i].end {
                i += 1;
            }

            lqseqs_i = get_lqseqs_next_idx_by_lable(lqseqs, lqseqs_i, lable);
        } else {
            consensusbases.push(consensus[i]);
            i += 1;
        }
    }
    consensusbases
}

fn reupdate_consensus_with_lqseqs(
    lqseqs: &mut [LqSeqs],
    consensus: Vec<ConsensusBase>,
    kmer_info: &KmerInfo,
    min_kmer_count: u16,
    select_by_count: bool
) -> Vec<ConsensusBase> {
    let ksize = kmer_info.ksize as usize;
    let mut kmer_hash: HashMap<u64, u16> = HashMap::default();
    let mut i = 0;
    let mut lqseqs_i = get_lqseqs_next_idx_by_lable(lqseqs, lqseqs.len(), LQSEQS_LABLE_RECH);
    while i < consensus.len() {
        if lqseqs_i < lqseqs.len() && consensus[i].pos == lqseqs[lqseqs_i].start {
            let start = if i > ksize { i - ksize + 1 } else { 0 };
            for seq in lqseqs[lqseqs_i].seqs.iter().map(|x| &x.seq) {
                let mut p = i;
                let iter = consensus[start..p]
                    .iter()
                    .map(|x| x.base)
                    .chain(seq.chars());
                while consensus[p].pos <= lqseqs[lqseqs_i].end {
                    p += 1;
                }
                let mut end = p + ksize - 1;
                if end > consensus.len(){
                    end = consensus.len();
                }
                let iter = iter.chain(consensus[p..end].iter().map(|x| x.base));
                for kmer in iter2kmer(iter.map(|x| x as u8), ksize) {
                    kmer_hash.insert(kmer_info.to_dump_key(kmer), 0);
                }
            }

            while i < consensus.len() && consensus[i].pos <= lqseqs[lqseqs_i].end {
                i += 1;
            }

            lqseqs_i = get_lqseqs_next_idx_by_lable(lqseqs, lqseqs_i, LQSEQS_LABLE_RECH);
        } else {
            i += 1;
        }
    }

    kmer_info.retrieve_kmers_count_from_dump(&mut kmer_hash);

    i = 0;
    lqseqs_i = get_lqseqs_next_idx_by_lable(lqseqs, lqseqs.len(), LQSEQS_LABLE_RECH);
    while i < consensus.len() {
        if lqseqs_i < lqseqs.len() && consensus[i].pos == lqseqs[lqseqs_i].start {
            let start = if i > ksize { i - ksize + 1 } else { 0 };
            let lqseq_end = lqseqs[lqseqs_i].end;
            for seq in lqseqs[lqseqs_i].seqs.iter_mut() {
                let mut p = i;
                let iter = consensus[start..p]
                    .iter()
                    .map(|x| x.base)
                    .chain(seq.seq.chars());
                while consensus[p].pos <= lqseq_end {
                    p += 1;
                }
                let mut end = p + ksize - 1;
                if end > consensus.len(){
                    end = consensus.len();
                }
                let iter = iter.chain(consensus[p..end].iter().map(|x| x.base));
                seq.kscore = iter2kmer(iter.map(|x| x as u8), ksize)
                    .map(|x| {
                        let kmer = kmer_info.to_dump_key(x);
                        kmer_hash
                            .get(&kmer)
                            .map_or(0, |&v| if v <= min_kmer_count { 0 } else { v })
                    })
                    .min()
                    .unwrap_or(0);
            }
            while i < consensus.len() && consensus[i].pos <= lqseqs[lqseqs_i].end {
                i += 1;
            }

            lqseqs_i = get_lqseqs_next_idx_by_lable(lqseqs, lqseqs_i, LQSEQS_LABLE_RECH);
        } else {
            i += 1;
        }
    }

    for lqseq in lqseqs.iter_mut().filter(|x| x.has_lable(LQSEQS_LABLE_RECH)) {
        let mut c = 0;
        let mut valid_count = 0;
        for (p, seq) in lqseq.seqs.iter().enumerate() {
            // lqseq.seqs has been sorted by count
            if seq.kscore != 0  {
                if c == 0 || seq.order == 0 {//in case ref is prefer
                    c = p + 1;
                }

                valid_count += 1;
            }
        }

        //two or more valid kmers, need to recheck using longer kmer databases.
        if valid_count > 1 {
            lqseq.set_lable(LQSEQS_LABLE_TEMP);
        }

        if c != 0 {
            lqseq.sudoseed.clear();
            lqseq.sudoseed.push_str(&lqseq.seqs[c - 1].seq);
        } else if !select_by_count {// keep the reference sequence unchanged if all lqseqs are invalid.
            lqseq.sudoseed.clear();
            let mut i = 0;
            for (p, seq) in lqseq.seqs.iter().enumerate(){
                if seq.order == 0 {
                    i = p;
                    break;
                }
            }
            lqseq.sudoseed.push_str(&lqseq.seqs[i].seq);
        }
    }
    // display_lqseqs_vec(&lqseqs);
    let consensus = update_consensus_with_lqseqs(lqseqs, consensus, LQSEQS_LABLE_RECH);

    //clean LQSEQS_LABLE_TEMP flag
    for lqseq in lqseqs.iter_mut().filter(|x| x.has_lable(LQSEQS_LABLE_RECH)){
        if lqseq.has_lable(LQSEQS_LABLE_TEMP){
            lqseq.unset_lable(LQSEQS_LABLE_TEMP);
        }else {
            lqseq.unset_lable(LQSEQS_LABLE_RECH);
        }
    }

    consensus
}

fn generate_lqseqs_from_tags_kmer(
    alignseqs: &mut [AlignSeq],
    mut lqseqs: Vec<LqSeqs>,
    consensus: Vec<ConsensusBase>,
    opt: &Opt,
    out_cns: bool,
) -> Option<Vec<ConsensusBase>> {
    let mut align_bases: Vec<AlignBase> = Vec::with_capacity(1_000_000);

    let kmer_info = &opt.yak[0];
    let ksize = kmer_info.ksize as u64;
    let shift: u64 = 2 * (ksize - 1);
    let mask: u64 = (1 << (2 * ksize)) - 1;
    let mut kmers: [u64; 2] = [0, 0];
    let mut l;

    let mut j;
    let mut s = lqseqs.len() - 1;
    for (idx, _align_bases) in alignseqs
        .iter()
        .enumerate()
        .filter(|(_, x)| !x.align_bases.is_empty())
    {
        //find start index
        while s > 0 && lqseqs[s].start < _align_bases.aln_t_s {
            s -= 1;
        }
        if lqseqs[s].start < _align_bases.aln_t_s || lqseqs[s].end > _align_bases.aln_t_e {
            continue;
        }

        //find end index
        j = s;
        while j > 0 && lqseqs[j].end <= _align_bases.aln_t_e {
            j -= 1;
        }
        if lqseqs[j].end > _align_bases.aln_t_e {
            j += 1;
        }

        align_bases.clear();
        let mut p = 0;
        let mut align_base: AlignBase = Default::default();
        while _align_bases.get_align_tag(&mut p, &mut align_base) {
            align_bases.push(align_base);
            if align_base.t_pos > lqseqs[j].end + ksize as u32 {
                // out range
                break;
            }
        }

        for lqseq in &mut lqseqs[j..=s] {
            if lqseq.seqs.len() >= LQSEQ_MAX_CAN_COUNT {
                continue;
            }

            l = 0;
            kmers[0] = 0;
            kmers[1] = 0;
            let mut seq = String::with_capacity(10);
            for align_base in &align_bases[lqseq.start as usize - _align_bases.aln_t_s as usize..] {
                if align_base.t_pos >= lqseq.start && align_base.q_base != 4 {
                    if align_base.t_pos <= lqseq.end {
                        seq.push(SEQ_NUM[align_base.q_base as usize] as char);
                    }

                    if l < ksize {
                        kmers[0] = (kmers[0] << 2 | align_base.q_base as u64) & mask;
                        kmers[1] = (kmers[1] >> 2) | (3 ^ align_base.q_base as u64) << shift;
                        l += 1;
                    }

                    if align_base.t_pos > lqseq.end && l >= ksize {
                        break;
                    }
                }
            }

            let kmer = if l >= ksize {
                if kmers[0] < kmers[1] {
                    kmers[0]
                } else {
                    kmers[1]
                }
            } else {
                INVALID_KMER
            };
            if !seq.is_empty() {
                seq.shrink_to_fit();
                lqseq.seqs.push(LqSeq {
                    order: idx as u32,
                    kscore: 0,
                    kmer: if kmer != INVALID_KMER {
                        kmer_info.to_dump_key(kmer)
                    } else {
                        INVALID_KMER
                    },
                    seq,
                });
            }
        }
    }
    
    retrieve_kmer_count(&mut lqseqs, kmer_info, opt.min_kmer_count);
    // display_lqseqs_vec(&lqseqs);
    if out_cns{
        // display_lqseqs_vec(&lqseqs);
        fill_seed_lqseqs(&mut lqseqs, opt.select_by_count);
        let mut consensus = update_consensus_with_lqseqs(&lqseqs, consensus, LQSEQS_LABLE_SUCC);

        for kmer_info in &opt.yak{
            consensus = reupdate_consensus_with_lqseqs(&mut lqseqs, consensus, kmer_info, opt.min_kmer_count, opt.select_by_count);
        }

        Some(consensus)
    }else {
        mark_hete_lqseqs(&mut lqseqs);
        let invalid_ids = phase_reads_by_lqseqs(&lqseqs, opt.model == "ref");
        for id in invalid_ids {
            alignseqs[id as usize].align_bases = Vec::new();
        }
        None
    }
}

fn generate_cns_from_best_score_lq<'a>(
    msas: &'a [Msa],
    alignseqs: &mut [AlignSeq],
    mut global_best_kmer: &'a Kmer,
    opt: &Opt,
    out_cns: bool,
) -> Option<Vec<ConsensusBase>> {
    let mut lqseqs: Vec<LqSeqs> = Vec::with_capacity(100);
    let mut consensusbases: Vec<ConsensusBase> = Vec::with_capacity(msas.len());

    let hq_min_qv = 95;
    let lq_min_length = 2;
    let mut has_lq = false;
    let mut lq_s = usize::MAX;
    let mut lq_e = 0;
    let mut p = 0;

    let (_, mut base2, mut base3) = global_best_kmer.bases(msas.len() as u32 - 1);
    loop {
        if base3.q_base != 4 {
            let coverage = msas[base3.t_pos as usize].coverage();
            let qv = global_best_kmer.count as i64 * 100 / coverage;

            // println!("{:?} {} {} {} {}", base3.t_pos, qv, coverage,  global_best_kmer.to_string(), global_best_kmer.count);

            consensusbases.push(ConsensusBase {
                pos: base3.t_pos,
                // qv: qv as u8 as char,
                base: SEQ_NUM[base3.q_base as usize] as char,
            });

            if coverage < 2 {
                has_lq = false;
                lq_s = usize::MAX;
            } else if qv < hq_min_qv {
                if lq_s == usize::MAX {
                    lq_s = p;
                }
                lq_e = p;
                has_lq = true;
            } else if has_lq
                && p - lq_e > 2 * lq_min_length
                && consensusbases[p - 1].pos != consensusbases[p - 2].pos
                && consensusbases[p - 1].base != consensusbases[p - 2].base
            {
                lq_e = p - 2;
                lq_s = if lq_s > lq_min_length {
                    lq_s - lq_min_length
                } else {
                    1
                };
                while lq_s > 1
                    && (consensusbases[lq_s - 1].pos == consensusbases[lq_s].pos
                        || consensusbases[lq_s - 1].base == consensusbases[lq_s].base)
                {
                    lq_s -= 1;
                }
                let lqseqs_index = lqseqs.len();
                if lqseqs_index >= 1 && consensusbases[lq_s].pos >= lqseqs[lqseqs_index - 1].start {
                    lqseqs[lqseqs_index - 1].start = consensusbases[lq_e].pos;
                } else {
                    lqseqs.push(LqSeqs {
                        end: consensusbases[lq_s].pos,
                        start: consensusbases[lq_e].pos,
                        lable: 0,
                        ..Default::default()
                    });
                }
                has_lq = false;
                lq_s = usize::MAX;
            }
            p += 1;
        }

        if base2.is_head() {
            break;
        }
        global_best_kmer = &msas[base2.t_pos as usize].kmers[global_best_kmer.besti.get() as usize];
        (_, base2, base3) = global_best_kmer.bases(base2.t_pos);
    }
    // drop(msas);

    consensusbases.reverse();
    if lqseqs.is_empty() {
        Some(consensusbases)
    } else {
        generate_lqseqs_from_tags_kmer(alignseqs, lqseqs, consensusbases, opt, out_cns)
    }
}

fn get_cns_from_align_tags(
    msas: &[Msa],
    alignseqs: &mut [AlignSeq],
    opt: &Opt,
    out_cns: bool,
) -> Option<Vec<ConsensusBase>> {
    let mut global_best_kmer: &Kmer = &Default::default();

    for (p, msa) in msas.iter().enumerate() {
        for kmer in &msa.kmers {
            let (base1, base2, _) = kmer.bases(p as u32);
            let coverage = msa.coverage();
            let mut besti = 0;
            let kmer_score = if base2.is_head() {
                10 * kmer.count as i64 - 4 * coverage
            } else {
                let mut kmer_score = i64::MIN >> 1; //prevent overflow
                for (pi, pkmer) in msas[base2.t_pos as usize].get(base1, base2) {
                    let base1 = pkmer.bases(base2.t_pos).0;
                    // this can prevent the later backtracking algorithm from stopping at the start
                    // mapping position of reads, not the start position of the reference.
                    if base2.t_pos >= 3 && base1.is_head() {
                        continue;
                    }
                    let score = pkmer.score.get() + 10 * kmer.count as i64 - 4 * coverage;
                    if score > kmer_score || (score == kmer_score && base1.q_base != 4) {
                        kmer_score = score;
                        besti = pi;
                    }
                }
                kmer_score
            };
            kmer.score.set(kmer_score);
            kmer.besti.set(besti);

            if p == msas.len() - 1 && kmer_score >= global_best_kmer.score.get() {
                global_best_kmer = kmer;
            }
        }
    }

    generate_cns_from_best_score_lq(msas, alignseqs, global_best_kmer, opt, out_cns)
}

fn main() {
    let opt = Opt::from_args();
    let sec_seqs = if opt.use_secondary {
        //secondary alignments do not contain the SEQ field, so retrieve sequences first.
        retrieve_secondary_seq_from_bam(&opt.bam, opt.thread)
    } else {
        HashMap::default()
    };

    thread::scope(|work| {
        let opt = &opt;
        let (in_s, in_r) = bounded(opt.thread + 1);
        let (ou_s, ou_r) = bounded(opt.thread + 1);

        work.spawn(move |_| {
            // input thread
            let mut records = parse_path(&opt.fa).unwrap();
            while let Some(record) = records.iter_record().unwrap() {
                assert!(
                    record.len() < u32::MAX as usize,
                    "{} is too long!",
                    record.head()
                );
                in_s.send((record.head().to_string(), record.seq().to_string()))
                    .unwrap();
            }
        });

        work.spawn(move |_| {
            // work thread
            thread::scope(|scoped| {
                (0..opt.thread).for_each(|_| {
                    let in_r = in_r.clone();
                    let ou_s = ou_s.clone();
                    let sec_seqs = &sec_seqs;
                    scoped.spawn(move |_| {
                        while let Ok((tid, tseq)) = in_r.recv() {
                            let mut aln = Alignment::new();
                            let mut alignseqs: Vec<AlignSeq> = Vec::with_capacity(65536);
                            let mut msas: Vec<Msa> = vec![Default::default(); tseq.len()];
                            aln.aln_t_e = tseq.len() as u32;
                            aln.aln_q_e = tseq.len() as u32;
                            aln.q_aln_str = tseq.to_owned(); //TODO use cow
                            aln.t_aln_str = tseq.to_owned(); //TODO use cow
                            alignseqs.push(AlignSeq::new(&aln));
                            if tseq.len() >= 50_000 {
                                aln.clear();
                                aln.shrink_to(50_000);
                            }

                            let mut bam = IndexedReader::from_path(&opt.bam).unwrap();
                            bam.fetch((&tid, 0, tseq.len() as u32))
                                .expect("Faield random access BAM/SAM!");

                            let mut r = Record::new();
                            let (mut pre_tid, mut pre_pos) = (0, 0);
                            while let Some(ret) = bam.read(&mut r) {
                                ret.expect("BAM/SAM parsing failed!");
                                assert!(
                                    r.tid() > pre_tid || r.reference_start() >= pre_pos,
                                    "Unsorted input file!"
                                );

                                let rlen = r.seq_len_from_cigar(true);
                                if r.flags() & 0x404 != 0
                                    || r.mapq() as i16 <= opt.min_map_qual
                                    || rlen <= opt.min_read_len
                                    || (r.is_secondary() && !opt.use_secondary)
                                    || (r.is_supplementary() && !opt.use_supplementary)
                                    || (r.reference_end() - r.reference_start()
                                        < max(
                                            opt.min_map_len as i64,
                                            (rlen as f32 * opt.min_map_fra) as i64,
                                        ))
                                {
                                    continue;
                                }

                                aln.clear();
                                aln.aln_t_s = r.reference_start() as u32;
                                if opt.use_secondary && r.is_secondary() {
                                    let temp;
                                    let qseq = if r.is_reverse() {
                                        let qseq = &sec_seqs[r.qname()];
                                        temp = reverse_complement_seq_u8(qseq, qseq.len());
                                        &temp
                                    } else {
                                        &sec_seqs[r.qname()]
                                    };
                                    aln.fill_with_cigar(
                                        &r.cigar(),
                                        &tseq[aln.aln_t_s as usize..],
                                        qseq,
                                    );
                                } else {
                                    aln.fill_with_cigar(
                                        &r.cigar(),
                                        &tseq[aln.aln_t_s as usize..],
                                        &r.seq(),
                                    );
                                }
                                let is_clip =
                                    aln.aln_q_e - aln.aln_q_s + opt.max_clip_len < rlen as u32;
                                aln.trim(8);

                                if aln.aln_len() <= opt.min_map_len {
                                    continue;
                                }
                                // println!("{} {}",String::from_utf8(r.qname().to_vec()).unwrap(), alignseqs.len());

                                let mut alignseq = AlignSeq::new(&aln);
                                if is_clip {
                                    if tseq.len() < 500_000{// short references usually contain lots of incorrectly alignments with clipping
                                        continue;
                                    }
                                    alignseq.set_lable(); //should call unset_lable before use the aln_t_s field
                                }
                                alignseqs.push(alignseq);
                                pre_tid = r.tid();
                                pre_pos = r.reference_start();
                            }
                            filter_alignseqs_by_clip(&mut alignseqs);

                            let mut i = 0;
                            let cns = loop {
                                update_msas(&mut msas, &alignseqs);
                                shrink_msas(&mut msas);
                                sort_msas(&mut msas);
                                if i + 1 == opt.iter_count {
                                    break get_cns_from_align_tags(
                                        &msas,
                                        &mut alignseqs,
                                        opt,
                                        true,
                                    );
                                } else {
                                    get_cns_from_align_tags(&msas, &mut alignseqs, opt, false);
                                    clear_msas(&mut msas);
                                }
                                i += 1;
                            };
                            ou_s.send((tid, cns)).unwrap();
                        }
                    });
                });
            })
            .unwrap();
        });

        work.spawn(move |_| {
            // output thread
            while let Ok((tid, cns)) = ou_r.recv() {
                display_consensusbase_vec(&tid, &cns.unwrap(), opt.uppercase, opt.out_pos);
                //safely unwrap
            }
        });
    })
    .unwrap();

    eprintln!("{}", resource_str());
}
