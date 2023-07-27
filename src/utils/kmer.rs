//most of function were adopted from yak:https://github.com/lh3/yak
use fxhash::FxHashSet as HashSet;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::BufReader;
use std::io::Read;

const YAK_MAGIC: &[u8] = "YAK\u{2}".as_bytes();
const YAK_COUNTER_BITS: u32 = 10;

pub const SEQ_NUM: [u8; 128] = [
    // translate ACGTU-NM to 01233456
    //A, C,  G,  T,   -,  N, M
    65, 67, 71, 84, 45, 78, 77, 4, 4, 4, 4, 4, 4, 4, 4, 4, //0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //48-63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 6, 5, 4, //64-79
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //80-95
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 6, 5, 4, //96-111
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //112-127
];



#[derive(Clone, Debug)]
struct Kmer {
    hash: u64,
}

impl Hash for Kmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.kmer().hash(state);
    }
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.kmer() == other.kmer()
    }
}

impl Eq for Kmer {}

impl Kmer {
    fn new(hash: u64) -> Kmer {
        Kmer {
            hash
        }
    }

    fn kmer(&self) -> u64 {
        self.hash >> YAK_COUNTER_BITS
    }

    fn count(&self) -> u16 {
        (self.hash & ((1 << YAK_COUNTER_BITS) - 1)) as u16
    }
}

#[derive(Clone, Debug)]
pub struct KmerInfo {
    pub ksize: u32,
    kmask: u64,
    pre: u32,
    pmask: u64,
    path: String,
    sets: Vec<HashSet<Kmer>>
}

impl KmerInfo {
    pub fn new(path: &str) -> Self {
        let mut f = File::open(path).unwrap();
        let mut buffer = [0; 16];
        f.read_exact(&mut buffer).unwrap();
        assert_eq!(
            &buffer[..4],
            YAK_MAGIC,
            "The input binary k-mer dump file is incompatible."
        );
        let (ksize, pre, counter_bits) = unsafe {
            let (n1, n2, n3) = buffer[4..].align_to::<u32>();
            assert!(
                n1.is_empty() && n3.is_empty() && n2.len() == 3,
                "Failed to parse the dump file header."
            );
            (n2[0], n2[1], n2[2])
        };

        assert_eq!(counter_bits, YAK_COUNTER_BITS, "different YAK_COUNTER_BITS");

        KmerInfo {
            ksize,
            kmask: (1 << (2 * ksize as u64)) - 1,
            pre,
            pmask: (1 << pre as u64) - 1,
            path: path.to_string(),
            sets: vec![HashSet::default(); 1 << pre],
        }
    }

    pub fn to_hash(&self, kmer: u64) -> u64 {
        if self.ksize < 32 {
            //Here we don't hash kmers first, so it easier to convert kmers into seqs for debugging.
            yak_hash64(kmer, self.kmask)
        } else {
            //for ksize >= 32, kmer has been hashed
            kmer
        }
    }

    //clear_count: set count as zero
    pub fn insert(&mut self, hash: u64, clear_count: bool) -> bool {
        self.sets[(hash & self.pmask) as usize].insert(Kmer::new(
            if clear_count{
                hash >> YAK_COUNTER_BITS << YAK_COUNTER_BITS
            }else {
                hash
            }
        ))
    }

    pub fn get(&self, hash: u64) -> Option<u16> {
        self.sets[(hash & self.pmask) as usize].get(&Kmer::new(hash)).map(|x| x.count())
    }

    //return true if successfully replaced
    pub fn replace(&mut self, hash: u64) -> bool {
        self.sets[(hash & self.pmask) as usize].replace(Kmer::new(hash)).is_some()
    }

    pub fn retrieve_kmers(&mut self, min_count: u16) -> Vec<u32> {
        let f = File::open(&self.path).unwrap();
        let mut reader = BufReader::new(f);
        let mut buffer = [0; 16];
        reader.read_exact(&mut buffer).unwrap();
        let max_count: u64 = (1 << YAK_COUNTER_BITS) - 1;

        let mut hist = vec![0; max_count as usize + 1];
        for i in 0..1 << self.pre {
            let mut buffer = [0; 8];
            reader.read_exact(&mut buffer).unwrap();
            let size = u32::from_ne_bytes(
                buffer[4..]
                    .try_into()
                    .expect("Failed to parse size from the dump file header."),
            );
            let sets = &mut self.sets[i];
            let mut buffer = [0u8; std::mem::size_of::<u64>()];
            for _j in 0..size {
                let res = reader.read_exact(&mut buffer);
                match res {
                    Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => break,
                    _ => {}
                }
                res.expect("Failed to parse the dump file");
                let hash = u64::from_ne_bytes(buffer);
                let count = (hash & max_count) as u16;
                hist[count as usize] += 1;
                if count < min_count{
                    continue;
                }
                let kmer = Kmer::new(hash);
                if sets.contains(&kmer){
                    sets.replace(kmer);
                }
            }
        }
        hist
    }

    pub fn load_kmers(&mut self, min_count: u16) -> Vec<u32> {
        let f = File::open(&self.path).unwrap();
        let mut reader = BufReader::new(f);
        let mut buffer = [0; 16];
        reader.read_exact(&mut buffer).unwrap();
        let max_count: u64 = (1 << YAK_COUNTER_BITS) - 1;

        let mut hist = vec![0; max_count as usize + 1];
        for i in 0..1 << self.pre {
            let mut buffer = [0; 8];
            reader.read_exact(&mut buffer).unwrap();
            let size = u32::from_ne_bytes(
                buffer[4..]
                    .try_into()
                    .expect("Failed to parse size from the dump file header."),
            );
            let sets = &mut self.sets[i];
            sets.reserve(size as usize);
            let mut buffer = [0u8; std::mem::size_of::<u64>()];
            for _j in 0..size {
                let res = reader.read_exact(&mut buffer);
                match res {
                    Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => break,
                    _ => {}
                }
                res.expect("Failed to parse the dump file");
                let hash = u64::from_ne_bytes(buffer);
                let count = (hash & max_count) as u16;
                hist[count as usize] += 1;
                if count < min_count {
                    continue;
                }
                let t = sets.insert(Kmer::new(hash));
                assert_eq!(t,  true);
            }
        }
        hist
    }

    //estimated kmer total count
    pub fn estimated_len(&self) -> usize {
        std::fs::metadata(&self.path).expect("failed to read read metadata").len() as usize / std::mem::size_of::<u64>()
    }

    pub fn clear(&mut self) {
        for set in &mut self.sets{
            set.clear();
        }
    }
}

fn yak_hash64(key: u64, mask: u64) -> u64 // invertible integer hash function
{
    let mut key = (!key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    key
}

fn yak_hash64_64(key: u64) -> u64 {
    let mut key = !key + (key << 21);
    key = key ^ key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key = key ^ key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key = key ^ key >> 28;
    key = key + (key << 31);
    key
}

fn yak_hash_long(x: [u64; 4]) -> u64 {
    let j = if x[1] < x[3] { 0 } else { 1 };
    yak_hash64_64(x[j << 1]) + yak_hash64_64(x[j << 1 | 1])
}

pub fn seq2kmer(seq: &str, ksize: usize) -> impl Iterator<Item = u64> + '_ {
    iter2kmer(seq.chars().map(|x| x as u8), ksize)
}

pub fn iter2kmer<'a>(
    mut iter: impl Iterator<Item = u8> + 'a,
    ksize: usize,
) -> Box<dyn Iterator<Item = u64> + 'a> {
    let mut l = 0;
    if ksize < 32 {
        let shift: u64 = 2 * (ksize as u64 - 1);
        let mask: u64 = (1 << (2 * ksize as u64)) - 1;
        let mut kmer: [u64; 2] = [0, 0];
        Box::new(std::iter::from_fn(move || {
            loop {
                if let Some(c) = iter.next() {
                    let c = SEQ_NUM[c as usize] as u64;
                    if c < 4 {
                        kmer[0] = (kmer[0] << 2 | c) & mask; // forward k-mer
                        kmer[1] = (kmer[1] >> 2) | (3 ^ c) << shift; // reverse k-mer
                        l += 1;
                    } else {
                        l = 0;
                    }

                    if l >= ksize {
                        return if kmer[0] < kmer[1] {
                            Some(kmer[0])
                        } else {
                            Some(kmer[1])
                        }; // strand
                    }
                } else {
                    return None;
                }
            }
        }))
    } else {
        let shift: u64 = ksize as u64 - 1;
        let mask: u64 = (1 << (ksize as u64)) - 1;
        let mut kmer: [u64; 4] = [0, 0, 0, 0];
        Box::new(std::iter::from_fn(move || loop {
            if let Some(c) = iter.next() {
                let c = SEQ_NUM[c as usize] as u64;
                if c < 4 {
                    kmer[0] = (kmer[0] << 1 | (c & 1)) & mask;
                    kmer[1] = (kmer[1] << 1 | (c >> 1)) & mask;
                    kmer[2] = kmer[2] >> 1 | (1 - (c & 1)) << shift;
                    kmer[3] = kmer[3] >> 1 | (1 - (c >> 1)) << shift;
                    l += 1;
                } else {
                    l = 0;
                    kmer.fill(0);
                }

                if l >= ksize {
                    return Some(yak_hash_long(kmer));
                }
            } else {
                return None;
            }
        }))
    }
}

pub fn kmer2seq(kmer: u64, ksize: usize) -> String {
    let mut seq = String::new();
    for i in (0..ksize).rev() {
        let c = kmer >> (i * 2) & 3;
        seq.push(SEQ_NUM[c as usize] as char)
    }
    seq
}
