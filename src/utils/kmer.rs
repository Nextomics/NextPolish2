//most of function were adopted from yak:https://github.com/lh3/yak
use fxhash::FxHashMap as HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;

const YAK_MAGIC: &[u8] = "YAK\u{2}".as_bytes();

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
pub struct KmerInfo {
    pub ksize: u32,
    l_pre: u32,
    counter_bits: u32,
    path: String,
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
        let (ksize, l_pre, counter_bits) = unsafe {
            let (n1, n2, n3) = buffer[4..].align_to::<u32>();
            assert!(
                n1.is_empty() && n3.is_empty() && n2.len() == 3,
                "Failed to parse the dump file header."
            );
            (n2[0], n2[1], n2[2])
        };

        KmerInfo {
            ksize,
            l_pre,
            counter_bits,
            path: path.to_string(),
        }
    }

    pub fn to_dump_key(&self, kmer: u64) -> u64 {
        if self.ksize < 32 {
            let mask: u64 = (1 << (2 * self.ksize as u64)) - 1;
            yak_hash64(kmer, mask) >> self.counter_bits
        }else { //for ksize >= 32, kmer has been hashed
            kmer >> self.counter_bits
        }
    }

    pub fn retrieve_kmers_count_from_dump(&self, map: &mut HashMap<u64, u16>) -> Vec<u32> {
        let f = File::open(&self.path).unwrap();
        let mut reader = BufReader::new(f);
        let mut buffer = [0; 16];
        reader.read_exact(&mut buffer).unwrap();
        let max_count: u64 = (1 << self.counter_bits) - 1;

        let mut hist = vec![0; max_count as usize + 1];
        for _i in 0..1 << self.l_pre {
            let mut buffer = [0; 8];
            reader.read_exact(&mut buffer).unwrap();
            // let capacity = u32::from_ne_bytes(buffer[..4].try_into().expect("Failed to parse capacity from the dump file header."));
            let size = u32::from_ne_bytes(
                buffer[4..]
                    .try_into()
                    .expect("Failed to parse size from the dump file header."),
            );
            let mut buffer = [0u8; std::mem::size_of::<u64>()];
            for _j in 0..size {
                let res = reader.read_exact(&mut buffer);
                match res {
                    Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => break,
                    _ => {}
                }
                res.expect("Failed to parse the dump file");
                let kmer_count = u64::from_ne_bytes(buffer);
                let (kmer, count) = (kmer_count >> self.counter_bits, kmer_count & max_count);
                hist[count as usize] += 1;
                if let Some(x) = map.get_mut(&kmer) {
                    //some kmer=>hash collisions here
                    let count = count as u16;
                    if *x < count {
                        *x = count;
                    }
                }
                // println!("{:?} {:?}", kmer, count);
            }
        }
        hist
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
    let j = if x[1] < x[3] {0} else {1};
    return yak_hash64_64(x[j << 1 | 0]) + yak_hash64_64(x[ j << 1 | 1]);
}

pub fn seq2kmer(seq: &str, ksize: usize) -> impl Iterator<Item = u64> + '_ {
    iter2kmer(seq.chars().map(|x| x as u8), ksize)
}

pub fn iter2kmer<'a>(mut iter: impl Iterator<Item = u8> + 'a, ksize: usize) -> Box<dyn Iterator<Item = u64> + 'a> {
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
    }else {
        let shift: u64 = ksize as u64 - 1;
        let mask: u64 = (1 << (ksize as u64)) - 1;
        let mut kmer: [u64; 4] = [0, 0, 0, 0];
        Box::new(std::iter::from_fn(move || {
            loop {
                if let Some(c) = iter.next() {
                    let c = SEQ_NUM[c as usize] as u64;
                    if c < 4 {
                        kmer[0] = (kmer[0] << 1 | (c&1))  & mask;
                        kmer[1] = (kmer[1] << 1 | (c>>1)) & mask;
                        kmer[2] = kmer[2] >> 1 | (1 - (c&1))  << shift;
                        kmer[3] = kmer[3] >> 1 | (1 - (c>>1)) << shift;
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
