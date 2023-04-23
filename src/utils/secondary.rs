use crossbeam_channel::bounded;
use crossbeam_utils::thread;
use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use rust_htslib::bam::{self, IndexedReader, Read, Record};
use std::str::FromStr;

fn retrieve_secondary_ids(path: &str, thread: usize) -> HashSet<Vec<u8>> {
    //TODO only retrieve secondary ids from opt.genome
    thread::scope(|work| {
        let (in_s, in_r) = bounded(thread + 1);

        work.spawn(move |_| {
            // input thread
            let bam = bam::Reader::from_path(path).expect("Faield random access BAM/SAM!");
            let header = bam::Header::from_template(bam.header());
            for (_, records) in header.to_hashmap().iter().filter(|(k, _)| *k == "SQ") {
                for record in records {
                    in_s.send((
                        record["SN"].to_owned(),
                        u32::from_str(&record["LN"]).unwrap(),
                    ))
                    .unwrap();
                }
            }
        });

        let ret = work.spawn(move |_| {
            // work thread
            thread::scope(|scoped| {
                let mut handles = Vec::with_capacity(thread);
                (0..thread).for_each(|_| {
                    let mut ids = HashSet::default();
                    let in_r = in_r.clone();
                    let handle = scoped.spawn(move |_| {
                        while let Ok((tid, tlen)) = in_r.recv() {
                            let mut bam = IndexedReader::from_path(path).unwrap();
                            bam.fetch((&tid, 0, tlen))
                                .expect("Faield random access BAM/SAM!");

                            let mut r = Record::new();
                            while let Some(ret) = bam.read(&mut r) {
                                ret.expect("BAM/SAM parsing failed!");
                                if r.is_secondary() {
                                    ids.insert(r.qname().to_vec());
                                }
                            }
                        }
                        ids
                    });
                    handles.push(handle);
                });
                let mut ids = HashSet::default();
                for handle in handles {
                    ids.extend(handle.join().unwrap());
                }
                ids
            })
            .unwrap()
        });
        ret.join().unwrap()
    })
    .unwrap()
}

pub fn reverse_complement_seq_u8<T>(seq: &T, len: usize) -> Vec<u8>
where
    T: std::ops::Index<usize, Output = u8>,
{
    (0..len)
        .rev()
        .map(|i| match seq[i] {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            n => n,
        })
        .collect()
}

pub fn retrieve_secondary_seq_from_bam(path: &str, thread: usize) -> HashMap<Vec<u8>, Vec<u8>> {
    let secondary_ids = retrieve_secondary_ids(path, thread);

    thread::scope(|work| {
        let (in_s, in_r) = bounded(thread + 1);

        work.spawn(move |_| {
            // input thread
            let bam = bam::Reader::from_path(path).expect("Faield random access BAM/SAM!");
            let header = bam::Header::from_template(bam.header());
            for (_, records) in header.to_hashmap().iter().filter(|(k, _)| *k == "SQ") {
                for record in records {
                    in_s.send((
                        record["SN"].to_owned(),
                        u32::from_str(&record["LN"]).unwrap(),
                    ))
                    .unwrap();
                }
            }
        });

        let ret = work.spawn(move |_| {
            // work thread
            thread::scope(|scoped| {
                let mut handles = Vec::with_capacity(thread);
                (0..thread).for_each(|_| {
                    let secondary_ids = &secondary_ids;
                    let mut seqs = HashMap::default();
                    let in_r = in_r.clone();
                    let handle = scoped.spawn(move |_| {
                        while let Ok((tid, tlen)) = in_r.recv() {
                            let mut bam = IndexedReader::from_path(path).unwrap();
                            bam.fetch((&tid, 0, tlen))
                                .expect("Faield random access BAM/SAM!");

                            let mut r = Record::new();
                            while let Some(ret) = bam.read(&mut r) {
                                ret.expect("BAM/SAM parsing failed!");
                                if secondary_ids.contains(r.qname())
                                    && !(r.is_secondary() || r.is_supplementary())
                                {
                                    let seq = if r.is_reverse() {
                                        let seq = r.seq();
                                        reverse_complement_seq_u8(&seq, seq.len())
                                    } else {
                                        r.seq().as_bytes()
                                    };
                                    assert!(seqs.insert(r.qname().to_vec(), seq).is_none());
                                }
                            }
                        }
                        seqs
                    });
                    handles.push(handle);
                });
                let mut seqs = HashMap::default();
                for handle in handles {
                    seqs.extend(handle.join().unwrap());
                }
                seqs
            })
            .unwrap()
        });
        ret.join().unwrap()
    })
    .unwrap()
}

// fn main() {
// 	let bam = "hifi.fasta.bam";
// 	let seqs = retrieve_secondary_seq_from_bam(bam, 10);
// 	for (k, v) in seqs{
// 		println!(">{}\n{}", String::from_utf8(k).unwrap(), String::from_utf8(v).unwrap());
// 	}
// }
