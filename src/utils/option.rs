use clap::{builder::PossibleValue, value_parser, Arg, ArgAction, ArgMatches, Command};
use libc_stdhandle::stdout;
use path_absolutize::*;
use std::{
    ffi::CString,
    io::{Error, ErrorKind, Result},
    path::Path,
};

use crate::KmerInfo;

const VERSION: &str = include_str!(concat!(env!("OUT_DIR"), "/VERSION"));

#[derive(Clone, Debug)]
pub struct Option {
    pub yak: Vec<KmerInfo>,
    pub bam: String,
    pub fa: String,
    pub model: String,           //-m
    pub uppercase: bool,         //-u
    pub out_pos: bool,           
    pub min_kmer_count: u16,     //-k
    pub thread: usize,           //-t
    pub iter_count: usize,       //-i
    pub min_read_len: usize,     //-l
    pub use_supplementary: bool, //-s
    pub use_secondary: bool,     //-S
    pub use_all_reads: bool,     //-r
    pub min_map_qual: i16,       //-q // user may set this as -1
    pub min_map_len: usize,      //-a
    pub min_map_fra: f32,        //-a
    pub max_clip_len: u32,       //-c
    pub min_base_cov: usize,
}

impl Option {
    pub fn new() -> Option {
        Option::default()
    }

    pub fn from_args() -> Option {
        let opt = Option::default();
        let args = Command::new("nextPolish2")
            .version(VERSION)
            .about("repeat-aware polishing genome assemblies using HiFi and short reads")
            .arg_required_else_help(true)
            .arg(
                Arg::new("bam")
                    .value_name("HiFi.map.bam")
                    .value_parser(to_abspath_string)
                    .required(true)
                    .help("HiFi-to-ref mapping file in sorted BAM format."),
            )
            .arg(
                Arg::new("fa")
                    .value_name("genome.fa[.gz]")
                    .value_parser(to_abspath_string)
                    .required(true)
                    .help("genome assembly file in [GZIP] FASTA format."),
            )
            .arg(
                Arg::new("yak")
                    .value_name("short.read.yak")
                    .value_parser(to_abspath_string)
                    .required(true)
                    .action(ArgAction::Append)
                    .help("one or more k-mer dataset in yak format."),
            )
            .arg(
                Arg::new("out")
                    .short('o')
                    .long("out")
                    .value_name("FILE")
                    .default_value("stdout")
                    .value_parser(|x: &str| {
                        if x != "stdout" {
                            freopen_stdout(x)
                        } else {
                            Ok(())
                        }
                    })
                    .help("output file."),
            )
            .arg(
                Arg::new("uppercase")
                    .short('u')
                    .long("uppercase")
                    .help("output in uppercase sequences.")
                    .action(ArgAction::SetTrue),
            )
            .arg(
                Arg::new("out_pos")
                    .long("out_pos")
                    .help("output each base and its position.")
                    .hide_short_help(true)
                    .action(ArgAction::SetTrue),
            )
            .arg(
                Arg::new("min_kmer_count")
                    .short('k')
                    .long("min_kmer_count")
                    .value_name("INT")
                    .default_value(opt.min_kmer_count.to_string())
                    .value_parser(value_parser!(u16))
                    .help("filter kmers in k-mer dataset with count <= INT."),
            )
            .arg(
                Arg::new("thread")
                    .short('t')
                    .long("thread")
                    .value_name("INT")
                    .default_value(opt.thread.to_string())
                    .value_parser(value_parser!(usize))
                    .help("number of threads."),
            )
            .arg(
                Arg::new("iter_count")
                    .short('i')
                    .long("iter_count")
                    .value_name("INT")
                    .default_value(opt.iter_count.to_string())
                    .value_parser(value_parser!(usize))
                    .help("number of iterations to attempt phasing."),
            )
            .arg(
                Arg::new("model")
                    .short('m')
                    .long("model")
                    .value_name("STR")
                    .ignore_case(true)
                    .default_value(&opt.model)
                    .help("phasing model.")
                    .value_parser([
                        PossibleValue::new("ref").help("output the same haplotype phase blocks as the reference"),
                        PossibleValue::new("len").help("output longer haplotype phase blocks"),
                    ]),
            )
            .arg(
                Arg::new("min_read_len")
                    .short('l')
                    .long("min_read_len")
                    .value_name("INT")
                    .default_value(opt.min_read_len.to_string())
                    .value_parser(value_parser!(usize))
                    .help("filter reads with length <= INT."),
            )
            .arg(
                Arg::new("use_supplementary")
                    .short('s')
                    .long("use_supplementary")
                    .help("use supplementary alignments.")
                    .action(ArgAction::SetTrue),
            )
            .arg(
                Arg::new("use_secondary")
                    .short('S')
                    .long("use_secondary")
                    .help("use secondary alignments, consider setting `min_map_qual` to -1 when using this option.")
                    .action(ArgAction::SetTrue),
            )
            .arg(
                Arg::new("min_map_len")
                    .long("min_map_len")
                    .short('a')
                    .value_name("INT.FLOAT")
                    .default_value(opt.min_map_fra.to_string())
                    .value_parser(value_parser!(f32))
                    .help("filter alignments with alignment length <= min(INT, FLOAT * read_length)."),
            )
            .arg(
                Arg::new("min_map_qual")
                    .long("min_map_qual")
                    .short('q')
                    .value_name("INT")
                    .default_value(opt.min_map_qual.to_string())
                    .value_parser(value_parser!(i16))
                    .allow_negative_numbers(true)
                    .help("filter alignments with mapping quality <= INT."),
            )
            .arg(
                Arg::new("max_clip_len")
                    .long("max_clip_len")
                    .short('c')
                    .value_name("INT")
                    .default_value(opt.max_clip_len.to_string())
                    .value_parser(value_parser!(u32))
                    .help("filter alignments with unaligned length >= INT."),
            )
            .arg(
                Arg::new("use_all_reads")
                    .short('r')
                    .long("use_all_reads")
                    .help("use all unfiltered reads, reads with different haplotypes from the reference assembly are discarded by default.")
                    .action(ArgAction::SetTrue),
            )
            .arg(
                Arg::new("min_base_cov")
                    .long("min_base_cov")
                    .value_name("INT")
                    .default_value(opt.min_base_cov.to_string())
                    .value_parser(value_parser!(usize))
                    .help("minimum depth to correct a raw base.")
                    .hide(true),
            )
            .get_matches();

        opt.update(args)
    }

    fn update(self, mut args: ArgMatches) -> Option {
        //safely unwrap, becasue the default values have been set
        let min_map_len = args.remove_one::<f32>("min_map_len").unwrap();
        let mut yaks: Vec<KmerInfo> = args.remove_many::<String>("yak").expect("Missing yak file!").map(|x| KmerInfo::new(&x)).collect();
        yaks.sort_by_key(|k| k.ksize);

        Option {
            bam: args.remove_one::<String>("bam").expect("Missing bam file!"),
            fa: args
                .remove_one::<String>("fa")
                .expect("Missing fasta file!"),
            yak: yaks,
            model: args.remove_one::<String>("model").unwrap(),
            uppercase: args.get_flag("uppercase"),
            out_pos: args.get_flag("out_pos"),
            use_all_reads: args.get_flag("use_all_reads"),
            min_kmer_count: args.remove_one::<u16>("min_kmer_count").unwrap(),
            thread: args.remove_one::<usize>("thread").unwrap(),
            iter_count: args.remove_one::<usize>("iter_count").unwrap(),
            min_read_len: args.remove_one::<usize>("min_read_len").unwrap(),
            use_supplementary: args.get_flag("use_supplementary"),
            use_secondary: args.get_flag("use_secondary"),
            min_map_len: min_map_len as usize,
            min_map_fra: min_map_len.fract(),
            min_map_qual: args.remove_one::<i16>("min_map_qual").unwrap(),
            max_clip_len: args.remove_one::<u32>("max_clip_len").unwrap(),
            min_base_cov: args.remove_one::<usize>("min_base_cov").unwrap(),
        }
    }
}

impl Default for Option {
    fn default() -> Self {
        Option {
            yak: Vec::new(),
            bam: String::new(),
            fa: String::new(),
            model: "ref".to_string(),
            uppercase: false,
            out_pos: false,
            min_kmer_count: 1,
            thread: 1,
            iter_count: 2,
            min_read_len: 1000,
            use_supplementary: false,
            use_secondary: false,
            use_all_reads: false,
            min_map_len: 500,
            min_map_fra: 500.5,
            min_map_qual: 1,
            max_clip_len: 100,
            min_base_cov: 1,
        }
    }
}

fn to_abspath_string(path: &str) -> Result<String> {
    let path = Path::new(path)
        .absolutize()
        .expect("Failed convert input file to abspath!");
    if path.exists() {
        Ok(path.to_string_lossy().to_string())
    } else {
        Err(Error::new(
            ErrorKind::NotFound,
            format!("{:?} does not exist!", path),
        ))
    }
}

fn freopen_stdout(path: &str) -> Result<()> {
    let path = Path::new(path)
        .absolutize()
        .expect("Failed convert input file to abspath!");
    if path.exists() {
        return Err(Error::new(
            ErrorKind::AlreadyExists,
            format!("{:?} already exists!", path),
        ));
    } else {
        let w = CString::new("w").unwrap();
        let p = CString::new(path.to_string_lossy().as_bytes()).unwrap();
        if unsafe { libc::freopen(p.as_ptr(), w.as_ptr(), stdout()) }.is_null() {
            return Err(Error::new(
                ErrorKind::Other,
                format!("Failed to freopen: {:?}", path),
            ));
        }
    }
    Ok(())
}
