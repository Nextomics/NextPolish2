use ctor::ctor;
use lazy_static::lazy_static;
use libc::{getrusage, rusage, RUSAGE_SELF};
use std::{env, fmt::Write, fs, mem::MaybeUninit, time::Instant};

lazy_static! {
    static ref NOW: Instant = Instant::now();
}

#[ctor]
fn init_NOW() {
    lazy_static::initialize(&NOW);
}

pub fn cputime() -> i64 {
    let r = unsafe {
        let mut r = MaybeUninit::<rusage>::uninit();
        getrusage(RUSAGE_SELF, r.as_mut_ptr());
        r.assume_init()
    };

    r.ru_utime.tv_sec + r.ru_stime.tv_sec
}

pub fn realtime() -> u64 {
    NOW.elapsed().as_secs()
}

pub fn peakrss() -> i64 {
    let r = unsafe {
        let mut r = MaybeUninit::uninit();
        getrusage(RUSAGE_SELF, r.as_mut_ptr());
        r.assume_init()
    };
    r.ru_maxrss
}

pub fn resource_str() -> String {
    let mut s = String::with_capacity(1024);
    let version_file = concat!(env!("OUT_DIR"), "/VERSION");
    if let Ok(v) = fs::read_to_string(version_file) {
        writeln!(&mut s, "Version: {}", v).unwrap();
    }
    s.push_str("CMD:");
    for arg in env::args() {
        write!(&mut s, " {}", arg).unwrap();
    }
    writeln!(
        &mut s,
        "\nReal time: {} sec; CPU: {} sec; Peak RSS: {:.3} GB",
        realtime(),
        cputime(),
        peakrss() as f64 / 1024.0 / 1024.0
    )
    .unwrap();
    s
}
