use std::env;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

// Copied from https://blog.biofan.org/2019/08/cargo-build-script/
fn get_git_version() -> String {
    let version = env::var("CARGO_PKG_VERSION").unwrap();

    let child = Command::new("git").args(["describe", "--always"]).output();
    match child {
        Ok(child) => {
            let buf = String::from_utf8(child.stdout).expect("failed to read stdout");
            version + "-" + &buf
        }
        Err(err) => {
            eprintln!("`git describe` err: {}", err);
            version
        }
    }
}

fn build_yak() {
    let c_dir = PathBuf::from("./yak");
    for file in (c_dir.read_dir().expect("failed to read dir")).flatten() {
        let fpath = file.path();
        if fpath
            .extension()
            .map_or_else(|| true, |x| x != "o" && x != "a")
        {
            println!("cargo:rerun-if-changed={:?}", fpath);
        }
    }

    let output = Command::new("make")
        .arg("all")
        .current_dir(&c_dir)
        .output()
        .expect("failed to execute process");
    if !output.status.success() {
        panic!("make error: {}", String::from_utf8_lossy(&output.stderr));
    }
}

fn main() {
    build_yak();
    let version = get_git_version();
    let mut f = File::create(Path::new(&env::var("OUT_DIR").unwrap()).join("VERSION")).unwrap();
    f.write_all(version.trim().as_bytes()).unwrap();
}
