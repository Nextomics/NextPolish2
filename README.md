# NextPolish2
Haplotype-aware polishing genome assemblies using HiFi and short reads

## Table of Contents

- [Installation](#install)
- [General usage](#usage)
- [Algorithm overview](#algorithm)
- [Getting help](#help)
- [Citation](#cite)
- [License](#license)
- [Limitations](#limit)
- [Benchmarking](#benchmark)

### <a name="install"></a>Installation

#### Dependencies

`NextPolish2` is written in rust, try below commands (no root required) or refer [here](https://www.rust-lang.org/tools/install) to install `Rust` first.
```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

#### Download and install

```sh
git clone --recursive git@github.com:Nextomics/NextPolish2.git
cd NextPolish2 && cargo build --release
```

#### Test

```sh
cd test && bash hh.sh
```

### <a name="usage"></a>General usage

NextPolish2 takes a genome assembly file, a HiFi mapping file and one or more k-mer database files as input and generates the polished genome.

1. Prepare HiFi mapping file

```sh
# mapping
minimap2 -ax map-hifi -t 5 asm.fa.gz hifi.fasta.gz|samtools sort -o hifi.map.sort.bam
# indexing
samtools index hifi.map.sort.bam
```

2. Prepare k-mer database files. Here we only produce 21-mer and 51-mer databases, you can produce more k-mer databases with different k-mer size

```sh
# produce a 21-mer database
./yak/yak count -o k21.yak -k 21 <(zcat sr.R*.fastq.gz) <(zcat sr.R*.fastq.gz)
# produce a 51-mer database
./yak/yak count -o k51.yak -k 51 <(zcat sr.R*.fastq.gz) <(zcat sr.R*.fastq.gz) 
```

3. Run NextPolish2

```sh
./target/release/nextPolish2 -t 5 hifi.map.sort.bam asm.fa.gz k21.yak k51.yak > asm.np2.fa
```

#### More options

Use `./target/release/nextPolish2 -h` to see options.

### <a name="algorithm"></a>Algorithm overview

Coming soon...

### <a name="help"></a>Getting help

#### Help

   Feel free to raise an issue at the [issue page](https://github.com/Nextomics/NextPolish2/issues/new).

   ***Note:*** Please ask questions on the issue page first. They are also helpful to other users.
#### Contact
   
   For additional help, please send an email to huj\_at\_grandomics\_dot\_com.

### <a name="cite"></a>Citation

Coming soon...

### <a name="license"></a>License

NextPolish2 is only freely available for academic use and other non-commercial use.

### <a name="limit"></a>Limitations

Coming soon...

### <a name="benchmark"></a>Benchmarking

Coming soon...

### Star
You can track updates by tab the **Star** button on the upper-right corner at the [github page](https://github.com/Nextomics/NextPolish2).