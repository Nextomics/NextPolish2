# Benchmarking NextPolish2 using [data from Arabidopsis thaliana (Col-XJTU)](https://ngdc.cncb.ac.cn/gsa/browse/CRA004538)

## Requirement
- [fxTools](https://github.com/moold/fxTools) v0.1.0
- [hifiasm](https://github.com/chhylp123/hifiasm) v0.18.5
- [yak](https://github.com/lh3/yak) v0.1
- [nextPolish2](https://github.com/Nextomics/NextPolish2) v0.1.0
- [meryl](https://github.com/marbl/meryl) v1.4
- [merqury](https://github.com/marbl/merqury) v1.3
- [Racon + Merfin](https://github.com/arangrhie/T2T-Polish/tree/master/automated_polishing)

## Software commands

1. Download reads
```sh
#HiFi
wget https://download.cncb.ac.cn/gsa/CRA004538/CRR302668/CRR302668.fastq.gz

# ngs
wget https://download.cncb.ac.cn/gsa/CRA004538/CRR302670/CRR302670_f1.fastq.gz
wget https://download.cncb.ac.cn/gsa/CRA004538/CRR302670/CRR302670_r2.fastq.gz

```

2. Downsample HiFi reads
```sh
fxTools -s 4.68g CRR302668.fastq.gz > CRR302668.subsample.35x.fastq
```

3. Assemble reference
```sh
hifiasm -o CRR302668.asm -t 30 CRR302668.subsample.35x.fastq
awk '{if ($1~/^S/ && length($3)>=1000000){print ">"$2;print $3}}' CRR302668.asm.bp.p_ctg.gfa > CRR302668.asm.bp.p_ctg.gfa.fa
```

4. Polishing with Racon + Merfin
```sh
k=19
thread=5
# Construct k-mer db
meryl count k=${k} threads=${thread} CRR302670*.fastq.gz output reads.meryl

# Collect histogram for GenomeScope
meryl histogram reads.meryl > reads.hist

# Exclude frequency = 1 k-mers
meryl greater-than 1 reads.meryl output reads.gt1.meryl

bash automated-polishing.sh ${thread} 1 CRR302668.asm.bp.p_ctg.gfa.fasta CRR302668.subsample.35x.fastq reads.gt1.meryl racon.meryl

```

5. Polishing with NextPolish2
```sh
thread=5
yak count -t ${thread} -k 21 -b 37 -o sr.k21.yak <(zcat CRR302670*.fastq.gz) <(zcat CRR302670*.fastq.gz)
yak count -t ${thread} -k 31 -b 37 -o sr.k31.yak <(zcat CRR302670*.fastq.gz) <(zcat CRR302670*.fastq.gz)

HiFi_mapping_file=racon.meryl.iter_1.winnowmap.bam # output by `Racon + Merfin`
samtools index ${HiFi_mapping_file}

nextPolish2 -t ${thread} -r ${HiFi_mapping_file} CRR302668.asm.bp.p_ctg.gfa.fa sr.k21.yak sr.k31.yak -o CRR302668.asm.bp.p_ctg.gfa.np2.fa
```

6. Performance assessment
```sh
genome=CRR302668.asm.bp.p_ctg.gfa.fa #CRR302668.asm.bp.p_ctg.gfa.np2.fa, racon.meryl.iter_1.consensus.fasta
sh $MERQURY/merqury.sh reads.meryl ${genome} out
```