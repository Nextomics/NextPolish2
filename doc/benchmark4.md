# Benchmarking NextPolish2 using data from [Homo sapiens (CHM13)](https://github.com/marbl/CHM13/blob/master/Sequencing_data.md)

## Requirement
- [hifiasm](https://github.com/chhylp123/hifiasm) v0.18.5
- [yak](https://github.com/lh3/yak) v0.1
- [nextPolish2](https://github.com/Nextomics/NextPolish2) v0.1.0
- [meryl](https://github.com/marbl/meryl) v1.4
- [merqury](https://github.com/marbl/merqury) v1.3
- [Racon + Merfin](https://github.com/arangrhie/T2T-Polish/tree/master/automated_polishing)

## Software commands

1. Download reads
```sh
#HiFi data from chm13_4cell
SRR11292120,SRR11292121,SRR11292122,SRR11292123

#NGS data from chm13_4cell
SRR1997411,SRR3189741,SRR3189742,SRR3189743
```

2. Assemble reference
```sh
#primary assembly
hifiasm -o chm13_4cell -t 32 chm13_4cell.hifi.fasta.gz
awk '{if ($1~/^S/ && length($3)>=1000000){print ">"$2;print $3}}' chm13_4cell.p_ctg.gfa > chm13_4cell.p_ctg.gfa.fa
```

3. Polishing with Racon + Merfin
```sh
k=21
thread=5
# Construct k-mer db
meryl count k=${k} threads=${thread} chm13_4cell.ngs.fastq.gz output reads.meryl

# Collect histogram for GenomeScope
meryl histogram reads.meryl > reads.hist

# Exclude frequency = 1 k-mers
meryl greater-than 1 reads.meryl output reads.gt1.meryl

#polish chm13_4cell.p_ctg.gfa.fa
bash automated-polishing.sh ${thread} 1 chm13_4cell.p_ctg.gfa.fa chm13_4cell.hifi.fasta.gz reads.gt1.meryl racon.meryl

```

4. Polishing with NextPolish2
```sh
thread=5
yak count -t ${thread} -k 21 -b 37 -o sr.k21.yak chm13_4cell.ngs.fastq.gz
yak count -t ${thread} -k 31 -b 37 -o sr.k31.yak chm13_4cell.ngs.fastq.gz

HiFi_mapping_file=racon.meryl.iter_1.winnowmap.bam # output by `Racon + Merfin`
samtools index ${HiFi_mapping_file}

genome_file=chm13_4cell.p_ctg.gfa.fa
nextPolish2 -t ${thread} ${HiFi_mapping_file} ${genome_file} sr.k21.yak sr.k31.yak -o ${genome_file}.np2.fa

```

5. Performance assessment
```sh
genome_file=chm13_4cell.p_ctg.gfa.fa  # chm13_4cell.p_ctg.gfa.fa.np2.fa,  racon.meryl.iter_1.consensus.fasta
sh $MERQURY/merqury.sh reads.meryl ${genome_file} out
```