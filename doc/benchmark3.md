# Benchmarking NextPolish2 using data from [Homo sapiens (HG002)](https://github.com/genome-in-a-bottle/giab_data_indexes)

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
#HiFi data from hg002
wget https://sra-pub-src-2.s3.amazonaws.com/SRR10382245/m64011_190830_220126.fastq.1
wget https://sra-pub-src-2.s3.amazonaws.com/SRR10382244/m64011_190901_095311.fastq.1
wget https://sra-pub-src-2.s3.amazonaws.com/SRR10382249/m64012_190920_173625.fastq.1
wget https://sra-pub-src-2.s3.amazonaws.com/SRR10382248/m64012_190921_234837.fastq.1

#NGS data from hg002
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz

#NGS data from hg003
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG003/HG003_HiSeq30x_subsampled_R1.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG003/HG003_HiSeq30x_subsampled_R2.fastq.gz

#NGS data from hg004
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG004/HG004_HiSeq30x_subsampled_R1.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/parents/ILMN/downsampled/HG004/HG004_HiSeq30x_subsampled_R2.fastq.gz
```

3. Assemble reference and Trio Binning
```sh
#primary assembly
hifiasm -o hg002 hg002.hifi.fasta.gz
awk '{if ($1~/^S/ && length($3)>=1000000){print ">"$2;print $3}}' hg002.p_ctg.gfa > hg002.p_ctg.gfa.fa

#haplotype-resolved assembly
yak count -k31 -b37 -t16 -o mat.yak hg004.ngs.fastq.gz
yak count -k31 -b37 -t16 -o pat.yak hg003.ngs.fastq.gz
hifiasm -o hg002.dip hg002.hifi.fasta.gz -1 mat.yak -2 pat.yak
awk '{if ($1~/^S/ && length($3)>=1000000){print ">"$2;print $3}}' hg002.dip.hap1.p_ctg.gfa > hg002.dip.hap1.p_ctg.gfa.fa
awk '{if ($1~/^S/ && length($3)>=1000000){print ">"$2;print $3}}' hg002.dip.hap2.p_ctg.gfa > hg002.dip.hap2.p_ctg.gfa.fa

#classify and bin HiFi reads
yak triobin pat.yak mat.yak hg002.hifi.fasta.gz > triobin.out
#extract HiFi reads from hg003/paternal  
awk '$2=="p"||$2=="a"||$2==0' yak.triobin.out |cut -f 1 > hg003.read.list
fxTools getseq -r hg003.read.list hg002.hifi.fasta.gz > hg003.hifi.fasta

#extract HiFi reads from hg004/maternal
awk '$2=="m"||$2=="a"||$2==0' yak.triobin.out |cut -f 1 > hg004.read.list
fxTools getseq -r hg004.read.list hg002.hifi.fasta.gz > hg004.hifi.fasta
```

4. Polishing with Racon + Merfin
```sh
k=21
thread=5
# Construct k-mer db
meryl count k=${k} threads=${thread} hg002.ngs.fastq.gz output reads.meryl

# Collect histogram for GenomeScope
meryl histogram reads.meryl > reads.hist

# Exclude frequency = 1 k-mers
meryl greater-than 1 reads.meryl output reads.gt1.meryl

#polish hg002.p_ctg.gfa.fa
bash automated-polishing.sh ${thread} 1 hg002.p_ctg.gfa.fa hg002.hifi.fasta.gz reads.gt1.meryl racon.meryl

#polish hg002.dip.hap1.p_ctg.gfa.fa
bash automated-polishing.sh ${thread} 1 hg002.dip.hap1.p_ctg.gfa.fa hg004.hifi.fasta reads.gt1.meryl racon.meryl

#polish hg002.dip.hap2.p_ctg.gfa.fa
bash automated-polishing.sh ${thread} 1 hg002.dip.hap2.p_ctg.gfa.fa hg003.hifi.fasta reads.gt1.meryl racon.meryl
```

5. Polishing with NextPolish2
```sh
thread=5
yak count -t ${thread} -k 21 -b 37 -o sr.k21.yak hg002.ngs.fastq.gz
yak count -t ${thread} -k 31 -b 37 -o sr.k31.yak hg002.ngs.fastq.gz

HiFi_mapping_file=racon.meryl.iter_1.winnowmap.bam # output by `Racon + Merfin`
samtools index ${HiFi_mapping_file}

genome_file=hg002.p_ctg.gfa.fa #hg002.dip.hap1.p_ctg.gfa.fa, hg002.dip.hap2.p_ctg.gfa.fa
nextPolish2 -t ${thread} -r ${HiFi_mapping_file} ${genome_file} sr.k21.yak sr.k31.yak -o ${genome_file}.np2.fa

```

6. Performance assessment
```sh
wget https://obj.umiacs.umd.edu/marbl_publications/merqury/HG002/HG002.k21.meryl.tar.gz
wget https://obj.umiacs.umd.edu/marbl_publications/merqury/HG002/HG002.mat.hapmers.meryl.tar.gz
wget https://obj.umiacs.umd.edu/marbl_publications/merqury/HG002/HG002.pat.hapmers.meryl.tar.gz
tar -vxzf HG002.k21.meryl.tar.gz
tar -vxzf HG002.mat.hapmers.meryl.tar.gz
tar -vxzf HG002.pat.hapmers.meryl.tar.gz

genome_file=hg002.p_ctg.gfa.fa #hg002.dip.hap1.p_ctg.gfa.fa, hg002.dip.hap2.p_ctg.gfa.fa, hg002.dip.hap1.p_ctg.gfa.fa.np2.fa, hg002.dip.hap1.p_ctg.gfa.fa.np2.fa, hg002.dip.hap2.p_ctg.gfa.fa.np2.fa,  racon.meryl.iter_1.consensus.fasta
sh $MERQURY/merqury.sh HG002.k21.meryl/ HG002.mat.hapmers.meryl/ HG002.pat.hapmers.meryl/ ${genome_file} out
```