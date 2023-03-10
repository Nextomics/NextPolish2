# Benchmarking NextPolish2

## Requirement
- [pIRS](https://github.com/galaxy001/pirs) v2.0.0
- [pbsim3](https://github.com/yukiteruono/pbsim3) v3.0.0
- [samtools](https://github.com/samtools/samtools) v1.9
- [ccs](https://ccs.how/) v6.4.0
- [art_illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) v2.5.8
- [hifiasm](https://github.com/chhylp123/hifiasm) v0.18.5
- [yak](https://github.com/lh3/yak) v0.1
- [minimap2](https://github.com/lh3/minimap2) v2.24
- [nextPolish2](https://github.com/Nextomics/NextPolish2) v0.1.0
- [meryl](https://github.com/marbl/meryl) v1.4
- [merqury](https://github.com/marbl/merqury) v1.3

## Software commands

1. Download reference
```sh
wget https://www.arabidopsis.org/download_files/Genes/Col-CEN%20genome%20assembly%20release/ColCEN.fasta
```

2. Simulate diploid genome
```sh
pirs diploid --snp-rate 0.005 --indel-rate 0.002 --random-seed=1 ColCEN.fasta -o alta.simu
cat ColCEN.fasta alta.simu.snp.indel.inversion.fa > alta.sim.diploid.fa
```

3. Simulate PacBio HiFi data
```sh
pbsim --strategy wgs --method errhmm --errhmm ERRHMM-SEQUEL.model --depth 30 --genome alta.sim.diploid.fa --length-mean 13000 --prefix alta.sim --id-prefix alta.sim --length-sd 2000 --pass-num 10
ls *sam|while read line;do samtools view -b $line > $line.bam && conda run ccs $line.bam $line.hifi.bam && samtools fasta $line.hifi.bam > $line.fa; done;
cat alta.sim.hifi*.sam.fa > alta.sim.diploid.hifi.fa
```

4. Simulate Illumina data
```sh
art_illumina -ss HS25 -i ColCEN.fasta -p -l 150 -f 50 -m 300 -s 10 -o alta.sim.hap1.R
#output: alta.sim.hap1.R1.fq alta.sim.hap1.R2.fq 

art_illumina -ss HS25 -i alta.simu.snp.indel.inversion.fa -p -l 150 -f 50 -m 300 -s 10 -o alta.sim.hap2.R
#output: alta.sim.hap2.R1.fq alta.sim.hap2.R2.fq 
```

5. Assemble reference
```sh
hifiasm -o alta.sim.diploid.hifiasm -t 20 alta.sim.diploid.hifi.fa
awk '{if ($1~/^S/ && length($3)>=1000000){print ">"$2;print $3}}' alta.sim.diploid.hifiasm.p_ctg.gfa > alta.sim.diploid.hifiasm.p_ctg.gfa.fa
```

6. Polishing with NextPolish2
```sh
thread=5
yak count -t ${thread} -k 21 -o sr.k21.yak <(cat alta.sim.hap*.fq)
yak count -t ${thread} -k 51 -o sr.k51.yak <(cat alta.sim.hap*.fq)
minimap2 -ax map-hifi -t ${thread}  alta.sim.diploid.hifiasm.p_ctg.gfa.fa alta.sim.diploid.hifi.fa |samtools sort --threads ${thread} -o alta.sim.diploid.hifi.fa.bam -
samtools index alta.sim.diploid.hifi.fa.bam
nextPolish2 -t ${thread} -k 0 alta.sim.diploid.hifi.fa.bam alta.sim.diploid.hifiasm.p_ctg.gfa.fa sr.k21.yak sr.k51.yak -o alta.sim.diploid.hifiasm.p_ctg.gfa.polish.fa
```

7. Phasing assessment with hap-mers
```sh
meryl k=19 count alta.sim.hap1.R*.fq output alta.sim.hap1.meryl
meryl k=19 count alta.sim.hap2.R*.fq output alta.sim.hap2.meryl
meryl k=19 count alta.sim.hap*.fq output alta.sim.hap12.meryl

asm=alta.sim.diploid.hifiasm.p_ctg.gfa.fa
mkdir ${asm}.merqury
$MERQURY/merqury.sh alta.sim.hap12.meryl alta.sim.hap1.meryl alta.sim.hap2.meryl ${asm} ${asm}.merqury/out

asm=alta.sim.diploid.hifiasm.p_ctg.gfa.polish.fa
mkdir ${asm}.merqury
$MERQURY/merqury.sh alta.sim.hap12.meryl alta.sim.hap1.meryl alta.sim.hap2.meryl ${asm} ${asm}.merqury/out

```