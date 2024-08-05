# Benchmarking NextPolish2 overcorrection using data from [Homo sapiens (HG005)](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/)

## Requirement
- [fastp](https://github.com/OpenGene/fastp) v0.23.4
- [hifiasm](https://github.com/chhylp123/hifiasm) v0.18.5
- [yak](https://github.com/lh3/yak) v0.1
- [meryl](https://github.com/marbl/meryl) v1.4
- [winnowmap](https://github.com/marbl/Winnowmap) v2.03
- [minimap2](https://github.com/lh3/minimap2) v2.26
- [paftools.js](https://github.com/lh3/minimap2/tree/master/misc) v2.26
- [nextPolish2](https://github.com/Nextomics/NextPolish2) v0.2.1
- [deepvariant](https://github.com/google/deepvariant) v1.6.1
- [merqury](https://github.com/marbl/merqury) v1.3

## Software commands

1. Reads download and filtering
```sh
#Download HiFi data
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/8ced218e-9699-458e-9708-eef6969a8065--EXTRAMURAL_SAMPLES/HG005/PacBio_HiFi/PBmixSequel788_1_A01_PBXX_30hours_19kbV2PD_70pM_HumanHG005_CCS/m64109_200304_195708.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/8ced218e-9699-458e-9708-eef6969a8065--EXTRAMURAL_SAMPLES/HG005/PacBio_HiFi/PBmixSequel789_2_B01_PBXX_30hours_19kbV2PD_70pM_HumanHG005_CCS/m64109_200311_013444.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/8ced218e-9699-458e-9708-eef6969a8065--EXTRAMURAL_SAMPLES/HG005/PacBio_HiFi/PBmixSequel789_1_A01_PBXW_30hours_15kbV2PD_70pM_HumanHG005_CCS/m64109_200309_192110.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/8ced218e-9699-458e-9708-eef6969a8065--EXTRAMURAL_SAMPLES/HG005/PacBio_HiFi/PBmixSequel840_1_A01_PCCD_30hours_15kbV2PD_70pM_HumanHG005_CCS/m64017_200723_190224.fastq.gz

# HiFi data statistics:
# [length stat]
# Types Count (#) Length (bp)
# N10    425798      21313   
# N20    891131      20026   
# N30    1382540     19067   
# N40    1896930     18264   
# N50    2433397     17514   
# N60    2993357     16768   
# N70    3577790     16081   
# N80    4186930     15431   
# N90    4822776     14722   

# Min.      -         45   
# Max.      -        31877  
# Ave.      -        17358
# Total  5527654  95951017712

# Download NGS data
wget ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_250bps_fastqs/150506_HG005_Homogeneity_04_FCB-22365345/*/*.fastq.gz

# Quality control and filtering.
fastp -t 32 -f 5 -t 5 --cut_front --cut_tail -i R1.fastq.gz -I R2.fastq.gz -o R1.clean.fastq.gz -O R2.clean.fastq.gz

# NGS data statistics after quality control filtering:
# [length stat]
# Types Count (#) Length (bp)
# N10   78414841      234    
# N20   156829681     234    
# N30   235244521     234    
# N40   313659361     234    
# N50   392074202     234    
# N60   470801088     233    
# N70   549552473     233    
# N80   628511063     232    
# N90   708163275     227    

# Min.      -         15     
# Max.      -         234    
# Ave.      -         229    
# Total 800733216 183490726099
```
***Note: As detailed in the NextPolish2 documentation, the performance of NextPolish2 relies heavily on the quality of short reads. We strongly recommend performing quality control and filtering on NGS data.***

2. Assemble reference
```sh
# primary assembly
hifiasm -o hg005 -t 48 m64017_200723_190224.fastq.gz m64109_200304_195708.fastq.gz m64109_200309_192110.fastq.gz m64109_200311_013444.fastq.gz
awk '/^S/{print ">"$2;print $3}' hg005.bp.p_ctg.gfa > hg005.bp.p_ctg.gfa.fa
```

3. Polishing with `NextPolish2`
```sh
# mapping using winnowmap
meryl count k=15 output merylDB hg005.bp.p_ctg.gfa.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
winnowmap -t 5 -W repetitive_k15.txt -ax map-pb hg005.bp.p_ctg.gfa.fa m64017_200723_190224.fastq.gz m64109_200304_195708.fastq.gz m64109_200309_192110.fastq.gz m64109_200311_013444.fastq.gz|samtools sort -o hifi.map.sort.bam -

# indexing
samtools index hifi.map.sort.bam

# Prepare k-mer dataset files
yak count -t 5 -k 21 -b 37 -o sr.k21.yak <(zcat *clean.fastq.gz) <(zcat *clean.fastq.gz)
yak count -t 5 -k 31 -b 37 -o sr.k31.yak <(zcat *clean.fastq.gz) <(zcat *clean.fastq.gz)

#Run NextPolish2
nextPolish2 -t 5 hifi.map.sort.bam hg005.bp.p_ctg.gfa.fa sr.k21.yak sr.k31.yak -o hg005.bp.p_ctg.gfa.polish.fa
```

4. Performance assessment using `merqury`
```sh
meryl k=21 count  *clean.fastq.gz output sr.k21.meryl
sh $MERQURY/merqury.sh sr.k21.meryl hg005.bp.p_ctg.gfa.fa hg005.bp.p_ctg.gfa.fa.out
sh $MERQURY/merqury.sh sr.k21.meryl hg005.bp.p_ctg.gfa.polish.fa hg005.bp.p_ctg.gfa.polish.fa.out

# Here is the QV result:
# hg005.bp.p_ctg.gfa  267547  3094873399  53.8544 4.11676e-06
# hg005.bp.p_ctg.gfa.polish   177940  3094828306  55.6257 2.73798e-06
```

5. Performance assessment using `deepvariant`
```sh
# hg005.bp.p_ctg.gfa.fa
minimap2 -ax map-hifi -t 60 hg005.bp.p_ctg.gfa.fa m64017_200723_190224.fastq.gz m64109_200304_195708.fastq.gz m64109_200309_192110.fastq.gz m64109_200311_013444.fastq.gz |samtools sort -o hg005.bp.p_ctg.gfa.fa.bam -
samtools index hg005.bp.p_ctg.gfa.fa.bam
docker run --privileged -u $(id -u):$(id -g) -v `pwd`:/input -v `pwd`:/output google/deepvariant:"1.6.1" /deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=hg005.bp.p_ctg.gfa.fa --reads=hg005.bp.p_ctg.gfa.fa.bam --output_vcf=/output/hg005.bp.p_ctg.gfa.fa.vcf  --output_gvcf=/output/hg005.bp.p_ctg.gfa.fa.gvcf --num_shards=32 --logging_dir=/output/hg005.bp.p_ctg.gfa.fa.logs --intermediate_results_dir=/output/hg005.bp.p_ctg.gfa.fa.tmp --dry_run=false

# hg005.bp.p_ctg.gfa.polish.fa
minimap2 -ax map-hifi -t 60 hg005.bp.p_ctg.gfa.polish.fa m64017_200723_190224.fastq.gz m64109_200304_195708.fastq.gz m64109_200309_192110.fastq.gz m64109_200311_013444.fastq.gz |samtools sort -o hg005.bp.p_ctg.gfa.polish.fa.bam -
samtools index hg005.bp.p_ctg.gfa.polish.fa.bam
docker run --privileged -u $(id -u):$(id -g) -v `pwd`:/input -v `pwd`:/output google/deepvariant:"1.6.1" /deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=hg005.bp.p_ctg.gfa.polish.fa --reads=hg005.bp.p_ctg.gfa.polish.fa.bam --output_vcf=/output/hg005.bp.p_ctg.gfa.polish.fa.vcf  --output_gvcf=/output/hg005.bp.p_ctg.gfa.polish.fa.gvcf --num_shards=32 --logging_dir=/output/hhg005.bp.p_ctg.gfa.polish.fa.logs --intermediate_results_dir=/output/hg005.bp.p_ctg.gfa.polish.fa.tmp --dry_run=false

# Here we count homozygous high-quality variants as potential errors:
grep -v '#' hg005.bp.p_ctg.gfa.fa.vcf|awk '$6>=20' |awk '$7=="PASS"' |grep "1/1"|wc -l
grep -v '#' hg005.bp.p_ctg.gfa.polish.fa.vcf|awk '$6>=20' |awk '$7=="PASS"' |grep "1/1"|wc -l

# Here is the results of homozygous high-quality variant counts:
# hg005.bp.p_ctg.gfa.fa.vcf  14155
# hg005.bp.p_ctg.gfa.polish.fa.vcf  6955 
```

6. Performance assessment using `paftools.js`
```sh
# download the high-confidence variants of the Genome in a Bottle (GIAB) sample HG005 (GRCh37 as ref)
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/latest/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz && gzip -d HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz

# mapping to GRCh37
minimap2 -cx asm10 --cs human_g1k_v37.fasta hg005.bp.p_ctg.gfa.fa -o hg005.bp.p_ctg.gfa.fa.hg19.paf
minimap2 -cx asm10 --cs human_g1k_v37.fasta hg005.bp.p_ctg.gfa.polish.fa -o hg005.bp.p_ctg.gfa.polish.fa.hg19.paf

# variants calling
sort -k6,6 -k8,8n hg005.bp.p_ctg.gfa.fa.hg19.paf| paftools.js call -L 1000000 - > hg005.bp.p_ctg.gfa.fa.hg19.paf.vcf
sort -k6,6 -k8,8n hg005.bp.p_ctg.gfa.polish.fa.hg19.paf| paftools.js call -L 1000000 - > hg005.bp.p_ctg.gfa.polish.fa.hg19.paf.vcf

# Here, we count variants that are not in the high-confidence variants (HG005_GRCh37_1_22_v4.2.1_benchmark.vcf) as potential errors. 
python3 e_use_vcf.py HG005_GRCh37_1_22_v4.2.1_benchmark.vcf hg005.bp.p_ctg.gfa.fa.hg19.paf.vcf
python3 e_use_vcf.py HG005_GRCh37_1_22_v4.2.1_benchmark.vcf hg005.bp.p_ctg.gfa.polish.fa.hg19.paf.vcf

# Here is the results:
# hg005.bp.p_ctg.gfa.fa.hg19.paf.vcf:
# Total matched variants: 2520470, Total potential errors: 771270, potential error rate (%): 0.029927

# hg005.bp.p_ctg.gfa.polish.fa.hg19.paf.vcf:
# Total matched variants: 2522817, Total potential errors: 762908, potential error rate (%): 0.029575

```

The source of e_use_vcf.py:
```python
#!/usr/bin/env python3

import sys

def read_vcf(infile):
    snps = {}
    with open(infile) as IN:
        for line in IN:
            if line.startswith("#"):
                continue
            lines = line.strip().split()
            n, p = lines[0], int(lines[1])
            if n not in snps:
                snps[n] = set()
            snps[n].add(p)
    return snps

hq_vcf_file = sys.argv[1]
eva_pafcall_file = sys.argv[2]

hq_snps = read_vcf(hq_vcf_file)

total_potential_error_var = 0
total_correct_var = 0
total_len = 0
with open(eva_pafcall_file) as IN:
    for line in IN:
        if line.startswith("#"):
            continue
        lines = line.strip().split()
        l, n, s, e = lines[0], lines[1], int(lines[2]), int(lines[3])
        if n not in hq_snps:
            continue
        if l == "R":
            total_len += e - s 
            continue
        elif l == "V":
            if e not in hq_snps[n]:
                total_potential_error_var += 1
            else:
                total_correct_var += 1
print ("Total matched variants: %d, Total potential errors: %d, potential error rate (%%): %f" % (\
    total_correct_var, total_potential_error_var, total_potential_error_var/total_len * 100))
```

