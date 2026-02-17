# Metagenomic Scripts
Custom perl scripts used to calculate diversity indices, functional similarities and other related stuff
 
## Filter out human contamination from metagenomic reads

1. First, check your Illumina reads using **FASTQC**, and determine your read length and how many low quality bases are both on the front and tail of the reads.

2. Run **fastp**, we will assumme that both the front and tail have 5 low quality bases, and that our reads are 150 bp, using 16 threads. 

```
fastp --in1 [R1] --in2 [R2] --out1 [trimmed_R1] --out2 [trimmed_R2] -z 9 -V -f 5 -F 5 -D --dup_calc_accuracy 6 -l 140 -c -w 16
```

3. Map the trimmed reads to the GRCh38 human genome reference using **bwa-mem2**

```
bwa-mem2 mem -o [mapping.sam] -t 16 [GRCh38 reference] [trimmed_R1] [trimmed_R2]
```

4. Filter out the human contamination using **exclude_human_illumina.pl**

```
perl /media/databases/exclude_human_illumina.pl [mapping.sam]
```
5. Your reads are now ready for further analysis and/or assembly. Delete [mapping.sam]
