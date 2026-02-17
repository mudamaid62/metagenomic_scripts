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

## Assemble your reads using PLASS and filter out partial proteins

1. Assuming your reads are already trimmed and filtered out of human contamination, run **PLASS**

```
plass assemble [R1] [R2] [assembly.faa] [temp_dir] --remove-tmp-files
```

2. Filter out partial proteins using **plass_orf_filter.pl**

```
perl plass_orf_filter.pl [assembly.faa] [PREFIX] [minimun protein length] [maximum protein length]
```

3. After the script is done, you will have 3 output in your current directory. [PREFIX]_complete_orfs.faa is the one appropiate for further usage.
	- [PREFIX]_complete_orfs.faa
	- [PREFIX]_semi-partial_orfs.faa
	- [PREFIX]_partial_orfs.faa

## Assemble your reads using metaSPAdes and filter out short contigs

1. Run **metaSPAdes** using your reads

```
metaspades.py -o [metaspades_assembly] -1 [R1] -2 [R2] -t 16 -m [RAM in Gb]
```

2. Filter out short contigs using **filter_by_length.pl**

```
perl filter_by_length.pl metaspades_assembly/contigs.fasta [minimun contig length] [maximun contig length] > [filtered_assembly.fa]
```

## Run SingleM to obtain community profiles, and then, use ESCGs to calculate alpha-diversity metrics

1. Run **SingleM** in pipe mode to obtain the community profile and OTU table

```
singlem pipe -1 [R1] -2 [R2] -p [taxonomic_profile] --otu-table [OTU_table] --threads 16
```

2. Summarise the profile output in tables by taxonomic levels, to then plot community stack plots. You can summarise multiple profile outputs at once.

```
singlem summarise --input-taxonomic-profiles [taxonomic_profiles] --output-species-by-site-relative-abundance-prefix [PREFIX]
```

3. Summarise the OTU table output by ESCG, to then calculate alpha-diversity metrics. You can summarise multiple OTU tables at once.

```
singlem summarise --input-otu-tables [OTU_tables] --cluster --unifrac-by-otu [PREFIX]
```

4. The previous command generates 59 summarised OTU tables using the PREFIX, one for each ESCG. Then, run **get_indices_from_unifrac_tables.pl** for each ESCG

```
for i in PREFIX*; do perl get_indices_from_unifrac_tables.pl "$i" > "$i"_alpha_diversity_indices; done
```

