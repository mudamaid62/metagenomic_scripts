# Metagenomic Scripts
Custom perl scripts used to calculate diversity indices, functional similarities and other related stuff
 
### Filter out human contamination from metagenomic reads

1. First, check your Illumina reads using **FASTQC**, determine your read length and how many low quality bases are both on the front and tail of the reads.

2. Run **fastp**, we will assumme that both the front and tail have 5 low quality bases, and that our reads are 150 bp, using 16 threads. 

```
fastp --in1 [R1] --in2 [R2] --out1 [trimmed_R1] --out2 [trimmed_R2] -z 9 -V -f 5 -F 5 -D --dup_calc_accuracy 6 -l 140 -c -w 16
```

3. Map the trimmed reads to the GRCh38 human genome reference using **bwa-mem2**

```
bwa-mem2 mem -o [mapping.sam] -t 16 [GRCh38 reference] [trimmed_R1] [trimmed_R2]
```

4. Filter out the human contamination using **exclude_human_illumina.pl** and **SAMTools**

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

3. After the script is done, you will have 3 output files in your current directory. [PREFIX]_complete_orfs.faa is the one appropiate for further usage.
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

2. Summarise the profile output in tables by taxonomic level, to then plot community stack plots. You can summarise multiple profile outputs at once.

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

## Cluster a protein dataset, recover the clusters containing relevant proteins, and then, retrieve the proteins in those clusters

1. Start with a protein multifasta file of relevant proteins to search, e.g. the BLDB database. Get the protein names

```
grep ">" BLDB.faa > BLDB_list
sed -i 's/>//g' BLDB_list
```

2. Cluster your relevant proteins and some problem proteins, e.g. PLASS assembled metagenomic proteins, using **MMSeqs2**

```
mmseqs createdb BLDB.faa PLASS.faa all_DB
mmseqs linclust all_DB all_DB_clu tmp --alignment-mode 3 --min-seq-id [Minimum identity] -c [Minimun bidirectional coverage] --cov-mode 0 --cluster-mode 2 --threads 16 --realign --remove-tmp-files
mmseqs createsubdb all_DB_clu all_DB all_DB_clu_rep
mmseqs createtsv all_DB all_DB all_DB_clu all_DB_clu.tsv --full-header
```

3. Recover relevant clusters using **get_relevant_clusters_list_v2.pl**

```
perl get_relevant_clusters_list_v2.pl all_DB_clu.tsv BLDB_list > relevant_clusters
````

4. Recover the proteins contained in the relevant clusters using **get_esm_fasta.pl** and **MMSeqs2**

```
cut -f 2 relevant_clusters > relevant_proteins
perl get_esm_fasta.pl all_DB relevant_proteins temp_dir > relevant_proteins.faa
```

## Quantify genes/proteins in copies/cell and get related taxonomies

1. Run **CAT** with metaSPAdes assembled contigs to get contig taxonomic classifications

```
CAT_pack contigs -c metaspades_contigs.fa -d cat_database/db/ -t cat_database/tax/ --no_stars -n 16 --sensitive --block_size 2 --tmpdir temp -o contigs_CAT
```

2. Use **ARGs_OAP** to build a database from a proteins-to-quantify multifasta 

```
args_oap make_db -i proteins.faa
```

3. Run **ARGs_OAP** to get reads that map to the proteins to quantify. 

- Read files must be in a single directory and be called sample_name_1.fastq.gz and sample_name_2.fastq.gz
- Create a structure file for your proteins (tab-separated), where the first column has the header Protein and all protein names, and any ammount of extra columns (at least one) where cluster types are specicified
> e.g.

| Protein | Cluster |
| ----------- | ----------- |
| prot_A | cluster_A |
| prot_B | cluster_A |
| prot_C | cluster_B |
| prot_D | cluster_B |

```
args_oap stage_one -i [reads_directory] -o [args_oap_out] -t 16 -f fastq --database proteins.faa
args_oap stage_two -i [args_oap_out]  -t 16 --database proteins.faa --structure1 proteins_structure.txt
```

4. Run **fasta_to_fastq.pl** to recover the reads mapping to the proteins from the args_oap output

```
perl fasta_to_fastq.pl args_oap_out/extracted.filtered.fa [RAT_1.fastq] [RAT_2.fastq] [RAT_single.fastq]
```

5. Get taxonomic classifications from the paired end reads using **RAT**

```
CAT_pack reads -c metaspades_contigs.fa -t cat_database/tax/ -m cr -o [RAT_paired] -1 [RAT_1.fastq] -2 [RAT_2.fastq] -d cat_database/db/ --no_stars -n 16 --sensitive --block_size 6 --tmpdir temp --c2c contigs_CAT.contig2classification.txt
```

6. Use the modified CAT_pack script included in this repo to run RAT using single end reads

```
CAT_pack reads -c metaspades_contigs.fa -t cat_database/tax/ -m cr -o [RAT_single] -1 [RAT_single.fastq] -d cat_database/db/ --no_stars -n 16 --sensitive --block_size 6 --tmpdir temp/ --c2c contigs_CAT.contig2classification.txt
```

7. Run **SingleM** in microbial_fraction mode using all the metagenomic reads (NOT THE RAT READS) to get the estimated prokaryote genome size and number of prokaryotic bases

```
singlem pipe -1 [sample_name_1.fastq.gz] -2 [sample_name_2.fastq.gz] -p [taxonomic_profile] --otu-table [OTU_table] --threads 16
singlem microbial_fraction -1 [sample_name_1.fastq.gz] -2 [sample_name_2.fastq.gz] -p [taxonomic profile] > sample_name_smf
```
8. Create a **MMSeqs2** database for your proteins

```
mmseqs createdb proteins.faa proteins_DB
```

9. Run **MMSeqs2** using the RAT reads as queries against the proteins

```
mmseqs easy-search [RAT_1.fastq] proteins_DB [RAT_1.m8] tmp --alignment-mode 3 -s 7 --format-output "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen" --remove-tmp-files
mmseqs easy-search [RAT_2.fastq] proteins_DB [RAT_2.m8] tmp --alignment-mode 3 -s 7 --format-output "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen" --remove-tmp-files
mmseqs easy-search [RAT_single.fastq] proteins_DB [RAT_single.m8] tmp --alignment-mode 3 -s 7 --format-output "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen" --remove-tmp-files
```

10. Concatenate files

```
cat RAT_*.m8 > reads.m8
cat *read2classification.txt > read2classification.txt
```

11. Run **get_abundance_and_taxonomy.pl** to quantify proteins in copies/cell and get taxonomic classifications for each protein based on the LCA of all reads that map to them. The RAT_otu_table output is useful if you want to calculate alpha-diversity metrics for each protein afterwards.

```
perl get_abundance_and_taxonomy.pl [sample_name_smf] [read2classification.txt] [reads.m8] [RAT_otu_table] > RAT_abundance_and_tax
```
