#!/bin/bash

### Cluster sequences obtained by metabarcoding using Usearch ###

# Usearch needs to be installed (https://www.drive5.com/usearch/) and loaded into environment
# Sequences are clustered within each individual (sample)
# Raw sequences need to be filtered and cleaned prior to processing
#
# Creates an OTU table for each individual separately

###################################################################

# Go to directory with filtered fastq 
cd /path/fastqFiles  

###################################################################

# 1. Merge forward and reverse
mkdir merged

for f in *_R1_*; do
usearch -fastq_mergepairs $f -fastqout /merged/$f -relabel @
done

# 2. Make fasta
mkdir merged_fasta
cd merged

for f in *; do
usearch -fastq_filter $f -fastaout ../merged_fasta/$f;
done

# 3. Filter for unique sequences
mkdir ../uniques
cd ../merged_fasta

for f in *; do
usearch -derep_fulllength $f -fastaout ../uniques/uniques_$f -relabel OTU_ind -sizeout;
done

# 4. Discard singletons
mkdir ../woSingletons_uniques
cd ../uniqueds

for f in *; do
usearch -sortbysize $f -fastaout ../woSingletons_uniques/woSingletons_$f -minsize 2;
done

# 5. Cluster uniques without singletons [adjust loop in naming]
mkdir ../clusters
cd ../woSingletons_uniques

for f in *; do
usearch -cluster_otus $f  -otus ../clusters/otus97_$f.fa -otu_radius_pct 3;
done

# 6 Mapping
mkdir ../maps_individuals
cd ../merged_fasta

for f in *; do
usearch -usearch_global $f -db ../clusters/otus97_woSingletons_uniques_$f.fa -strand plus -id 0.97 -otutabout ../maps_individuals/otuTable97_$f;
done

##################################################################
# 7. Identify clusters with BLAST using a species database
# Create blast species [database] and provide path after -dp
# adjust [label]

cd ../clusters/

for f in  otus97_woSingletons_uniques_[label]*; do
blastn -db [path/database] \
-query $f \
-out ../blast_results/blast_$f \
-perc_identity 90 \
-max_target_seqs 1 \
-outfmt 6 ;
done

# OTU table across individuals needs to be compiled by hand based on the blast results and sequence alignments!
