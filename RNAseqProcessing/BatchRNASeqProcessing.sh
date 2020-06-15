#!/bin/bash
# Group RNA-seq processing script version 1.0
# Author: Steven Dea
# Date Last Modified: 5/27/20

# Required tools and packages:  conda - multiqc, htseq-count, stringtie, gffcompare
#                               Griffith lab - stringtie_expression_matrix.pl
#                               Propietary - IndividualRNASeqProcessingUnpaired, Get_Gene_Names.py

# This script will iterate through all fastq files present in the directory it is called in
# and run the individual RNA-seq script.
# Once all files have been processed, it will perform multiqc on all fastqc files,
# Merge all HISAT2 aligned BAM files, and generate expression matrices for FPKM and TPM
# using stringtie_expression_matrix.pl from Griffith Lab
# This script has 5 command line arguments:
# 1. Path to the .gtf file for Hg38 ($1)
# 2. Path to the .fa reference genome file for Hg38 ($2)
# 3. # of threads/cores to use for processing ($3)
# 4. If data is paired end or not

# Begin timing of the script
start=`date +%s`

# Assigns first positional parameter to the variable gtf
gtf=${1?Error: no .gtf file specified}

# Assigns second positional parameter to the variable ref
ref=${2?Error: no reference genome .fa specified}

# Assigns third positional parameter to variable cores, if none is specified default to 1
cores=${3:-1}

# Assigns fourth positional parameter to variable of paired, either 1 or 2.
paired=${4?Error: Please specify whether or not data is paired or not. 1: Not paired, 2: Paired.}

if [[ $paired -eq 1 ]];
then
    # Performs batch alignment for unpaired data
    BatchAlignmentUnpaired.sh $gtf $ref $cores
elif [[ $paired -eq 2 ]];
then
    BatchAlignmentPaired.sh $gtf $ref $cores
else
    echo "paired positional parameter must be either 1 or 2"
    exit 1
fi

# Performs batch post alignment analysis
BatchPostAlignment.sh $gtf $ref $cores

# End timing of the script
end=`date +%s`
runtime=$((end-start))
printf "\nTotal runtime of batch process script: $runtime seconds\n"