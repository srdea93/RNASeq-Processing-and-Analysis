#!/bin/bash

# This script only handles the alignment process using HISAT2
# Assigns first positional paramter to the variable gtf
gtf=${1?Error: no .gtf file specified}

# Assigns second positional paramter to the variable ref
ref=${2?Error: no reference genome .fa specified}

# Assigns third positional parameter to variable cores, if none is specified default to 1
cores=${3:-1}

# Run Separate_Paired_End_Reads.py in directory containing .fq or .fastq files
Separate_Paired_End_Reads.py .

# Iterate through each directory
for d in */
do
    shopt -s nullglob
    paired_array=($d/*)
    fq1=${paired_array[0]}
    fq2=${paired_array[1]}
    echo $PWD
    IndividualRNASeqProcessingPaired.sh $fq1 $fq2 $gtf $ref $cores
done
