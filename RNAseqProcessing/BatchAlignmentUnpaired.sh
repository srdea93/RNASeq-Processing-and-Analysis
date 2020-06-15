#!/bin/bash

# This script only handles the alignment process using HISAT2
# Assigns first positional paramter to the variable gtf
gtf=${1?Error: no .gtf file specified}

# Assigns second positional paramter to the variable ref
ref=${2?Error: no reference genome .fa specified}

# Assigns third positional parameter to variable cores, if none is specified default to 1
cores=${3:-1}

# Assigns fourth positional parameter to variable of paired, either 1 or 2.
#paired=${4?Error: Please specify whether or not data is paired or not. 1: Not paired, 2: Paired.}

# Iterate through current working directory fq files
for fq in *.fastq
do
  IndividualRNASeqProcessingUnpaired.sh $fq $gtf $ref $cores
done