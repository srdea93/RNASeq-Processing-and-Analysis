#!/bin/bash
# Individual RNA-seq processing - paired script version 1.0
# Author: Steven Dea
# Date Last Modified: 6/8/20
# Required tools and packages:  conda - fastqc, hisat2, samtools, stringtie


# This script will run RNA-seq data pipeline within a project directory. It will only perform the following steps
# on a single fastq file. It is incorporated into an outer running script that iterates through all files.
# It will run fastqc, hisat2 alignment, samtools sam to bam/sort, and stringtie to generate expression estimates
# It will have four command line arguments:
# 1. Fastq R1 file to be analyzed ($1)
# 2. Fastq R2 file to be analyzed ($2)
# 3. Path to the .gtf file for Hg38 - this file is used for stringtie and featurecounts and should contain gene names ($3)
# 4. Path to the .fa reference genome file for Hg38 - this file is used for HISAT2 alignment ($4)
# 5. # of threads/cores to use for processing ($5)
# positional parameter inputs: 1. fastq file 2. .gtf 3. reference genome indexed file
# command-line example: RNAseqScript fastq_file gtf_file reference_genome_indexed_file

# Begin timing of the script
start=`date +%s`

# Assign first positional paramter to variable fq
fq1=${1?Error: no fastq1 file specified}

fq2=${2?Error: no fastq2 file specified}

# Grab the base of the file name for output naming
base=`basename $fq1 _1.fastq`
base2=`basename $fq2 _2.fastq`

# Assigns third positional paramter to the variable gtf
gtf=${3?Error: no .gtf file specified}

# Assigns fourth positional paramter to the variable ref
ref=${4?Error: no reference genome specified}

# Assigns fifth positional parameter to variable cores, if none is specified default to 1
cores=${5:-1}

#printf "\nFastq file to be processed: $fq"
#printf "\nBasename of file to be processed: $base"
#printf "\nGTF file specified: $gtf"
#printf "\nREF file specified: $ref"
#printf "\nCores to be used: $cores\n"

# Generate first level output directories using -p -- only creates a directory if it does not already exist
mkdir -p -m777 FASTQ_DIR
mkdir -p -m777 ALIGNED_DIR
mkdir -p -m777 ALIGNED_DIR/fastqc
mkdir -p -m777 ALIGNED_DIR/$base
mkdir -p -m777 EXPRESSION_DIR
# mkdir -p -m777 EXPRESSION_DIR/$base/featurecounts
mkdir -p -m777 EXPRESSION_DIR/$base/stringtie
mkdir -p -m777 EXPRESSION_DIR/$base/htseq_counts
mkdir -p -m777 EXPRESSION_DIR/all_samples_ST
mkdir -p -m777 FASTQ_DIR/fastqc


adir=ALIGNED_DIR
adirfq=ALIGNED_DIR/fastqc
adirbase=ALIGNED_DIR/$base
edir=EXPRESSION_DIR
edirbase=EXPRESSION_DIR/$base
# edirbasefc=EXPRESSION_DIR/$base/featurecounts
edirbasest=EXPRESSION_DIR/$base/stringtie
edirbasehtc=EXPRESSION_DIR/$base/htseq_counts
edirasst=EXPRESSION_DIR/all_samples_ST
fqdir=FASTQ_DIR
fqdirfqc=FASTQ_DIR/fastqc


printf "\nProcessing file: $fq1 and $fq2"
printf "\nPerforming fastqc on file: $fq1 and $fq2"

# Fastqc on our fastq files
fastqc $fq1
fastqc $fq2
# Move fastqc files to the FASTQ_DIR fastqc subdirectory
printf "\nMoving fastqc files to FASTQ_DIR\n"
mv $base/*.html $fqdirfqc
mv $base/*.zip $fqdirfqc
mv $base2/*.html $fqdirfqc
mv $base2/*.zip $fqdirfqc


printf "\nPerforming HISAT2 Alignment"
printf "\nReference Index file: $ref"
printf "\nInput file 1: $fq1"
printf "\nInput file 2: $fq2"
printf "\nOutput file: ALIGNED_DIR/$base/$base.sam\n"
# Run HISAT2 Alignment on the samples
# ‘-p’ tells HISAT2 how many cores to use.
# ’–rna-strandness RF’ specifies strandness of RNAseq library. We will specify RF since the TruSeq strand-specific library was used to make these libraries.
# ’–rg-id $ID’ specifies a read group ID that is a unique identifier.
# ’–rg SM:$SAMPLE_NAME’ specifies a read group sample name. This together with rg-id will allow you to determine which reads came from which sample in the merged bam later on.
# ’–rg LB:$LIBRARY_NAME’ specifies a read group library name. This together with rg-id will allow you to determine which reads came from which library in the merged bam later on.
# ’–rg PL:ILLUMINA’ specifies a read group sequencing platform.
# ’–rg PU:$PLATFORM_UNIT’ specifies a read group sequencing platform unit. Typically this consists of FLOWCELL-BARCODE.LANE
# ’–dta’ Reports alignments tailored for transcript assemblers.
# ‘-x /path/to/hisat2/index’ The HISAT2 index filename prefix (minus the trailing .X.ht2) built earlier including splice sites and exons.
# ‘-1 /path/to/read1.fastq.gz’ The read 1 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
# ‘-2 /path/to/read2.fastq.gz’ The read 2 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
# ‘-S /path/to/output.sam’ The output SAM format text file of alignments.
hisat2 -p $cores --rg-id=$base --rg PL:ILLUMINA -x $ref -1 $fq1 -2 $fq2 -S $adirbase/$base.sam
printf "\nHISAT2 Alignment of $fq1 and $fq2 completed\n"

# SAM to BAM Conversion
# -@ sets the number of sorting and compression threads. Default is 1.
# -o write the final sorted output to a file instead of to a standard output
printf "\nSorting aligned $base.sam"
samtools sort -@ $cores -o $adirbase/$base.bam $adirbase/$base.sam
printf "\nSorting of $base .sam and .bam completed\n"

# Index sorted bam files to allow for visualization in IGV
printf "\nIndexing sorted $base.bam"
samtools index $adirbase/*.bam
printf "\nIndexing sorted $base.bam completed\n"


printf "\nPerforming Stringtie on $base.bam"
printf "\nTranscript Output: EXPRESSION_DIR/$base/stringtie/transcripts.gtf"
printf "\nGene Abundance Output: EXPRESSION_DIR/$base/stringtie/gene_abundances.tsv\n"
# Stringtie
# ‘-p ’ tells Stringtie how many cores to use
# ‘-G ' reference annotation to use for guiding the assembly process (GTF/GFF3)
# ‘-e’ only estimate the abundance of given reference transcripts (requires -G)
# ‘-B’ enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)
# ‘-o’ output path/file name for the assembled transcripts GTF (default: stdout)
# ‘-A’ output path/file name for gene abundance estimates
stringtie -p $cores -G $gtf -e -B -o $edirbasest/transcripts.gtf -A $edirbasest/gene_abundances.tsv $adirbase/$base.bam
printf "\nCopying files to all_samples_ST/$base"
cp -r $edirbasest $edirasst/$base
printf "\nStringtie on $base.bam completed\n"

# End timing of the script
end=`date +%s`
runtime=$((end-start))
printf "\nTotal runtime of script: $runtime seconds\n"

# End individual fastq processing - paired script




