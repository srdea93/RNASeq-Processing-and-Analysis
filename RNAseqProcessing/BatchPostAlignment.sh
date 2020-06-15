#!/bin/bash

# Handles all processes after HISAT2 alignment

# Assigns first positional paramter to the variable gtf
gtf=${1?Error: no .gtf file specified}

# Assigns second positional paramter to the variable ref
ref=${2?Error: no reference genome .fa specified}

# Assigns third positional parameter to variable cores, if none is specified default to 1
cores=${3:-1}

# Perform multiqc on all fastqc files
multiqc -o FASTQ_DIR FASTQ_DIR/fastqc

# Generate expression matrices for all analyzed files
mkdir -p -m777 EXPRESSION_DIR/Combined_Output
mkdir -p -m777 EXPRESSION_DIR/Combined_htseq
mkdir -p -m777 EXPRESSION_DIR/STmerged
mkdir -p -m777 EXPRESSION_DIR/Differential_Gene_Expression
exdir=EXPRESSION_DIR
exdirm=EXPRESSION_DIR/STmerged
exdirco=EXPRESSION_DIR/Combined_Output
exdirht=EXPRESSION_DIR/Combined_htseq
exdirdge=EXPRESSION_DIR/Differential_Gene_Expression

# htseq-count on all .sam files
# –format’ specify the input file format one of BAM or SAM. Since we have BAM format files, select ‘bam’ for this option.
# ’–order’ provide the expected sort order of the input file. Previously we generated position sorted BAM files so use ‘pos’.
# ’–mode’ determines how to deal with reads that overlap more than one feature. We believe the ‘intersection-strict’ mode is best.
# ’–stranded’ specifies whether data is stranded or not. The TruSeq strand-specific RNA libraries suggest the ‘reverse’ option for this parameter.
# ’–minaqual’ will skip all reads with alignment quality lower than the given minimum value
# ’–type’ specifies the feature type (3rd column in GFF file) to be used. (default, suitable for RNA-Seq and Ensembl GTF files: exon)
# ’–idattr’ The feature ID used to identify the counts in the output table. The default, suitable for RNA-SEq and Ensembl GTF files, is gene_id.
htseqlist=`ls -dm ALIGNED_DIR/*/*.bam | tr ',' ' '`
printf "\nPerforming htseq-count on \n$htseqlist"
htseq-count --format bam --order pos --mode intersection-strict --minaqual 1 --type exon --idattr gene_id $htseqlist $gtf > $exdirht/combined_htseq_counts.tsv

# Use Htseq_Counts_Header.py to add sample header to the .tsv file
Htseq_Counts_Header.py $htseqlist $exdirht/combined_htseq_counts.tsv
printf "\nhtseq-count completed"
printf "\nPerforming DESeq2 Analysis...\n"

../Differential $exdirht/combined_htseq_counts.tsv $exdirdge/
printf "\nDESeq2 Analysis Complete."

## Unnecessary
## Stringtie merge to combine all of the transcript.gtf files from the individual stringtie assembly runs for a consensus .gtf file
## use find . -name '*.gtf' to get all of the .gtf files from each stringtie sample
#printf "\nRunning stringtie --merge...\n"
#find $exdir -name '*.gtf' > $exdirm/stringtie_gtf_list.txt
#stringtie --merge -p $cores -G $gtf -o $exdirm/stringtie_merged.gtf $exdirm/stringtie_gtf_list.txt
#
## Count how many transcripts we have in our stringtie merge.gtf
#printf "\n# of transcripts in stringtie_merged.gtf\n"
#cat $exdirm/stringtie_merged.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l

## Unnecessary
## gffcompare to compare stringtie merge transcripts to known transcripts
#gffcompare -r $gtf -G -o $exdirm/merged $exdirm/stringtie_merged.gtf
#printf "\nMerged stats:\n"
#cat $exdirm/merged.stats

# Generate expression matrices for all analyzed files
printf "\nCombining expression matrices...\n"
exdiras=EXPRESSION_DIR/all_samples_ST
exdirlist=`ls -dm $exdiras/*/ | tr -d '[:space:]'`
echo $exdirlist

# TPM
stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs="$exdirlist" --transcript_matrix_file=$exdirco/transcript_tpm_combined.tsv --gene_matrix_file=$exdirco/gene_tpm_combined.tsv
# FPKM
stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs="$exdirlist" --transcript_matrix_file=$exdirco/transcript_fpkm_combined.tsv --gene_matrix_file=$exdirco/gene_fpkm_combined.tsv
# Transcript Coverage
stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs="$exdirlist" --transcript_matrix_file=$exdirco/transcript_coverage_combined.tsv --gene_matrix_file=$exdirco/gene_coverage_combined.tsv

# Run python Get_Gene_Names.py to query entrez database using mygene and retrieve gene names and symbols for the combined data files
# Important to change python script shebang depending on the python PATH the user is using
exdircolistgene=`ls -dm $exdirco/gene* | tr ',' ' '`
exdircolisttran=`ls -dm $exdirco/transcript* | tr ',' ' '`
# rewrite input files with output files
Get_Gene_Names.py $exdircolistgene refseq
Get_Gene_Names.py $exdircolisttran refseq