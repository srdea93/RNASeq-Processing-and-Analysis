# RNASeq-Processing-and-Analysis
RNASeqProcessing contains shell scripts and Python scripts to run RNASeq processing on raw .fastq files. It can run batch inputs of unpaired or 
paired end reads and run them through fastqc, multiqc, hisat2 alignment, stringtie, and htseq-counts to generate expression matrices
and a total raw count matrix.

DifferentialGeneExpression contains Python modules that can be run on htseq-counts matrix file for an experiment and return differentially
expressed gene matrices for all pairwise comparisons via DESeq2 as well as data visualizations of GO term enrichment analysis using
enrichr.
