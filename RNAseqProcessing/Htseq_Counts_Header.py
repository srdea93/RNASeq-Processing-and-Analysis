#!/usr/bin/env python
# This program will be able to take htseq_counts files and combine them into a single matrix
# Author: Steven Dea
# Version: 1.0
# Date last modified: 5/26/20

import pandas as pd
import argparse

def arg_parse():
    parser = argparse.ArgumentParser(description="Add sample header to Htseq_counts file")
    parser.add_argument('i', metavar='Input File', nargs='+', help="Combined input file")
    parser.add_argument('o', metavar='Output File', nargs=1, help="Combined output file")
    return parser

def add_header(input, output):
    # Read all htseq_counts.tsv files into data frames to be used
    sample_names_list = ["Ensembl Gene ID"]
    for file in input:
        # remove the .tsv from the file names
        sample_name = file[file.rindex('/')+1:-4]
        sample_names_list.append(sample_name)

    # Add header to the input data frame - input and output files come in as a list of strings
    if output[0][-4:] == ".tsv":
        outputdf = pd.read_csv(output[0], sep='\t', names=sample_names_list)
        outputdf.columns = sample_names_list
        outputdf.to_csv(output[0], sep='\t')
    else:
        outputdf = pd.read_csv(output[0], sep=',', names=sample_names_list)
        outputdf.columns = sample_names_list
        outputdf.to_csv(output[0], sep=',')

# Call parser
parser = arg_parse()
args = parser.parse_args()
# Set up input list
input = args.i
output = args.o

# Call combine_htseq_counts on input file list
add_header(input, output)

# END PROGRAM



