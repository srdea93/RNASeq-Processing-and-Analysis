#!/usr/bin/env python
# This program takes in a text file of genes of interest (txt file of gene symbols separated by a ',', a
# DESeq2 .tsv file, and an output file and extracts all information regarding the specified genes
#
# Command line Example:
# DEG_Extract_Genes_of_Interest.py path/to/input_DGE_File path/to/text_file_of_genes path/to/output_file

# Author: Steven Dea
# Version: 1.0
# Date last modified: 6/3/20

import pandas as pd
import argparse


def arg_parse():
    parser = argparse.ArgumentParser(description="Extract Genes of Interest from a DESeq2 .DGE file")
    parser.add_argument('i', metavar='Input File', nargs=1, help="DESeq2 input .DGE file")
    parser.add_argument('t', metavar='Text File', nargs=1,
                        help="Genes of interest .txt file (gene symbols separated by ','")
    parser.add_argument('o', metavar="Output File", nargs=1, help="Output Genes of Interest .tsv file")
    return parser


class GenesOfInterest:
    def __init__(self, text_file, DESeq2_file, output_file):
        self.text_file = text_file
        self.DESeq2_file = DESeq2_file
        self.output_file = output_file
        self.GOI_list = None
        self.DESeq2_df = None
        self.GOI_df = None

        self.read_txt()
        self.read_deseq2_file()
        self.extract_goi()
        self.write_goi_df(self.output_file)

    def read_txt(self):
        file = open(self.text_file, "r")
        file_str = file.read()
        self.GOI_list = file_str.split(",")

    def read_deseq2_file(self):
        if self.DESeq2_file[-4:] == ".tsv":
            df = pd.read_csv(self.DESeq2_file, sep="\t")
        else:
            df = pd.read_csv(self.DESeq2_file, sep=",")
        # drop extra index column and all values with NaN
        df.drop(df.columns[0], axis=1, inplace=True)
        df.dropna(inplace=True)
        # handle data with p-adj and p-values values of 0 by replacing them with the next lowest p - just as significant
        df['padj'].replace({0:1.0e-300}, inplace=True)
        df['pvalue'].replace({0:1.0e-300}, inplace=True)
        self.DESeq2_df = df

    def extract_goi(self):
        self.GOI_df = self.DESeq2_df[self.DESeq2_df['Gene Symbol'].isin(self.GOI_list)]

    def write_goi_df(self, output_file):
        if output_file[-4:] == ".tsv":
            self.GOI_df.to_csv(output_file, sep='\t')
        else:
            self.GOI_df.to_csv(output_file, sep=',')

def main():
    # Call parser
    parser = arg_parse()
    args = parser.parse_args()
    # Set up input list
    input_file = args.i[0]
    text_file = args.t[0]
    output_file = args.o[0]

    genes = GenesOfInterest(text_file, input_file, output_file)


# **************************** TEST PROCESS ****************************


# **************************** RUN PROCESS ****************************
if __name__ == "__main__":
    main()