#!/usr/bin/env python
# Make sure the shebang is the working Python installation.
# This python script is designed to be run on the command line with a list of input files as an argument.
# It will call mygene to find gene symbols and names from gene IDs and add columns containing both to
# the resulting data frame and output as a tsv or csv.
# A little too hard-coded and only will work for gene/transcript expression matrices from stringtie
# Author: Steven Dea
# Version 1.0
# Last updated: 5/26/2020

import mygene
import pandas as pd
import argparse
def arg_parse():
    parser = argparse.ArgumentParser(description="Process all input files and output file")
    parser.add_argument('i', metavar='input files', nargs='+', help="List of input files to "
                                                         "be run through Get_Gene_Names")
    parser.add_argument('scope', metavar='reference scope', nargs=1, help="Scope of which gene reference to search. "
                                                                          "Examples: 'refseq' ,'ensembl.gene'")

    return parser

def get_gene_names(input, scope):
    "This function takes all input .csv or .tsv files and retrieves gene names and ids using mygene"

    input1 = input[0].replace(',','')
    # Loads the input file into a pandas data frame
    df1 = pd.read_csv(input1, sep='\t')

    # Generate a gene ID list from the data frame
    # Handle if we are working with genes or transcripts
    if "gene" in input[0]:
        gene_ID_list = df1['Gene_ID'].tolist()
        symbol = "Gene Symbol"
        name = "Gene Name"
    else:
        gene_ID_list = df1['Transcript_ID'].tolist()
        symbol = "Transcript Symbol"
        name = "Transcript Name"

    # Initialize mygene
    mg = mygene.MyGeneInfo()

    # Create two lists to store all names and symbols
    gene_name_list = []
    gene_symbol_list = []

    # Iterate through all gene IDs to get gene names and symbols
    mg_list = mg.querymany(gene_ID_list[:], scopes=scope, returnall=True)
    for mg in mg_list['out']:
        try:
            gene_name_list.append(mg['name'])
            gene_symbol_list.append(mg['symbol'])
    # If no name or symbol is found, insert - into the data frame
        except:
            gene_name_list.append(" - ")
            gene_symbol_list.append(" - ")
            # print("No value found.")

    # Add gene symbol and gene name to the each input dataframe on columns 1 and 2
    for file in input:
        file = file.replace(',','')
        df = pd.read_csv(file, sep='\t')
        df.insert(1, symbol, gene_symbol_list)
        df.insert(2, name, gene_name_list)
        # output the new dataframe as a tsv or csv
        if (file[-4:] == ".tsv"):
            df.to_csv(file, sep='\t')
        else:
            df.to_csv(file)

# Call parser
parser = arg_parse()
args = parser.parse_args()
# Set up input list
input = args.i
scope = args.scope

# Call get_gene_names on input file list
get_gene_names(input, scope)

# END PROGRAM