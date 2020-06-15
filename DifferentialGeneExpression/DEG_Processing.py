#!/usr/bin/env python
# This program uses rpy2 and pandas to perform differential gene expression
# on htseq-count data provided from batch RNA-seq processing. This program
# implements rpy2 as an interface between Python and R to perform DGE using
# DESeq2. Portions of this code come from
# https://gist.github.com/wckdouglas/3f8fb27a3d7a1eb24c598aa04f70fb25
# Author: Steven Dea
# Version: 1.0
# Date last modified: 5/28/20

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula
pandas2ri.activate()
from rpy2.robjects.packages import importr
deseq = importr('DESeq2')
import pandas as pd
import sys
import DEG_Formatting
import argparse
import itertools


to_dataframe = ro.r('function(x) data.frame(x)')


def arg_parse():
    parser = argparse.ArgumentParser(description="Perform DESeq2 analysis on htseq-counts file")
    parser.add_argument('i', metavar='Input File', nargs=1, help="Htseq-counts input file")
    parser.add_argument('o', metavar="Output File", nargs=1, help="Output DESeq2 directory")
    return parser


# Differential Gene Expression class to store all matrices of interest to perform DESeq2
class DGE:
    """This class describes a Differential Gene Expression object to handle all DESeq2 differential gene expressions"""
    def __init__(self, count_matrix, design_matrix, design_formula, gene_column='Ensembl Gene ID'):
        try:
            assert gene_column in count_matrix.columns, "Wrong gene id column name"
            gene_id = count_matrix[gene_column]
        except AttributeError:
            sys.exit("Wrong Pandas Dataframe")

        self.dds = None
        self.deseq_result = None
        self.resLFC = None
        self.comparison = None
        self.normalized_count_matrix = None
        self.gene_column = gene_column
        self.gene_id = count_matrix[self.gene_column]
        self.count_matrix = ro.conversion.py2rpy(count_matrix.drop(gene_column, axis=1))
        self.design_matrix = ro.conversion.py2rpy(design_matrix)
        self.design_formula = Formula(design_formula)

    # Runs DESeq2 analysis and generates dds and normalized count matrix
    def run_deseq2(self, **kwargs):
        """Runs DESeq2 analysis on 2 conditions"""
        self.dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix,
                                                colData=self.design_matrix,
                                                design=self.design_formula)
        self.dds = deseq.DESeq(self.dds, **kwargs)
        # Doesn't work for some reason, but unsure if this is necessary
        # self.normalized_count_matrix = deseq.counts(self.dds, normalized=True)

    # Returns DESeq2 results
    def get_deseq_result(self, **kwargs):
        """Returns the results of the DESeq2 analysis"""
        self.comparison = deseq.resultsNames(self.dds)

        self.deseq_result = deseq.results(self.dds, **kwargs)
        self.deseq_result = to_dataframe(self.deseq_result)
        self.deseq_result = ro.conversion.rpy2py(self.deseq_result) # back to pandas data frame
        self.deseq_result[self.gene_column] = self.gene_id.values


# Main function to execute when run as a script - allows for importing of this module to other scripts
def main():
    # When running DESeq2, we must compare only 1 condition to another. I have decided to compare all

    # Call parser
    parser = arg_parse()
    args = parser.parse_args()
    # Set up input list
    input = args.i[0]
    output = args.o[0]

    total_design_matrix = DEG_Formatting.get_total_design_matrix(input)
    total_count_matrix = DEG_Formatting.get_total_count_data(input)

    all_groups = DEG_Formatting.get_groups(total_design_matrix)
    # sample_groups = DGE_Formatting.get_group_by_condition(sample_condition_start, sample_condition_end, all_groups)
    # control_groups = DGE_Formatting.get_group_by_condition(control_condition_start, control_condition_end, all_groups)

    # For loop to split design matrix into a design matrix of just 1 sample group and 1 control group - iterate through all
    # sample groups for all control groups - do the same for the count matrix and split it same as design matrix
    # For loop inside once matrices are set for DGE
    group_combinations = itertools.combinations(all_groups, 2)

    # Get gene_name_df from htseq-counts file
    htseq_counts_df = DEG_Formatting.get_df_from_file(input)
    gene_name_df = DEG_Formatting.get_gene_name(htseq_counts_df, 'Ensembl Gene ID', 'refseq')

    # runs duplicates
    for i, j in group_combinations:
        if (i == j):
            pass
        else:
            print("Performing DGE on Sample: " + str(i) + ", vs. Sample: " + str(j))
            sub_count_matrix = DEG_Formatting.get_sub_count_matrix(i, j, total_count_matrix)
            sub_design_matrix = DEG_Formatting.get_sub_design_matrix(i, j, total_design_matrix)

            # test
            print(sub_count_matrix)
            print(sub_design_matrix)

            # design formula seems to not read the input as integers for pancreas data
            dds = DGE(count_matrix=sub_count_matrix, design_matrix=sub_design_matrix, design_formula='~ Group_Number',
                      gene_column='Ensembl Gene ID')
            dds.run_deseq2()
            dds.get_deseq_result()
            res = dds.deseq_result

            res_df = pd.DataFrame(res)
            # add gene names and symbols to output data frame
            # HAVE TO RESET INDEXES OF BOTH DATA FRAMES - weird incomaptibility of pandas to concat horizontally
            output_df = pd.concat([res_df.reset_index(drop=True), gene_name_df.reset_index(drop=True)], axis=1)

            output_file = output + "SampleGroup_" + str(i) + "_vs_SampleGroup_" + str(j) + "_DESeq2_DGE.tsv"
            DEG_Formatting.write_df_to_file(output_df, output_file)

# **************************** TEST PROCESS ****************************



# **************************** RUN PROCESS ****************************

if __name__ == "__main__":
    main()

