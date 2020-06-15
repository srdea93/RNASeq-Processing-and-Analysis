#!/usr/bin/env python
# This module defines a Differential Expressed Genes object that contains all of the functions to create DEG objects
# That can be used in other modules. It has functions to create itself using an input file and a set_all() function
# that populates the object full of data from the input file. It has functions to generate and save different types of
# plots such as volcano, and distribution plots. It also generates DEG counts data frame that can be used in conjunction
# with other DEG objects to get larger data statistics on multiple sample group DGE analysis
# Author: Steven Dea
# Version: 1.0
# Date last modified: 6/9/20

import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
import bioinfokit
from bioinfokit import analys,visuz
import shutil, os
import argparse


def arg_parse():
    parser = argparse.ArgumentParser(description="Perform Analysis on DESeq2 .DGE file")
    parser.add_argument('i', metavar='Input File', nargs=1, help="Input DESeq2 .DGE file")
    parser.add_argument('o', metavar='Output Directory', nargs=1, help="Output Directory")
    parser.add_argument('l2fct', metavar='Log2 Fold Change Threshold', nargs=1, help="Log2 Fold Change Threshold"
                                                                                     "for DESeq2. Default = 2")
    parser.add_argument('adjpt', metavar='Adjusted p-value Threshold', nargs=1, help="Adjusted p-value Threshold"
                                                                                     "for DESeq2. Default = 0.01")
    return parser


class DEGObject:
    def __init__(self, input_file, output_dir=os.getcwd(), l2fc_thr=2, pv_thr=0.01):
        self.input_file = input_file
        self.output_dir = output_dir
        self.name = None
        self.df = None
        self.l2fc_thr = l2fc_thr
        self.pv_thr = pv_thr
        self.deg_df = None
        self.up_counts = None
        self.up_list = None
        self.down_counts = None
        self.down_list = None
        self.total_counts = None
        self.count_df = None
        self.l2fc_df = None

        self.set_all()

        if not os.path.exists(output_dir + "pairwise_summaries/"):
            os.mkdir(output_dir + "pairwise_summaries/")
        else:
            pass

        self.output_dir = output_dir + "pairwise_summaries/"

    def set_name(self):
        """Sets the name of the DEG obj based on the input file"""
        split_file = self.input_file.split('/')
        file_name = split_file[-1]
        file_name = file_name.split('_')
        class_name = str(file_name[1] + "_" + str(file_name[2] + "_" + str(file_name[4])))
        self.name = class_name

    def set_df(self):
        """Takes an input file as a .tsv or .csv and returns it as a pandas data frame

            Arguments:
                input_file -- input file in the form of a .tsv or a .csv"""

        if self.input_file[-4:] == ".tsv":
            df = pd.read_csv(self.input_file, sep="\t")
        else:
            df = pd.read_csv(self.input_file, sep=",")
        # drop extra index column and all values with NaN
        df.drop(df.columns[0], axis=1, inplace=True)
        df.dropna(inplace=True)
        # handle data with p-adj and p-values values of 0 by replacing them with the next lowest p - just as significant
        df['padj'].replace({0:1.0e-300}, inplace=True)
        df['pvalue'].replace({0:1.0e-300}, inplace=True)
        self.df = df

    def set_deg_df(self):
        """Sets the differentially expressed gene data frame based on the log2 Fold Change threshold
        and the adjusted p-value threshold. It then sorts the resulting DEG data frame in descending order by
        log2 Fold Change"""
        # How many differentially expressed genes are identified from DESeq2: log2fold change & p-adj
        log2foldchange_df = self.df[(abs(self.df['log2FoldChange']) >= self.l2fc_thr)]
        log2foldchange_df = log2foldchange_df[log2foldchange_df['padj'] <= self.pv_thr]
        self.deg_df = log2foldchange_df.sort_values('log2FoldChange', ascending=False)

    def set_count_df(self):
        """"Sets the count data frame for down regulated genes and up regulated genes for use in
        later visualizations with all other sample group comparisons"""
        self.up_counts = self.deg_df[self.deg_df['log2FoldChange'] < 0].count()[0]
        self.up_list = self.deg_df[self.deg_df['log2FoldChange'] < 0]['Gene Symbol']
        self.down_counts = self.deg_df[self.deg_df['log2FoldChange'] > 0].count()[0]
        self.down_list = self.deg_df[self.deg_df['log2FoldChange'] > 0]['Gene Symbol']
        self.total_counts = self.up_counts + self.down_counts
        data = {"Upregulated Gene Counts": [self.up_counts], "Downregulated Gene Counts": [self.down_counts],
                "Total Counts" : [self.total_counts]}
        deg_counts = pd.DataFrame(data, index=[self.name])
        self.count_df = deg_counts

    def set_all(self):
        """Calls all other set methods to set them in one."""
        self.set_name()
        self.set_df()
        self.set_deg_df()
        self.set_count_df()

    # Need to move the output of this using shutil because it doesn't allow to specify a output destination
    def volcano_plot_padj(self, num_DEGs=50):
        """Generates a volcano plot based on the log2 Fold Change threshold and adjusted p-value threshold.
        Automatically populates the plot with the top 50 DEG based on p-adjusted values."""
        # Volcano plot https://reneshbedre.github.io/blog/volcano.html
        # sort the data frame by padj and then get the top DEGs passed into the function
        plt.figure()
        padj_df = self.df = self.df.sort_values('padj', ascending=True)
        top_padj = padj_df[abs(padj_df['log2FoldChange']) >= self.l2fc_thr]
        top_padj = tuple(top_padj['Gene Symbol'][:num_DEGs].tolist())
        visuz.gene_exp.volcano(d=self.df, lfc='log2FoldChange', lfc_thr=self.l2fc_thr,
                               pv='pvalue', pv_thr=self.pv_thr, r=400, geneid='Gene Symbol',
                               genenames=top_padj, gstyle=2, sign_line=True, dim=(15, 15))
        output_dir = self.output_dir + self.name + "/" + self.name + "_padj_volcano.png"
        shutil.move("volcano.png", output_dir)
        plt.close()

    # Need to move the output of this using shutil because it doesn't allow to specify a output destination
    def volcano_plot_l2fc(self, num_DEGs=50):
        """Generates a volcano plot based on the log2 Fold Change threshold and adjusted p-value threshold.
        Automatically populates the plot with the top 50 DEG based on p-adjusted values."""
        # Volcano plot https://reneshbedre.github.io/blog/volcano.html
        # sort the data frame by padj and then get the top DEGs passed into the function
        plt.figure()
        l2fc_df = self.df = self.df.sort_values('log2FoldChange', ascending=False)
        top_l2fc = l2fc_df[abs(l2fc_df['padj']) >= self.pv_thr]
        top_l2fc = tuple(top_l2fc['Gene Symbol'][:num_DEGs].tolist())
        visuz.gene_exp.volcano(d=self.df, lfc='log2FoldChange', lfc_thr=self.l2fc_thr,
                               pv='pvalue', pv_thr=self.pv_thr, r=400, geneid='Gene Symbol',
                               genenames=top_l2fc, gstyle=2, sign_line=True, dim=(15, 15))
        output_dir = self.output_dir + self.name + "/" + self.name + "_l2fc_volcano.png"
        shutil.move("volcano.png", output_dir)
        plt.close()

    def l2fc_dist_plot(self):
        """Generates a log2 Fold Change distribution plot"""
        # seaborn distribution plot of log2foldchange
        plt.figure()
        sns.set(style="whitegrid")
        sns.set_color_codes("pastel")
        l2fc = sns.distplot(self.df['log2FoldChange'], hist_kws=dict(edgecolor="k", linewidth=1))
        l2fc.set(ylabel="Percentage", title="Distribution of log2 fold change")
        l2fc_fig = l2fc.get_figure()
        l2fc_fig.savefig(self.output_dir + self.name + "/" + self.name + '_l2fc.png', dpi=500)
        plt.close()

    def padj_dist_plot(self):
        """Generates a adjusted p-value distribution plot"""
        # seaborn distribution plot of p-adj
        plt.figure()
        sns.set(style="whitegrid")
        sns.set_color_codes("pastel")
        padj = sns.distplot(self.df['padj'], hist_kws=dict(edgecolor="k", linewidth=1))
        padj.set(ylabel="Percentage", title="Distribution of p-adj values")
        plt.xlim(-.2, 1.2)
        padj_fig = padj.get_figure()
        padj_fig.savefig(self.output_dir + self.name + "/" + self.name + '_padj.png', dpi=500)
        plt.close()

    def generate_all_plots(self):
        if not os.path.exists(self.output_dir + self.name):
            os.mkdir(self.output_dir + self.name)
        else:
            pass
        # generate all plots with default values
        self.volcano_plot_padj()
        self.volcano_plot_l2fc()
        self.l2fc_dist_plot()
        self.padj_dist_plot()

    def write_df(self, df):
        """Writes the specified data frame to an output file"""
        """write a pandas data frame to an output file

            Arguments:
                df: output DESeq2 data frame to be written"""
        df.to_csv(self.output_dir + self.name + "/" + self.name + "_l2fc_DEG_matrix.tsv", sep='\t')

    def print_dfs(self):
        """Prints all of the processed data frames in the object"""
        print("DEG DF")
        print(self.deg_df)
        print("Count DF")
        print(self.count_df)

def main():
    # Call parser
    parser = arg_parse()
    args = parser.parse_args()
    # Set up input and output
    input_file = args.i[0]
    output_dir = args.o[0]
    try:
        l2fc = args.l2fct[0]
        padj = args.padjt[0]
    except:
        print("No values or incorrect values passed for thresholds. Resorting to defaults.")
        l2fc = 2
        padj = 0.01

    deg_obj = DEGObject(input_file, output_dir, l2fc, padj)
    deg_obj.generate_all_plots()
    deg_obj.write_df(deg_obj.deg_df)

# **************************** TEST PROCESS ****************************


# **************************** RUN PROCESS ****************************
if __name__ == "__main__":
    main()
