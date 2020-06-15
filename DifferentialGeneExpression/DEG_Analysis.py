#!/usr/bin/env python
# This module uses the DEG_Object defined in another module to generate multiple DEG objects from all of the
# DESeq2 matrices generated from each group comparison. It will first iterate through a directory to find all DESeq2
# files and turn them into DEG_Objects. It will then extract data from each object regard differentially expressed genes
# and generate a DEG count graph, a Venn diagram showing all DEGs found in each sample grouping, and perform
# GO-enrichment analysis.

# Examples of command line usage:
# DEG_Analysis.py path/to/DGE_Directory log2foldchange_threshold_value padj_threshold_value [gene sets to analyze]
# DEG_Analysis.py DGE_subset_test/
# DEG_Analysis.py DGE_subset_test/ 2 0.01 Human_Gene_Atlas GO_Biological_Process_2018

# Run this program from the command line: only REQUIRED input needed is the path to the directory containing the DGE.tsv
# files from DESeq2 output to be analyzed as a batch. However, can specify threshold values for log2 fold change,
# adjusted p-value, and the gene sets to be analyzed. Otherwise will use defaults.

# Author: Steven Dea
# Version: 1.0
# Date last modified: 6/9/20

import pandas as pd
import seaborn as sns
import os, fnmatch, errno
from matplotlib import pyplot as plt
from venn import venn
import DEG_Object
import argparse
import gseapy as gp


def arg_parse():
    parser = argparse.ArgumentParser(description="Perform group analysis of DESeq2 .DGE files")
    parser.add_argument('--i', metavar='Input Directory', nargs=1, help="Path to Directory of DESeq2 .DGE files")
    parser.add_argument('-l2fct', '--l2fct', metavar='Log2 Fold Change Threshold', nargs='?', default=2,
                        help="Log2 Fold Change Threshold for DESeq2. Default = 2")
    parser.add_argument('-adjpt', '--adjpt', metavar='Adjusted p-value Threshold', nargs='?', default=0.01,
                        help="Adjusted p-value Threshold for DESeq2. Default = 0.01")
    parser.add_argument('-gs', '--gs', metavar='Gene Sets', nargs='*',
                        default=['Human_Gene_Atlas','GO_Biological_Process_2018','GO_Cellular_Component_2018',
                                 'GO_Molecular_Function_2018'],
                        help="Gene Sets to analyze for enrichment - Examples: 'Human_Gene_Atlas', "
                             "'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018'")

    return parser


# Finds all files matching a specific pattern within a directory and returns a list of strings
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


class DEGAnalysis:
    def __init__(self, input_dir, gene_sets=('Human_Gene_Atlas', 'GO_Biological_Process_2018'), l2fc=2, padj=0.01):
        self.input_dir = input_dir
        self.l2fc = l2fc
        self.padj = padj
        self.gene_sets = gene_sets
        self.DEG_object_list = None
        self.DEG_total_count_list = None
        self.DEG_all_counts_df = None
        self.set_all()

    def set_all(self):
        """Method to set all of the member variables of the DEGAnalysis object. Only require the input directory
        to set the rest of the member data."""
        DEG_file_list = sorted(find("*DGE.tsv", self.input_dir))
        DEG_total_count_list = []
        DEG_count_df_list = []
        DEG_object_list = []
        # iterate through all files and generate DEG_objects for them
        for file in DEG_file_list:
            deg_obj = DEG_Object.DEGObject(file, output_dir=self.input_dir, l2fc_thr=self.l2fc, pv_thr=self.padj)
            deg_obj.set_all()
            deg_obj.generate_all_plots()
            deg_obj.write_df(deg_obj.deg_df)
            DEG_total_count_list.append(deg_obj.total_counts)
            DEG_count_df_list.append(deg_obj.count_df)
            DEG_object_list.append(deg_obj)

        DEG_count_df = (pd.concat(DEG_count_df_list)).sort_index()

        self.DEG_object_list = DEG_object_list
        self.DEG_total_count_list = DEG_total_count_list
        self.DEG_all_counts_df = DEG_count_df

    def horizontal_barplot(self):
        """Generates a horizontal bar plot of all of the DGE files in the input directory. The total count is
        represented by the entirety of the bar, while the up regulated gene counts are represented by the color blue
        and the down regulated genes are represented by the color red"""
        plt.figure()
        sns.set(style="whitegrid")
        # Initialize the matplotlib figure
        f, ax = plt.subplots()

        counts = self.DEG_all_counts_df

        # Plot the total Counts
        sns.set_color_codes("pastel")
        count_hor_barplot = sns.barplot(x="Total Counts", y=counts.index, data=counts,
                    label="Down-regulated Gene Counts", color="r")

        # Plot the crashes where alcohol was involved
        sns.set_color_codes("muted")
        sns.barplot(x="Upregulated Gene Counts", y=counts.index, data=counts,
                    label="Up-regulated Gene Counts", color="b")

        # Add a legend and informative axis label
        ax.legend(loc="best", frameon=True, fontsize=5)
        ax.set(title="DEG Counts", ylabel="Comparisons",
               xlabel="Differentially Expressed Gene Counts per Sample Group Comparison")
        sns.despine(left=True, bottom=True)
        count_hor_barplot_figure = count_hor_barplot.get_figure()
        count_hor_barplot_figure.savefig(self.input_dir + "horizontal_barplot_DEG.png", dpi=500)
        plt.close()

    def venn_diagram(self, deg_objects_input_list):
        """Generates 3 venn diagrams for all of the DGE files in the input directory. The first venn diagram is the
        shared genes between each sample comparison for all comparisons. The second is the shared up regulated genes
        and the third is the shared down regulated genes

        CANNOT run on groups of more than 6 inputs"""
        # check if length of input list is > 6
        if len(deg_objects_input_list) > 6 or len(deg_objects_input_list) < 2:
            print("Error, venn diagram input list must be between 2 and 6 sample comparison groups")
            pass
        else:
            # Produce a venn diagram for all DEGs
            # extract deg index numbers from each of the deg_objects
            deg_index_list = []
            deg_up_index_list = []
            deg_down_index_list = []
            for deg_obj in deg_objects_input_list:
                # Extract ALL DEGs
                deg_indexes = set(deg_obj.deg_df.index.tolist())
                deg_index_list.append(deg_indexes)

                # Extract up regulated DEGs
                deg_up = set(deg_obj.deg_df[deg_obj.deg_df['log2FoldChange'] < 0].index.tolist())
                deg_up_index_list.append(deg_up)

                # Extract down regulated DEGs
                deg_down = set(deg_obj.deg_df[deg_obj.deg_df['log2FoldChange'] > 0].index.tolist())
                deg_down_index_list.append(deg_down)

            # create a list of indexes from count_df
            group_list = self.DEG_all_counts_df.index.tolist()
            group_set = set(group_list)

            # Create a dictionary where each index of all counts data frame is paired with a list of deg indexes
            venn_dict_all = dict(zip(group_set, deg_index_list))
            # print(venn_dict_all)
            # call venn on this dict
            plt.figure()
            venn_diagram = venn(venn_dict_all, cmap='plasma')
            venn_diagram.set_title("All Differentially Expressed Genes")
            # venn_diagram_fig = venn_diagram.get_figure()
            plt.savefig(self.input_dir + "venn_DEG.png", dpi=500)
            plt.close()

            # Produce a venn diagram for all up regulated genes
            plt.figure()
            venn_dict_up = dict(zip(group_set, deg_up_index_list))
            venn_diagram_up = venn(venn_dict_up, cmap='plasma')
            venn_diagram_up.set_title("Up-Regulated Differentially Expressed Genes")
            # venn_diagram_up_fig = venn_diagram_up.get_figure()
            plt.savefig(self.input_dir + "venn_up_DEG.png", dpi=500)
            plt.close()

            # Produce a venn diagram for all down regulated genes
            venn_dict_down = dict(zip(group_set, deg_down_index_list))
            venn_diagram_down = venn(venn_dict_down, cmap='plasma')
            venn_diagram_down.set_title("Down-Regulated Differentially Expressed Genes")
            # venn_diagram_down_fig = venn_diagram_down.get_figure()
            plt.savefig(self.input_dir + "venn_down_DEG.png", dpi=500)
            plt.close()

    def enrichr_barplot(self, data_frame, gene_set, reg):
        """Generates a barplot from enrichr data of all combined enriched upregulated or down regulated
        GO terms from different GO term libraries (gene sets)"""
        sns.set_palette("husl")
        sns.set(font_scale=0.5)
        plt.figure()
        barplot = sns.barplot(x='Combined Score', y='Term', hue='Sample', data=data_frame)
        ax = barplot.axes
        ax.legend(loc="best", frameon=True, fontsize=5)
        ax.set(title=gene_set, ylabel="Term",
               xlabel="Combined Scores")
        plt.setp(ax.patches, linewidth=0)
        plt.savefig(self.input_dir + "enrichr/" + reg + "_barplot_" + gene_set, dpi=500, bbox_inches='tight',
                                 pad_inches=0.2)
        plt.close()

    def enrichr(self, gene_dict, gene_set, key, reg):
        """Perform enrichr analysis on a gene dictionary of sample group : enriched gene list, a GO term library
        of gene set, and a reg varible that is either 'upreg' or 'downreg'"""
        # check if gene list is empty
        if not gene_dict[key]:
            pass
        else:
            # run enrichr - if there are no genes enriched with the cutoff level, it will not generate an output
            enr = gp.enrichr(gene_list=gene_dict[key],
                                gene_sets=[gene_set],
                                organism='Human',
                                description=key + "_" + reg,
                                outdir=self.input_dir + 'enrichr/' + key + "/",
                                cutoff=0.1
                                )
            enr_df = enr.results.copy()
            return enr_df

    def go_gene_enrichment_analysis(self, gene_sets):
        """Performs GO enrichment analysis with an input list of GO term libraries to run the process for (i.e.
        Human Gene Atlas, etc). Outputs individual up and down reg enrichment .txt files and pdfs of top enriched
        GO terms that maintain an adjusted p-value of < 0.05 for each comparison group. Also outputs group summary
        horizontal bar plots between all groups with top 5 GO terms that maintain an adjusted p-value of < 0.05"""
        # get up regulated gene list for analysis
        deg_up_dict = {}
        deg_down_dict = {}
        for deg_obj in self.DEG_object_list:
            # Extract up regulated DEGs
            deg_up = deg_obj.deg_df[deg_obj.deg_df['log2FoldChange'] < 0]['Gene Symbol'].tolist()
            deg_up_dict.update({deg_obj.name: deg_up})

            # Extract down regulated DEGs
            deg_down = deg_obj.deg_df[deg_obj.deg_df['log2FoldChange'] > 0]['Gene Symbol'].tolist()
            deg_down_dict.update({deg_obj.name: deg_down})

        # iterate through each gene set library desired i.e. Human_Gene_Atlas
        for gene_set in gene_sets:
            enr_list_up = []
            enr_list_down = []
            # gseapy for enrichr analysis for up genes
            for key in list(deg_up_dict.keys()):
                enr_up_df = self.enrichr(deg_up_dict, gene_set, key, "upreg")
                if enr_up_df is None:
                    pass
                else:
                    df_up = enr_up_df
                    df_up = df_up.sort_values(by=['Adjusted P-value'], ascending=True)
                    df_up['Sample'] = key + "_upreg"
                    # Only select the first 5 terms to be plotted
                    df_up = df_up[:5]
                    df_up = df_up[['Term', 'Combined Score', 'Adjusted P-value', 'Sample']]
                    enr_list_up.append(df_up)
            # gseapy for enrichr analysis for down genes for human gene atlas
            for key in list(deg_down_dict.keys()):
                enr_down_df = self.enrichr(deg_down_dict, gene_set, key, "downreg")
                if enr_down_df is None:
                    pass
                else:
                    df_down = enr_down_df
                    df_down = df_down.sort_values(by=['Adjusted P-value'], ascending=True)
                    df_down['Sample'] = key + "_downreg"
                    # Only select the first 5 terms to be plotted
                    df_down = df_down[:5]
                    df_down = df_down[['Term', 'Combined Score', 'Adjusted P-value', 'Sample']]
                    enr_list_down.append(df_down)

            # concat all of the up regulated data frames
            barplot_df_up = pd.concat(enr_list_up)
            barplot_df_up = barplot_df_up[barplot_df_up['Adjusted P-value'] <= 0.05]

            # horizontal bar plot of all combined score data for
            self.enrichr_barplot(barplot_df_up, gene_set, "upreg")

            # concat all of the down regulated data frames
            barplot_df_down = pd.concat(enr_list_down)
            barplot_df_down = barplot_df_down[barplot_df_down['Adjusted P-value'] <= 0.05]

            # horizontal bar plot of all combined score data for
            self.enrichr_barplot(barplot_df_down, gene_set, "downreg")

    def write_up_and_down_list(self):
        """Writes output .txt files of up regulated gene lists and down regulated gene lists. It can be used
        to run online GO analysis to compliment the current GO analysis. (http://bioinformatics.sdstate.edu/go/)
        """
        # up regulated gene lists
        dirname = self.input_dir + "up_down_reg_lists/"
        if not os.path.exists(os.path.dirname(dirname)):
            try:
                os.makedirs(os.path.dirname(dirname))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        for deg_obj in self.DEG_object_list:
            output_file_up = dirname + deg_obj.name + "_up.txt" # need to write outputs
            with open(output_file_up, 'w+') as f_up:
                f_up.write('\n'.join(deg_obj.up_list))
            output_file_down = dirname + deg_obj.name + "_down.txt"  # need to write outputs
            with open(output_file_down, 'w+') as f_down:
                f_down.write('\n'.join(deg_obj.down_list))

    def print(self):
        """Prints some of the important member variables of the DEGAnalysis object"""
        print("Object List")
        print(self.DEG_object_list)
        print("Total Count List")
        print(self.DEG_total_count_list)
        print("All Counts Data Frame")
        print(self.DEG_all_counts_df)

    def print_deg_dfs(self):
        """Prints the deg_df for each of the deg objects stored in the DEG analysis object"""
        for deg_obj in self.DEG_object_list:
            print(deg_obj.deg_df)

def main():
    # Call parser
    parser = arg_parse()
    args = parser.parse_args()
    # Set up input list
    input_dir = args.i[0]

    try:
        l2fc = int(args.l2fct[0])
        padj = float(args.adjpt)
    except:
        print("No values or incorrect values passed for thresholds. Resorting to defaults.")
        l2fc = 2
        padj = 0.01
    try:
        gene_sets = args.gs
    except:
        print("No values specified for gene_sets. Resorting to defaults.")
        gene_sets = ['Human_Gene_Atlas', 'GO_Biological_Process_2018',
                     'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018']

    print("Log2 Fold Change Threshold: " + str(l2fc))
    print("Adjusted p-value Threshold: " + str(padj))
    print("\nGene Sets: ")
    for gene in gene_sets:
        print(gene)


    deg_analysis = DEGAnalysis(input_dir, gene_sets=gene_sets, l2fc=l2fc, padj=padj)
    deg_analysis.write_up_and_down_list()
    deg_analysis.horizontal_barplot()
    deg_analysis.go_gene_enrichment_analysis(gene_sets)


# **************************** TEST PROCESS ****************************



# **************************** RUN PROCESS ****************************
if __name__ == "__main__":
    main()
