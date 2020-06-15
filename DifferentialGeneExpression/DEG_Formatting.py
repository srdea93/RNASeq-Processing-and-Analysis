#!/usr/bin/env python
# This program is designed to format pandas data frames/matrices for input to DGE program
# Takes in an input htseq-count file and returns pandas data frames. Also allows for generation of
# sub-matrices for design_matrix and count_matrix

# IMPORTANT FORMATTING NOTE: fastq sample file names must start with #1_#2, #1 = the group number and #2 = sample number
# Author: Steven Dea
# Version: 1.0
# Date last modified: 6/1/20


import pandas as pd
import mygene


# Generates a pandas dataframe from an input file.
def get_df_from_file(input_file):
    """Takes an input file as a .tsv or .csv and returns it as a pandas data frame

    Arguments:
        input_file -- input file in the form of a .tsv or a .csv"""

    if input_file[-4:] == ".tsv":
        df = pd.read_csv(input_file, sep="\t")
    else:
        df = pd.read_csv(input_file, sep=",")
    return df


def get_groups(total_design_matrix):
    """Given an input of the total design matrix, extract the group numbers of the different sample conditions

    Arguments:
        total design matrix: input design matrix that comes from get_total_design_matrix()"""
    group_number_list = total_design_matrix['Group_Number'].tolist()
    group_list = []
    for group_number in group_number_list:
        if group_number not in group_list:
            group_list.append(group_number)
        else:
            pass
    # ******* name of the group must start with numbers******** ????
    group_list = list(map(int, group_list))
    return group_list


def get_group_by_condition(condition_start, condition_end, group_list):
    """Divide the all conditions list into sample conditions and control conditions based on input for use in
        iterating through a list of sample conditions and a list of control conditions to get all permutations
        for DESeq2 analysis

        Arguments:
            condition start number: start number of either all sample conditions or all control conditions
            condition end number: end number of either all sample conditions or all control conditions
            group list: list of all group names that is the output of get_groups()"""
    if type(condition_start) == str and type(condition_end == str):
        condition_start = int(condition_start)
        condition_end = int(condition_end)

    condition_list = []
    condition_range = list(range(condition_start, condition_end+1))
    for group in group_list:
        if group in condition_range:
            condition_list.append(group)
    return condition_list


# Generates a design matrix from an input file
def get_total_design_matrix(input_file):
    """Given an input file, generate a design matrix to be passed into DESeq2 analysis

    Arguments:
        input file: input file in the form of a .tsv or .csv"""
    df = get_df_from_file(input_file)

    # Get a list of all the conditions in the htseq-count file
    columns = df.columns
    groups = columns[2:].to_list()

    total_design_matrix = pd.DataFrame()
    group_column = []
    sample_column = []
    sample_name_column = []
    # split each condition using '_' as a delimiter
    for group in groups:
        group_list = group.split('_')
        sample_name_column.append(group)
        # Extract group number (index 0) and append to design_matrix group column
        group_column.append(group_list[0])
        # Extract sample number (index 1) and append to design_matrix sample column
        sample_column.append(group_list[1])

    # add group column and sample column to the design_matrix data frame
    total_design_matrix = total_design_matrix.assign(Sample_Name=sample_name_column)
    total_design_matrix = total_design_matrix.assign(Group_Number=group_column)
    total_design_matrix = total_design_matrix.assign(Sample_Number=sample_column)
    # Returns design_matrix
    return total_design_matrix


# Given an input of the sample condition and the control condition to be compared to, retrieve all replicates of
# the sample condition and all replicates of the control condition and output a sub_design_matrix
def get_sub_design_matrix(sample_condition, control_condition, total_design_matrix):
    """Sub-divide design matrix into smaller design matrices based on conditions for use in DESeq2 analysis

    Arguments:
        sample condition: input sample condition number that comes from individual samples from get_conditions()
        control condition: input control condition number that comes from individual samples from get_conditions()
        total design matrix: input design matrix that comes from get_total_design_matrix()"""
    # Condition parameters must be changed to strings to be compared to data in data frame
    sample_matrix = total_design_matrix.loc[total_design_matrix['Group_Number'] == str(sample_condition)]
    control_matrix = total_design_matrix.loc[total_design_matrix['Group_Number'] == str(control_condition)]
    # Concatenate the sample and control sub matrices together
    sub_design_matrix = pd.concat([sample_matrix, control_matrix], ignore_index=True)
    return sub_design_matrix


# Removes superfluous first index column and extra summary data at the end of the file
def get_total_count_data(input_file):
    """Given an input file return a total count matrix for use in DESeq2 analysis

    Arguments:
        input file: input file in the form of a .tsv or .csv """
    total_count_matrix = get_df_from_file(input_file)
    total_count_matrix = total_count_matrix.iloc[:-5,1:]
    # Returns count_matrix
    return total_count_matrix


# Given an input of the sample condition and the control condition to be compared to, retrieve all replicates of
# the sample condition and all replicates of the control condition and output a sub_count_matrix
def get_sub_count_matrix(sample_condition, control_condition, total_count_matrix):
    """Sub-divide count matrix into smaller count matrices based on conditions for use in DESeq2 analysis

    Arguments:
        sample condition: input sample condition number that comes from individual samples from get_conditions()
        control condition: input control condition number that comes from individual samples from get_conditions()
        total design matrix: input design matrix that comes from get_total_design_matrix()"""
    # Retrieve lists of column names that match the sample and control conditions
    columns = list(total_count_matrix.columns)
    sample_condition_columns = []
    control_condition_columns = []
    for column in columns:
        if column[0] == str(sample_condition):
            sample_condition_columns.append(column)
        elif column[0] == str(control_condition):
            control_condition_columns.append(column)
        else:
            pass

    # Combine gene ID, sample columns, and control columns (have to do multiple steps to allow concat of lists
    sub_matrix_columns = []
    sub_matrix_columns.append(columns[0])
    sub_matrix_columns += sample_condition_columns + control_condition_columns

    # Subset total_design_matrix by only the sample and control condition columns
    sub_count_matrix = total_count_matrix.loc[:, sub_matrix_columns]
    return sub_count_matrix


def get_gene_name(df, gene_id_column, scope):
    """This function takes in a data frame, a gene id column name, and a scope to find all gene names for the
    respective gene ids in the data frame. Has a lot of overlap with Get_Gene_Names but isn't as hardcoded.

    Arguments:
        df: data frame to iterate over
        gene_id_column: name of the column containing gene names
        scope: the scope that is to be searched (the database -- 'refseq', 'ensembl'"""
    # Initialize mygene
    mg = mygene.MyGeneInfo()

    # Create two lists to store all names and symbols
    gene_name_list = []
    gene_symbol_list = []

    gene_id_list = df[gene_id_column].tolist()

    # Iterate through all gene IDs to get gene names and symbols
    mg_list = mg.querymany(gene_id_list, scopes=scope, returnall=True)
    for mg in mg_list['out']:
        try:
            gene_name_list.append(mg['name'])
            gene_symbol_list.append(mg['symbol'])
        # If no name or symbol is found, insert - into the data frame
        except:
            gene_name_list.append(" - ")
            gene_symbol_list.append(" - ")
            # print("No value found.")

    # Add gene symbol and gene name to new data frame
    outputdf = pd.DataFrame()
    outputdf['Gene Symbol'] = gene_symbol_list
    outputdf['Gene Name'] = gene_name_list
    return outputdf


def write_df_to_file(df, output_file):
    """write a pandas data frame to an output file

        Arguments:
            df: output DESeq2 data frame to be written
            output file: file name"""
    if output_file[-4:] == ".tsv":
        df.to_csv(output_file, sep='\t')
    else:
        df.to_csv(output_file, sep=',')


# **************************** TEST PROCESS ****************************
