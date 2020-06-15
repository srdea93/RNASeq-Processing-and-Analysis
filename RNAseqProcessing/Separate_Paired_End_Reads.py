#!/usr/bin/env python
# Command line python script to separate all paired end fastq files in a directory into their own sub directories
# paired end reads MUST be in the format of: "filename_..._1.fastq" (for read 1) and "filename_..._2.fastq" (for read 2)

import argparse
import os
import fnmatch
import shutil
import itertools

def arg_parse():
    parser = argparse.ArgumentParser(description="Separate all paired end fastq files into subdirectories")
    parser.add_argument('i', metavar='Input Directory', nargs=1, help="Paired-end fastq file directory")
    return parser


# Finds all files matching a specific pattern within a directory and returns a list of strings
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def separate_paired_end_reads(input_dir):
    ext_name = ".fastq"
    # get a list of all the fastq files
    fastq_list1 = find("*1.fastq", input_dir)
    fastq_list2 = find("*2.fastq", input_dir)

    # if files are not named .fastq, check for .fq
    if not fastq_list1:
        ext_name = ".fq"
        fastq_list1 = find("*1.fq", input_dir)
        fastq_list2 = find("*2.fq", input_dir)

    # iterate through the lists and move files with same base name, only last character before .fastq is different
    if ext_name == ".fastq":
        for fastqr1 in fastq_list1:
            for fastqr2 in fastq_list2:
                if fastqr1[:-7] == fastqr2[:-7]:
                    # Make subdirectory
                    dir_name = fastqr1[:-8]
                    if not os.path.exists(dir_name):
                        os.makedirs(dir_name)
                    shutil.move(fastqr1, dir_name)
                    shutil.move(fastqr2, dir_name)
                    #break out of this iteration of for loop and start next iteration once match is found
                    break
                else:
                    pass
    # if file ext is .fq
    else:
        for fastqr1 in fastq_list1:
            for fastqr2 in fastq_list2:
                if fastqr1[:-4] == fastqr2[:-4]:
                    # Make subdirectory
                    dir_name = fastqr1[:-5]
                    if not os.path.exists(dir_name):
                        os.makedirs(dir_name)
                    shutil.move(fastqr1, dir_name)
                    shutil.move(fastqr2, dir_name)
                    #break out of this iteration of for loop and start next iteration once match is found
                    break
                else:
                    pass


def main():
    # Call parser
    parser = arg_parse()
    args = parser.parse_args()
    # Set up input list
    input = args.i[0]
    separate_paired_end_reads(input)


# **************************** RUN PROCESS ****************************

if __name__ == "__main__":
    main()


# **************************** TEST PROCESS ****************************