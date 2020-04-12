#!/usr/bin/env python3

import pandas as pd
import itertools
import pathlib
import sys
import os
import argparse

def read_clstr(cluster_file):
    """
    Parse through a cd-hit cluster file and create a dictionary, where the key
    is the cluster title and the values are lists of headers present in the cluster.
    Adapted from https://github.com/Y-Lammers/CD-HIT-Filter/blob/master/CD-HIT-Filter.py
    """

    # parse through the .clstr file and create a dictionary
    # with the sequences per cluster

    # open the cluster file and set the output dictionary
    cluster_file, cluster_dic = open(cluster_file), {}

    # parse through the cluster file and store the cluster name + sequences in the dictionary
    cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
    for cluster in cluster_groups:
        name = next(cluster).strip()
        seqs = [seq.split('>')[1].split('...')[0] for seq in next(cluster_groups)]
        cluster_dic[name] = seqs

    # return the cluster dictionary
    return cluster_dic

def cluster_dict_to_df(cluster_dict):
    """
    Takes in the cluster_dict, which is a dictionary
    of structure -
    cluster_#:[seq_name1, seq_name2,...]

    And each seq_name is structure
    TXNAMExxxSTART_END_STRAND


    This function will make a dataframe with the following
    columns:
    cluster, transcript_name, start, end
    """

    # Will collect each row as a list, so the future dataframe
    # contents will be a list of lists.
    info = []

    for cluster_name, orf_list in cluster_dict.items():

        cluster_number = int(cluster_name.split(" ")[1])

        # Add each ORF into the output dataframe
        for orf in orf_list:

            orf = orf.split("xxx")

            tx_name = orf[0]
            start = orf[1].split("_")[0]
            end = orf[1].split("_")[1]

            # Append as a row
            row = [cluster_number, tx_name, int(start), int(end)]
            info.append(row)

    # Now generate the actual dataframe from the list of lists,
    # and specify the column names
    out_df = pd.DataFrame(info, columns=['cluster_number', 'tx_name', 'start', 'end'])

    return out_df

def main():

    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to extract the start and end
            coodinates, on the transcript, of each ORF. So, it takes in a cdhit
            cluster file and extracts the start and end site for each individual
            ORF, and labels it by its cluster. Thus, the coordinates are
            relative to each individual transcript, NOT to the genome.

            The output .tsv file contains the following columns:
            cluster_number, tx_name (transcript name), start, end.
            """)

    # Required arguments
    parser.add_argument(
        '-c',
        '--cluster_file',
        type=str,
        required=True,
        help='''
        Path to the cd-hit cluster file.
        '''
    )
    parser.add_argument(
        '-o',
        '--output_file',
        type=str,
        required=True,
        help='''
        Path to the output tsv file. This file will contain information on the
        start, end, transcript, and cluster of each ORF.
        '''
    )

    args = parser.parse_args()

    # Define input variables
    cluster_file = args.cluster_file
    output_file = args.output_file

    # Read in cluster file
    print("{}: Starting script.".format(sys.argv[0]))
    print("{}: Reading in cluster file.".format(sys.argv[0]))
    cluster_dict = read_clstr(cluster_file)

    # Generate dataframe with start and stop coords for each
    # transcript and each cluster
    print("{}: Extracting ORF information to a dataframe.".format(sys.argv[0]))
    coords_df = cluster_dict_to_df(cluster_dict)

    # Save the output df...
    print("{}: Saving output tsv file.".format(sys.argv[0]))
    out_dir = os.path.dirname(output_file)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    coords_df.to_csv(output_file, sep='\t', index=False)
    print("{}: Finished.".format(sys.argv[0]))

if __name__ == '__main__':
    main()
