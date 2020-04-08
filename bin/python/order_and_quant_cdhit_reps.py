#!/usr/bin/env python3

from Bio import SeqIO
import itertools
import os
import pathlib
import argparse
import sys

def read_fasta_to_memory(input_fasta):
    """
    Reads fasta into a memory as a dictionary with header:sequence
    """
    fasta_dict = dict()
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        fasta_dict[seq_record.id] = str(seq_record.seq)
    return fasta_dict


def read_clstr(cluster_file, rep_list=set(), reps_only = False):
    """
    Parse through a cd-hit cluster file and create a dictionary, where the key
    is the cluster title and the values are lists of headers present in the cluster.

    This has the option to only report out the representatives in the
    rep_list. If you want all members of the cluster, just enter the cluster_file.

    Adapted from https://github.com/Y-Lammers/CD-HIT-Filter/blob/master/CD-HIT-Filter.py
    """

    # parse through the .clstr file and create a dictionary
    # with the sequences per cluster

    # open the cluster file and set the output dictionary
    cluster_file, cluster_dic = open(cluster_file), {}

    # make sure rep_list is a set
    rep_list = set(rep_list)

    # parse through the cluster file and store the cluster name + sequences in the dictionary
    cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
    if reps_only == True:
        for cluster in cluster_groups:
            name = next(cluster).strip()
            seqs = [seq.split('>')[1].split('...')[0] for seq in next(cluster_groups)]
            seqs = [seq for seq in seqs if seq in rep_list]
            cluster_dic[name] = seqs[0]
    else:
        for cluster in cluster_groups:
            name = next(cluster).strip()
            seqs = [seq.split('>')[1].split('...')[0] for seq in next(cluster_groups)]
            cluster_dic[name] = set(seqs)

    # return the cluster dictionary
    return cluster_dic

def count_number_in_cluster(cluster_to_header):
    """
    Takes in a dictionary of key:list, and returns
    a dictionary of key:len(list)
    """

    cluster_name_to_cluster_count = dict()
    for cluster, header_list in cluster_to_header.items():
        cluster_name_to_cluster_count[cluster] = len(header_list)

    return cluster_name_to_cluster_count

def get_cluster_name(header, cluster_to_header):
    """
    For a given header, determine which cluster it is in.
    """
    for cluster, header_list in cluster_to_header.items():
        if header in header_list:
            return cluster

def format_headers_and_make_fasta(cluster_to_header_reps_only, header_to_seq, cluster_name_to_cluster_count):
    """
    This function reports out the seqs in fasta format. It iterates over
    cluster_to_header_reps_only which is structured in order from cluster 0
    and upwards. This function also puts the cluster count at the end of each header.

    The output is a string that will be fasta format when written to an output file.
    """

    output_contents = ""
    for cluster, header in cluster_to_header_reps_only.items():
        seq = header_to_seq[header]
        count = cluster_name_to_cluster_count[cluster]
        formated_cluster = cluster[1:].replace(" ", "_")

        fasta_formatted = ">{}xxx{}xxx{}\n{}\n".format(header,formated_cluster, count, seq)
        output_contents += fasta_formatted


    return output_contents

def write_output(output, output_path):

    # Generate the output directory if necessary
    out_dir = os.path.dirname(output_path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as output_handle:
        output_handle.write(output)

def main():
    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is order the CD-hit representatives in
            order by cluster, and to append to the representative's header the
            number of sequences in each cluster.
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
        '-f',
        '--fasta',
        type=str,
        required=True,
        help='''
        Path to the input fasta. This should have all of the sequences clustered
        in the cd-hit cluster file. But, this can be nucleotide or protein...
        only thing that matters is the sequence headers are the same.
        '''
    )
    parser.add_argument(
        '-o',
        '--output_representative_fasta',
        type=str,
        required=True,
        help='''
        Path to the output fasta containing the ordered and labeled cluster
        representatives.
        '''
    )

    args = parser.parse_args()

    # Define input variables
    cluster_file = args.cluster_file
    fasta = args.fasta
    output_representative_fasta = args.output_representative_fasta

    #--------------------------------------------------------------------------#
    # Main
    #--------------------------------------------------------------------------#
    # Parse fasta and prodigal files
    print("{}: Starting.".format(sys.argv[0]))
    print("{}: Reading in fasta to memory.".format(sys.argv[0]))
    header_to_seq = read_fasta_to_memory(fasta)

    print("{}: Reading in cluster file to memory.".format(sys.argv[0]))
    cluster_to_header_reps_only = read_clstr(cluster_file, header_to_seq.keys(), reps_only = True)
    cluster_to_header_all = read_clstr(cluster_file)

    # Get count for each cluster
    print("{}: Getting cluster counts.".format(sys.argv[0]))
    cluster_name_to_cluster_count = count_number_in_cluster(cluster_to_header_all)

    # Generate the output fasta contents
    print("{}: Generating and writing output.".format(sys.argv[0]))
    fasta_contents = format_headers_and_make_fasta(cluster_to_header_reps_only,
                                                   header_to_seq,
                                                   cluster_name_to_cluster_count)

    # Write output
    write_output(fasta_contents, output_representative_fasta)
    print("{}: Finished.".format(sys.argv[0]))

if __name__ == '__main__':
    main()
