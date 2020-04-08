#!/usr/bin/env python3

from Bio import SeqIO
import itertools
import os
import pathlib
import argparse

#==============================================================================#
# Define functions
#==============================================================================#

def read_fasta_to_memory(input_fasta):
    """
    Reads fasta into a memory as a dictionary with header:sequence
    """
    fasta_dict = dict()
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        fasta_dict[seq_record.id] = str(seq_record.seq)
    return fasta_dict


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

def remove_clusters_with_reference_sequences(input_dict):
    """
    Input dictionary:
    - keys are cluster names
    - values are a list of fasta headers in said cluster.

    This function removes any cluster that contains a fasta header that
    starts with "REF" and returns the modified dictionary.
    """

    # Go through cluster list and remove clusters that have a reference file
    output_dict = dict()

    # Iterate through input
    for cluster_title, sequences in input_dict.items():

        # Copy all of them to the new dict
        output_dict[cluster_title] = sequences

        # But go back and remove clusters that have reference sequences
        for sequence in sequences:
            if sequence.startswith("REF"):
                del output_dict[cluster_title]
                break

    return output_dict

def write_to_fasta(cluster_dict, output_prefix, fasta_dict, min_number_of_sequences_in_cluster=True):
    # Generate the output directory if necessary
    out_dir = os.path.dirname(output_prefix)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Write clusters to output files
    for cluster_title, sequences in cluster_dict.items():

        if len(sequences) < min_number_of_sequences_in_cluster:
            continue

        # Clean up the cluster title
        cluster_title = cluster_title.strip(">").replace(" ", "_")

        # Define output path
        output_path = output_prefix + "_" + cluster_title + ".fasta"

        # Open output file and start writing sequences
        with open(output_path, "w") as output_handle:
            for sequence_header in sequences:
                seq = fasta_dict[sequence_header]
                output_handle.write(">" + sequence_header + "\n" + seq + "\n")

#==============================================================================#
# Main
#==============================================================================#

def main():
    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to take a cd-hit clstr file, as well
            as a mutlifasta, and split the multifasta into separate files based
            on the cd-hit clusters. Each output file will be named by the cd-hit
            cluster number.
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
        '--output_prefix',
        type=str,
        required=True,
        help='''
        Prefix, specifies where the output files will be saved and gives them
        a file prefix. Outfiles will be {PREFIX}_Cluster_NN.fasta
        '''
    )

    parser.add_argument(
        '-m',
        '--min_number_of_sequences_in_cluster',
        type=int,
        required=False,
        default=2,
        help='''
        Minimum number of sequences that need to be in a cluster for the cluster
        to be reported out.
        '''
    )

    args = parser.parse_args()

    # Define input variables
    cluster_file = args.cluster_file
    fasta = args.fasta
    output_prefix = args.output_prefix
    min_number_of_sequences_in_cluster = args.min_number_of_sequences_in_cluster

    #--------------------------------------------------------------------------#
    # Execute main
    #--------------------------------------------------------------------------#
    # Store the fasta contents into memory as a dictionary
    print("Reading fasta to memory.")
    fasta_dict = read_fasta_to_memory(fasta)

    # Read in cluster file
    print("Reading in cluster file.")
    cluster_dict = read_clstr(cluster_file)

    # Remove clusters with reference sequences
    print("Removing clusters that contain reference sequences.")
    cluster_dict = remove_clusters_with_reference_sequences(cluster_dict)

    # Write those to a fasta
    print("Writing clusters to separate fasta files.")
    write_to_fasta(cluster_dict, output_prefix, fasta_dict, min_number_of_sequences_in_cluster)

    print("Script finished.")

if __name__ == '__main__':
    main()
