#!/usr/bin/env python3

import pysam
import argparse
import pathlib
import os
import sys
import pandas as pd
from io import StringIO

#==============================================================================#
# Auxiliary functions
#==============================================================================#
def convert_cigar(cigar, cigar_key):
    result = []

    for tup in cigar:
        symbol = cigar_key[tup[0]]
        result.append((symbol, tup[1]))

    return result

def get_end_coordinate(cigar, start):

    current_pos = start

    for tup in cigar:

        if tup[0] == "I":
            continue
        elif tup[0] in {"M", "N", "D", "X", "="}:
            current_pos += tup[1]

    return current_pos

#==============================================================================#
# Main functions
#==============================================================================#
def read_genome_fasta(path):
    """
    Reads in fasta containing a single genome. It must therefore
    only have one header.
    """

    header_count = 0
    seq = ''
    with open(path) as infile_handle:
        for line in infile_handle:

            if line.startswith(">"):
                header_count += 1
                continue

            seq += line.rstrip("\n")

    if header_count > 1:
        raise ValueError("A single fasta entry is mandatory, this is a multifasta.")

    return seq

def get_junctions(cigar, start, min_intron_length=100):
    """
    This function reads in a cigar sequence and outputs all junctions
    (maked with N (deletion) or N (intron)) that are at least the
    required length. The junction position is based on the start
    position, meaning it is syncronized with the genome.

    The output is a list of tuples, where each tuple is
    (junction_start, junction_end).

    The resultant junctions are 0-indexed. Therefore, if a junction is
    i.e. (11, 30) of a sequence that starts at 1 and ends at 55,
    the resultant sequence ignoring the junction is simply

    seq[1:11] + seq[30:55]
    """


    junctions = []
    pos = start

    for entry in cigar:

        sign, length = entry

        # ignore softclipping and insertions
        if sign in {"S", "I"}:
            continue

        # matches, = (equals), and X (mismatch) advance the position
        elif sign in {"M", "=", "X"}:
            pos += length
            continue

        # if deletion or intron, action depends if it passes
        # the length requirement
        elif sign in {"D", "N"}:

            # If under the length requirement, don't consider it a
            # junction but do advance position
            if length < min_intron_length:
                pos += length
                continue

            # Otherwise, we have a junction
            elif length >= min_intron_length:

                # Add junction to output
                junc = (pos, pos+length) # this is 0-indexed
                junctions.append(junc)

                # still need to advance the position for the next item(s)
                pos += length

    return junctions

def parse_bam(bam_path, cigar_key, min_intron_length=100):
    """
    The purpose of this function is to parse in a bam
    and report out the start, junctions, and end of each
    alignment. Output is a dictionary of structure

    read_name:(start, junctions, end),
    where junctions is a list of structure:
    [(junc1_start, junc1_end), (junc2_start, junc2_end)...]
    """

    alignment_dict = dict()

    # parse bam
    for read in pysam.AlignmentFile(bam_path, 'rb'):

        if read.is_unmapped:
            continue

        # convert the cigar tuple list to their single-letter
        # designations.
        cigar = convert_cigar(read.cigar, cigar_key)

        # Get start and end coordinates
        start = read.reference_start
        end = get_end_coordinate(cigar, start)

        # Find junctions from cigar and start position
        junctions = get_junctions(cigar, start, min_intron_length)

        # Add information from this read to the output
        alignment_dict[read.query_name] = (start, junctions, end)


    return alignment_dict

def find_seq(genome, start, junctions, end):
    """
    This function takes in the start, junction, and ends
    of the desired sequence and extracts the non-junction
    sequence from the input genome sequence.
    """

    def get_regions(start, junctions, end):
        """
        This sub-function converts the start/junctions/end
        to a list of regions where there is sequence.
        """

        regions = []
        prev_junc_end = 0

        for i, junc in enumerate(junctions):
            junc_start, junc_end = junc

            # Handle first junction
            if i == 0:
                region = (start, junc_start)
                regions.append(region)
                prev_junc_end = junc_end

            # handle middle junctions
            if i > 0:
                region = (prev_junc_end, junc_start)
                regions.append(region)
                prev_junc_end = junc_end

            # handle the last junction
            if i == len(junctions) - 1:
                final_region = (junc_end, end)
                regions.append(final_region)

        return regions

    # if there are no junctions, just return start:end
    if junctions == []:
        return genome[start:end]

    # Get the regions where there IS sequence, aka
    # on either side of the junctions...
    regions = get_regions(start, junctions, end)

    # get the actual sequence from those coordinates
    seq = ''
    for i, region in enumerate(regions):

        # Need to add +1 to the region end, except not to the last last end
        # This fixes some really annoying indexing issue.
        if i != len(regions) -1:
            seq += genome[region[0]:region[1]]
        else:
            seq += genome[region[0]:region[1]]

    return seq

def gen_sequence_dict(genome, alignment_dict):
    """
    Takes in the alignment_dict, which is of structure
    read_header:(start, [(junc1_start, junc1_end)...], end)

    And returns a dictionary of structure
    read_header:sequence, where the sequence is from the coordinates
    for betweeen start and end but missing those sequences present
    in the junctions.
    """
    sequence_dict = dict()

    for header, coords in alignment_dict.items():

        start, junctions, end = coords

        seq = find_seq(genome, start, junctions, end)

        sequence_dict[header] = seq

    return sequence_dict


def dict_to_fasta(in_dict, path):
    """
    in_dict is dictionary of header:sequence,
    which will be written to the output path.
    """

    # Generate the output directory if necessary
    out_dir = os.path.dirname(path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Generate the output in fasta format
    fasta = ''
    for header, seq in in_dict.items():
        entry = ">{}\n{}\n".format(header, seq)
        fasta += entry

    # Write to output
    with open(path, "w") as outfile:
        outfile.write(fasta)

def gen_junction_df(alignment_dict, output_read_coordinates=False):
    """
    The alignment dict contains coordinates for the corrected
    read start, junctions, and end. This function extracts
    the junctions to a pandas dataframe with the columns
    'read_name', '5-end', '3-end'

    If output_read_coordinates == True, the output columns are
    read_name, 5-end, 3-end, read_start, read_end
    """

    # Initiate the output
    if output_read_coordinates == False:
        output = ["read_name", "5-end", '3-end']
    else:
        output = ["read_name", "5-end", '3-end', "read_start", "read_end"]

    output = '\t'.join(output) + '\n'

    for header, coords in alignment_dict.items():

        start, juncs, end = coords
        for junc in juncs:

            junc_start, junc_end = junc

            if output_read_coordinates == False:
                output += "{}\t{}\t{}\n".format(header, junc_start, junc_end)
            else:
                output += "{}\t{}\t{}\t{}\t{}\n".format(header, junc_start, junc_end, start, end)

    # return output as a pandas dataframe
    return pd.read_csv(StringIO(output), sep="\t")

def save_junction_df(junction_df, junction_df_path):

     # Generate the output directory if necessary
    out_dir = os.path.dirname(junction_df_path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    junction_df.to_csv(junction_df_path, sep='\t', index=False)

#==============================================================================#
# Main
#==============================================================================#
def main():
    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to generate synthetic transcripts from
            aligned direct-RNAseq data in bam format. Essentially, this script
            records the start and end of a given read's alignment, as well as
            any gaps/introns in the alignment that are at least a given size.
            The genome at these sites and in between these gaps is then
            extracted to act as synthetic reads. The rational is this approach
            will capture large-scale mRNA rearrangements but will correct
            nanopore read errors.
            """)

    # Required arguments
    parser.add_argument(
        '-g',
        '--genome_path',
        type=str,
        required=True,
        help='''
        Path to the genome that was mapped against, in fasta format.
        '''
    )

    parser.add_argument(
        '-b',
        '--bam_path',
        type=str,
        required=True,
        help='''
        Path to the input bam containig reads that were mapped to the specified
        genome. These are expected to be nanopore direct-RNAseq reads.
        '''
    )

    parser.add_argument(
        '-o',
        '--output_fasta',
        type=str,
        required=True,
        help='''
        Path to the resultant synthetic reads, which will be in fasta format.
        '''
    )

    parser.add_argument(
        '-j',
        '--junction_path',
        type=str,
        required=True,
        help='''
        Path to the a tab-delimited file containing information on each read's
        junctions. The columns are: read_name, 5'-end, 3'-end.
        '''
    )

    # Optional
    parser.add_argument(
        '-l',
        '--min_intron_length',
        type=int,
        required=False,
        default=1000,
        help='''
        Minimum intron (CIGAR labeled as N) or deletion (CIGAR labeled as D)
        that will be considered a junction.
        '''
    )
    parser.add_argument(
        '-r',
        '--output_read_coordinates',
        type=bool,
        required=False,
        default=False,
        help='''
        If False, the junction file has columns: read_name, 5'-end, 3'-end.
        If True, the columns are: read_name, 5'-end, 3'-end, read_start, read_end
        '''
    )

    args = parser.parse_args()

    # Define input variables
    genome_path = args.genome_path
    bam_path = args.bam_path
    output_fasta = args.output_fasta
    junction_path = args.junction_path
    min_intron_length = args.min_intron_length
    output_read_coordinates = args.output_read_coordinates

    # Constants
    #--------------------------------------------------------------------------#
    cigar_key = {
        0: "M",
        1: "I",
        2: "D",
        3: "N",
        4: "S",
        5: "H",
        6: "P",
        7: "=",
        8: "X",
        9: "B"
    }

    # Main
    #--------------------------------------------------------------------------#

    # Read in genome
    print("{}: Starting.".format(sys.argv[0]))
    print("{}: Reading in genome fasta.".format(sys.argv[0]))
    genome = read_genome_fasta(genome_path)

    # Parse bam
    print("{}: Reading in bam.".format(sys.argv[0]))
    alignment_dict = parse_bam(bam_path, cigar_key, min_intron_length)

    # Translates the start, junction, and end coordinates to sequences
    print("{}: Generating synthetic transcripts.".format(sys.argv[0]))
    sequence_dict = gen_sequence_dict(genome, alignment_dict)

    # Write out synthetic sequences as fasta
    print("{}: Writing out synthetic transcripts.".format(sys.argv[0]))
    dict_to_fasta(sequence_dict, output_fasta)

    # Generate junction df and save to output
    print("{}: Writing out junctions.".format(sys.argv[0]))
    junction_df = gen_junction_df(alignment_dict, output_read_coordinates)
    save_junction_df(junction_df, junction_path)

    print("{}: Finished..".format(sys.argv[0]))

if __name__ == '__main__':
    main()
