#!/usr/bin/env python3

from difflib import SequenceMatcher
from collections import Counter
import pathlib
import os
import argparse

#------------------------------------------------------------------------------#
# Main functions
#------------------------------------------------------------------------------#
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

def read_junction_file(junction_file_path):
    """
    Takes a tab-delimted junction file with the columns read_name, 5-end, and 3-end,
    and make each line a list of those three items. Ditches the column names.
    Returns a list of tuples - (read_name, 5-end, 3-end)
    """

    out_list = []

    with open(junction_file_path) as in_handle:
        for line in in_handle:
            line = line.rstrip("\n")

            if line.startswith("read_name"):
                continue

            line = line.split("\t")

            # Make the start end and positions integers
            line[1] = int(line[1])
            line[2] = int(line[2])

            # Append to output as a tuple
            out_list.append(tuple(line))

    return out_list

def find_seed_sequence(junction, n, ref):
    """
    Given a junction, which is a tuple or list of
    stucture: name, 5' location, 3' location

    Returns a tuple consisting of n/2 nucleotides on either
    side of each location.

    I.e. if n ==20, then this returns ten bases on either side
    of the junction. n must be even.
    """
    # Make sure it's even.
    if (n % 2) != 0:
        raise ValueError("n must be even. You input {}.".format(n))

    # calculate step
    step = n/2

    # Get n/2 bases on either side of each junction
    five_start = int(junction[1] - step)
    five_end = int(junction[1] + step)
    five_seq = ref[five_start:five_end]

    three_start = int(junction[2] - step)
    three_end = int(junction[2] + step)
    three_seq = ref[three_start:three_end]

    # Find longest continuous match between these sequences
    match = SequenceMatcher(None, five_seq, three_seq).find_longest_match(0, len(five_seq), 0, len(three_seq))
    five_seed_start = match[0] # position of match in the first sequence
    three_seed_start = match[1] # position of match in the second sequence
    seed_length = match[2] # length of match

    return five_seq[five_seed_start:five_seed_start+seed_length]

def assign_seeds(junction_list, positions, ref, n=30):
    """
    Junction list is a list of junctions, where each junction
    is of structure read_name, 5', 3'.

    Positions is a dictionary of structure
    position:(start, end), and is used to assign where the junction is.

    Junction position is determined by their 3' junction point, with the
    exception of ORF1a which is determined by its 5' junction point.

    To calcualte the seed, this function calls find_seed_sequence. What
    that function does is takes the n/2 nucleotides on either side of
    the start and end of the junction point (using the genome sequence),
    and finds the longest continuous 100% correct alignment between the
    start and end seqeunces.

    Therefore, this function assigns each seed based on where the junction
    is. Output is a dictionary of structure:
    position:[list of seeds]
    """

    seed_dict = dict()
    for junction in junction_list:

        five_prime = junction[1]
        three_prime = junction[2]

        # Find seed sequence of the junction
        seed = find_seed_sequence(junction, n, ref)

        # Determine which region it falls in
        position = ''

        # If the 5' junction is within ORF1a, will assign it to
        # ORF1a
        if five_prime in range(positions['ORF1a_internal'][0], positions['ORF1a_internal'][1]):
            position = "ORF1a_internal"

        # Otherwise, base it on the 3' coordiate. Because of this else here,
        # any ORF1a junctions that do end up somewhere else will still be
        # assigned as ORF1a.
        else:
            for location_name, coordinates in positions.items():
                coord_range = range(coordinates[0],coordinates[1])
                if three_prime in coord_range:
                    position = location_name

        # Append the seed to the seed_dict based on determined location
        if not position in seed_dict:
            seed_dict[position] = []
        seed_dict[position].append(seed)

    return seed_dict

def get_seed_percentages(seed_dict, top=2):
    """
    Takes the seed_dict, which is a dictionary of structure
    position: list of seed sequences

    And returns the top n (determined by the top input) seeds
    by percentage abundance. It then returns an output dictionary
    of the top seed's abundance, as well as 'other' percent abundance
    (such that percents for a given position sum to 100).

    output is a dictionary of structure:
    position: [(seed1, percentage_seed1), ..., (other, percentage_other)]
    """

    def find_other(in_list):
        """
        Given an input list containing tuples of structure
        some_string:number,

        returns 1 - sum(all numbers)
        """
        cumulative = 0
        for tup in in_list:
            cumulative += tup[1]

        return 100 - cumulative


    seed_percentage_dict = dict()

    for position, seed_list in seed_dict.items():
        c = Counter(seed_list)
        percentages = [(i, c[i] / len(seed_list) * 100.0) for i, count in c.most_common()]

        # Take top percentages
        out_percentages = percentages[:top]

        # Make entry "other" for the others
        out_percentages.append(("Other", find_other(out_percentages)))

        # append to output
        seed_percentage_dict[position] = out_percentages

    return seed_percentage_dict

def write_output(seed_percentages, output_path):
    """
    Writes out seed_percentages to a tab-delimited file with the
    columns:
        - ORF [i.e. ORF8]
        - position_type (cannonical or internal)
        - seed
        - percentage
    """

    # Format output
    output = 'ORF\tposition_type\tseed\tpercentage\trank\n'
    for position, percentage_tup_list in seed_percentages.items():

        # This is basically the end of the genome, and intergenic sequences.
        # Skip em.
        if position == "":
            continue

        ORF = position.split("_")[0]
        position_type = position.split("_")[1]
        rank = 0

        for tup in percentage_tup_list:

            rank += 1

            seed, percentage = tup[0], tup[1]

            if seed == "":
                seed = "NO_SEED"

            output += "{}\t{}\t{}\t{}\t{}\n".format(ORF,
                                                   position_type,
                                                   seed,
                                                   percentage,
                                                   rank)

     # Generate the output directory if necessary
    out_dir = os.path.dirname(output_path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as outfile:
        outfile.write(output)

#------------------------------------------------------------------------------#
# Constants
#------------------------------------------------------------------------------#
positions = {
    "N_cannonical":(28244,28280),
    "N_internal": (28280,29533),
    "ORF8_cannonical":(27872,27908),
    "ORF8_internal": (27908,28244),
    "ORF7A_cannonical":(27372,27408),
    "ORF7A_internal": (27408,27759),
    "M_cannonical":(26457,26493),
    "M_internal": (26523,27191),
    "E_cannonical":(26221,26257),
    "E_internal": (26257,26457),
    "ORF3A_cannonical":(25369,25405),
    "ORF3A_internal": (25405,26220),
    "S_cannonical":(21540,21576),
    "S_internal": (21576,25369),
    "ORF1a_internal": (266,13468)
}

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is detect homologous sequences at the 5'
            and 3' ends of each junction. It takes the 5' coordinate, identifies
            the genome sequence of 15bases on either side, and records that
            30 base sequence. It does the same for the 3' coordinate, resulting
            in two 30 base sequences. It then finds the longest continuous match
            between those two sequences, which is considered the "seed sequence",
            aka TRS aka homologous sequence.

            Then, for each category (explained below), this script determines
            what percentage of total seed sequences is each seed sequence. It
            outputs the top N seed sequences.

            Categories are detfined by a constant dictionary, 'positions', that
            defines coordinates. There are two types of positions. First, are
            positions called {ORF}_cannonical. These coordiantes include the
            canonical TRS site and 15 bases on either side. The other type of
            position is {ORF}_internal. These are cooridinates for an ORF, but
            these coordinates are superceded by any canonical coordiantes (i.e.
            if there is a TRS sequence for gene B that is located inside the
            5' gene A reading frame, the gene A 'internal' coordinates end
            early when it hits the 15 base buffer around gene B's TRS site).

            Which gene/category a junction falls into depends on it's 3'
            coordinate, with the exception of ORF1a. If junction has a 5'
            coordinate in the ORF1a it is assigned to ORF1a_internal.
            """)

    # Required arguments
    parser.add_argument(
        '-j',
        '--junction_file',
        type=str,
        required=True,
        help='''
        Path to file containing 5' and 3' junction coordinates.
        '''
    )
    parser.add_argument(
        '-o',
        '--output_file',
        type=str,
        required=True,
        help='''
        Path to the desired output file listing the top N TRS sequences for each
        category, as well as their percentages.
        '''
    )
    parser.add_argument(
        '-r',
        '--reference_fasta',
        type=str,
        required=True,
        help='''
        Path to the reference fasta.
        '''
    )

    # Optional arguments
    parser.add_argument(
        '-p',
        '--proximal_sequence',
        type=int,
        required=False,
        default=30,
        help='''
        Amount of sequence to consider on either side of each junction point.
        For example, an entry of 30 means it will look for 15 bases on either
        side of the 5' and 3' junctions, and compare those two 30 base sequences
        to find the homologous sequence. This must be even.
        '''
    )
    parser.add_argument(
        '-N',
        '--top',
        type=int,
        required=False,
        default=2,
        help='''
        The top N seed sequences for each category to output.
        '''
    )

    args = parser.parse_args()

    # Define input variables
    junction_file = args.junction_file
    output_file = args.output_file
    reference_fasta = args.reference_fasta
    proximal_sequence = args.proximal_sequence
    top = args.top

    # Read in reference fasta to a string
    ref = read_genome_fasta(reference_fasta)

    # read in junction file to a list
    junction_list = read_junction_file(junction_file)

    # Find seed sequences, and assign them to their position - cannonical/outside the gene
    # or inside the gene.
    seed_dict = assign_seeds(junction_list, positions, ref, n=proximal_sequence)

    # Find the percentages of each unique seed in the seed_dict, for each position
    seed_percentages = get_seed_percentages(seed_dict, top=top)

    # Write to output file
    write_output(seed_percentages, output_file)

if __name__ == '__main__':
    main()
