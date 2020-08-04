#!/usr/bin/env python3

from difflib import SequenceMatcher
from collections import Counter
import pathlib
import os
import argparse
import re
import random
import sys

def get_args():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to output the length distribution of
            homologous regions between 5' and 3' junction points.

            What this script does is find the n nucleotides flanking the 5'
            junction point, and the n nucleotides flanking the 3' junction point,
            and find the length of the longest homologous string between them.

            In addition, this script extracts all possible n-mers directly from
            the reference sequence, and determines the homology between random
            pairs of these reference-derived n-mers. This gives the random
            distribution of possible homology lengths. The number of random
            pairs, and the distance apart on the genome the member of each pair
            must be, can be specified.
            """)

    # Required
    parser.add_argument(
        '-j',
        '--junction_file_path',
        type=str,
        required=True,
        help='''
        Path to the junction file. File should be tab-delimited and of structure
        read_name, 5'-coordinate, 3'-coordinate.
        '''
    )
    parser.add_argument(
        '-g',
        '--genome',
        type=str,
        required=True,
        help='''
        Path to the reference genome.
        '''
    )
    parser.add_argument(
        '-o',
        '--output_length_distribution',
        type=str,
        required=True,
        help='''
        Path to the output directory containing the homology length distribution.
        '''
    )

    # Optional
    parser.add_argument(
        '-n',
        '--n',
        type=int,
        required=False,
        default=30,
        help='''
        How many nucleotides at each junction point to examin for homology. This
        must be even. If 30 is specified, for example, it will find 15 nucleotides
        on either side of the 5' junction point to be compared with the 15
        nucleotides on either side of the 3' junction point.
        '''
        )
    parser.add_argument(
        '-t',
        '--TRS',
        type=str,
        required=False,
        default="ACGAAC",
        help='''
        The TRS core sequence. This is used to identify if a junction falls near
        the core sequence. Currently, a junction is classified as being Canonical
        if it has a 5prime and 3prime end within 15 nucleotides of this TRS
        core sequence.
        '''
        )
    parser.add_argument(
        '-m',
        '--minimum_distance_apart',
        type=int,
        required=False,
        default=1000,
        help='''
        When calculating homology between random n-mers, this is the distance
        apart that the n-mers need to have originated from.
        '''
        )
    parser.add_argument(
        '-c',
        '--number_of_comparisons',
        type=int,
        required=False,
        default=100000,
        help='''
        When calculating homology between random n-mers, this is the number of
        random n-mer pairs to compare.
        '''
        )

    args = parser.parse_args()

    return args

#------------------------------------------------------------------------------#
# Functions
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

def define_canonical_trs_ranges(TRS, ref):
    """
    Find the locations of the TRS core sequence
    ACGAAC and returns a range of 10 nucleotides
    on either side. So, if a TRS is located
    from 21555 21561 this function reports back
    [(21545, 21571),...]
    """
    range_list = []
    for match in re.finditer(TRS, ref):
        start = match.start()
        end = match.end()

        range_start = start - 15
        range_end = end + 15

        range_list.append((range_start, range_end))

    return range_list

def find_seed_sequence(junction, n, ref):
    """
    Given a junction, which is a tuple or list of
    stucture: name, 5' location, 3' location
    Returns a tuple of (5' nucleotides, 3' nucleodies),
    where the nucleotides are the n/2 nucleotides on
    either side of the 5' and 3' junctions respectively.
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

    # Edge case (haha, get it?): If 5' end is less than 'step' nucleotides away from
    # front of ref sequence, need to start the five_seq at the start
    # of the reference sequence
    if junction[1] < step and five_seq == "":
        five_seq = ref[:junction[1] + int(step)]

    return (five_seq, three_seq)


def find_longest_match(seq1, seq2):
    """
    Given two sequences, finds the longest perfect continuous
    match between them
    """

    match = SequenceMatcher(None, seq1, seq2).find_longest_match(0, len(seq1), 0, len(seq2))
    seq1_seed_start = match[0] # position of match in the first sequence
    seq2_seed_start = match[1] # position of match in the second sequence
    seed_length = match[2] # length of match

    return seq1[seq1_seed_start:seq1_seed_start+seed_length]

def assess_junction_homology(junctions, n, ref):
    """
    Takes in the junctions list of tuples and, for each
    junction, determines 1) The seed sequence at the 5' end,
    2) the seed sequence at the 3' end, 3) the longest homologous
    sequence beteween the seed sequences, and 4) the length of that
    sequence.

    The 'seed sequences' are simply the n/2 bases on either side of
    each junction point.

    Output is a list of tuples, where each tuple contains:
    1) read name
    2) 5' junction location
    3) 3' junction location
    4) 5' seed sequence (n/2 nucleotides on either side of 5' junction location)
    5) 3' seed sequence (n/2 nucleotides on either side of 3' junction location)
    6) Longest homologous match between 5' and 3' seeds
    7) Length of that longest homologous match
    """

    output_list = []
    for junction in junctions:

        read_name, five, three = junction

        # Find the seed (flanking) sequencing
        five_seed, three_seed = find_seed_sequence(junction, n, ref)

        # Find homology between the flanking sequences
        homologous_sequence = find_longest_match(five_seed, three_seed)

        # Write to new output list of tuples
        output_list.append(
            (read_name,
             five,
             three,
             five_seed,
             three_seed,
             homologous_sequence,
             len(homologous_sequence)
            )
        )

    return output_list

def is_in_trs_range(pos, trs_ranges):
    """
    trs_ranges is a list of tupes, structure
    [(range1_start, range1_end), ...]

    Checks if pos is in one of the ranges in
    that list.
    """

    for trs_range in trs_ranges:
        r = range(trs_range[0], trs_range[1])
        if pos in r:
            return True
    return False

def split_canonical_and_noncanonical_homology_lengths(junctions, trs_ranges):
    """
    For each junction, determines if the 5' and 3' locations are
    within range of a TRS core sequence. If yes, then adds the
    homology length of that junction to the canonicals list, otherwise
    to the noncanonicals list.

    Output is a tuple of
    (canonical_junction_homology_list, noncanonical_junction_homology_list)
    """
    canonicals = []
    noncanonicals = []
    for junction in junctions:
        five = junction[1]
        three = junction[2]
        length = junction[6]

        if is_in_trs_range(five, trs_ranges) and is_in_trs_range(three, trs_ranges):
            canonicals.append(length)
        else:
            noncanonicals.append(length)

    return(canonicals, noncanonicals)

# For random n-mer homology calculations
#------------------------------------------------------------------------------#
def find_nmers(seq, n):
    """
    Reports out all nmers of the input seq, as well as their
    start locations. Output is

    [(nmer1, location1), (nmer2, location2), ... ]
    """
    nmers = []
    for i in range(len(seq)-30):

        nmer = seq[i:i+n]
        nmers.append((nmer, i))

    return nmers

def random_homology_comparisons(nmers, number_of_comparisons=100000, minimum_distance_apart=1000):
    """
    Given possible nmers, and thier locations, determines the length of the longest
    homologous sequence between a specified number of pairs of nmers. The analysis
    only occurs if the nmers are a sufficient distance appart.

    nmers is of structure:
    [(nmer1, location1), (nmer2, location2), ..]

    This function returns a list of homologous sequence lengths
    """

    i = 0
    match_list = []
    while i <= number_of_comparisons:

        nmer1 = random.choice(nmers)
        nmer2 = random.choice(nmers)

        # Make sure sufficiently far apart
        if not abs(nmer1[1] - nmer2[1]) >= minimum_distance_apart:
            continue

        # Find homology
        match = find_longest_match(nmer1[0], nmer2[0])
        match_list.append(len(match))
        i +=1

    return match_list


#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():

    # Get input arguments
    args = get_args()
    junction_file_path = args.junction_file_path
    genome = args.genome
    output_length_distribution = args.output_length_distribution
    n = args.n
    TRS = args.TRS
    minimum_distance_apart = args.minimum_distance_apart
    number_of_comparisons = args.number_of_comparisons

    # Main!
    #--------------------------------------------------------------------------#
    # Read in genome
    print("{}: Reading in reference genome.".format(sys.argv[0]))
    ref = read_genome_fasta(genome)

    # Find TRS core sequences, and 15nt on either side
    print("{}: Finding core TRS sequences.".format(sys.argv[0]))
    trs_ranges = define_canonical_trs_ranges(TRS, ref)

    # Read in junction file
    print("{}: Reading in junctions.".format(sys.argv[0]))
    junctions = read_junction_file(junction_file_path)

    # For each junction, find the 5' and 3' flanking sequences,
    # the longest homologous sequence, and the length of that
    # homologous sequence.
    print("{}: Finding homology between junction sites.".format(sys.argv[0]))
    junctions = assess_junction_homology(junctions, n, ref)

    # Generate separate lists of canonical and non-canonical
    # homology lengths
    print("{}: Splitting canonical and noncanonical.".format(sys.argv[0]))
    canonical_lengths, noncanonical_lengths = split_canonical_and_noncanonical_homology_lengths(junctions, trs_ranges)

    #--------------------------------------------------------------------------#
    # Find the homology length distribution of random n-mers from SARS2 genome
    # that are at least 1000bp away
    #--------------------------------------------------------------------------#
    print("{}: Assessment of random n-mer homology.".format(sys.argv[0]))

    # 1) Find all n-mers in SARS2 genome, along with their location
    # 2) Determine homology length between a specified number of random
    #    n-mer pairs, but the n-mer pairs must be >1000bp away from one
    #    another.

    nmers = find_nmers(ref, n)

    random_match_lengths = random_homology_comparisons(nmers, number_of_comparisons=number_of_comparisons, minimum_distance_apart=minimum_distance_apart)

    # Now, output the canonical_lengths, noncanonical_lengths, and random_match_lengths
    # to a tidy tsv w/ two columns : type (canonical/non-canonical-random), and length
    output = "type\tlengths\n"
    for match in canonical_lengths:
        entry = "canonical\t{}\n".format(match)
        output+=entry
    for match in noncanonical_lengths:
        entry = "noncanonical\t{}\n".format(match)
        output+=entry
    for match in random_match_lengths:
        entry = "random\t{}\n".format(match)
        output+=entry

    # Make output dir's if needed, and write output
    print("{}: Writing output.".format(sys.argv[0]))
    out_dir = os.path.dirname(output_length_distribution)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    with open(output_length_distribution, "w") as outfile:
        outfile.write(output)

    print("{}: Done!".format(sys.argv[0]))

if __name__ == '__main__':
    main()
