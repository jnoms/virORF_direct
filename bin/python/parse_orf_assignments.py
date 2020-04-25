#!/usr/bin/env python3

import argparse
import pathlib
import os
import sys
from collections import Counter
from io import StringIO
import pandas as pd


#==============================================================================#
# Auxiliary Functions
#==============================================================================#
def validate_fields(diamond_fields, required_fields):
    for field in required_fields:
        if not field in diamond_fields:
            raise ValueError("Missing a required diamond field, and possibly others: {}".format(field))

def get_most_fiveprime(orf_list):
    """
    Given a list of ORFs of structure START_END_STRAND,
    it returns the ORF that is closer to the 5' end (i.e.
    has a smaller START).
    """

    most_fiveprime = ''

    for orf in orf_list:

        start = int(orf.split("_")[0])
        if most_fiveprime == '':
            most_fiveprime = orf
        else:
            prev_start = int(most_fiveprime.split("_")[0])
            if start < prev_start:
                most_fiveprime = orf

    return most_fiveprime

def get_longest_match(match_list, diamond_fields):
    """
    match_list is a list of lists, where each sublist is
    the contents of a diamond line. This function returns
    the match (entire sublist) that has the longest alignment
    length.
    """

    longest = ""

    for match in match_list:
        if longest == "":
            longest = match
            continue

        prev_length = int(longest[diamond_fields.index("length")])
        curr_length = int(match[diamond_fields.index("length")])

        if curr_length > prev_length:
            longest = match

    return longest

def is_cannonical(match, diamond_fields, pr_length_dict, min_ok_nterm_extension=20):
    """
    This function takes in a match, consisting of diamond_fields,
    and checks if qstart == sstart and qend == send. If these do match,
    thus function returns True. If not, it returns False. It also
    makes sure the length is the same as in the pr_length_dict.

    This function takes in a match, consisting of diamond_fields,
    and checks:
    1) If the alignment length is the expected length reported in the
    pr_length_dict
    2) If the qstart==sstart and qend==send.

    If either is not true, this function returns  "False"

    !! Notable exception/ASSUMPTION !!
    If there is an N-terminal extension of 20AA or less, but the
    stop position is expected, will return TRUE. The reason is
    Prodigal is set to call from the first nTG. There are sometimes
    alternative upstream nTGs which are called, meaning they look
    noncannonical even though they are not.
    """

    subject = match[diamond_fields.index("sseqid")]

    # Need to check for an N-terminal extension that may suggest this is
    # simply translated from an early start codon. If this is the case,
    # the subject end should be the expected length and the query_start
    # should be under some small/allowable threshold.
    if subject in pr_length_dict:
        expected_len = pr_length_dict[subject]

        # Check that the subject alignment is as expected
        if int(match[diamond_fields.index("sstart")]) == 1 and int(match[diamond_fields.index("send")]) == expected_len:

            # Now, check if qstart is below the allowable threshold
            if int(match[diamond_fields.index("qstart")]) < min_ok_nterm_extension:
                return True

    # Check length (only if the ORF is in the pr_length_dict)
    if subject in pr_length_dict:
        expected_len = pr_length_dict[subject]
        if expected_len != int(match[diamond_fields.index("length")]):
            return False

    qstart = match[diamond_fields.index("qstart")]
    qend = match[diamond_fields.index("qend")]
    sstart = match[diamond_fields.index("sstart")]
    ssend = match[diamond_fields.index("send")]

    if qstart != sstart or qend != ssend:
        return False
    else:
        return True

def dict_value_counts(in_dict):
    """
    This function takes in a dictionary, and returns another dictionary of
    structure key:length of the input value. For example, if the in_dict
    values are a list, this returns the key:length of each list as a new
    dictionary.
    """

    out_dict = dict()

    for key, value in in_dict.items():
        out_dict[key] = len(value)

    return out_dict

def aggregate_cannonical_and_varaint_split_counts(cannonical_and_variant_split_counts):
    out_counts = dict()
    for ORF, count in cannonical_and_variant_split_counts.items():

        if ORF.endswith("_variant"):
            ORF = ORF.split("_variant")[0]

        out_counts[ORF] = out_counts.get(ORF, 0) + count

    return out_counts

#==============================================================================#
# Main Functions
#==============================================================================#
def parse_diamond(diamond_infile_path, diamond_fields, min_pident):
    """
    This function reads in a diamond infile. The query (qseqid) must
    be formatted appropriately, i.e. [TRANSCRIPT_NAME]xxx[ORF]
    but obviously without the brakets. ORF is of structure
    START_END_STRAND.

    This function make a dict within a dict, where the first dictinoary
    has keys that are each TRANSCRIPT_NAME, and the inner-dicts are based
    on the ORF. The value for each ORF is a list of each line it matched to.
    If the same orf from the same transcript matched to multiple things,
    this match list will have more than one inner list. The MATCH_LIST is thus
    a list of lists... Most ORFs will have only one match, so MATCH_LIST will
    usually be a list that contains one list.

    diamond_dict{TRANSCRIPT_NAME{ORF: MATCH_LIST}}
    MATCH_LIST = [[MATCH1_qseqid, MATCH1_sseqID, MATCH1_pident, ....], [MATCH2_qseqid, ...], ...]

    Thus, diamond_dict[TRANSCRIPT_NAME1] will give a subdict containing all
    ORF:MATCH_LIST for that TRANSCRIPT_NAME1. In turn, diamond_dict[TRANSCRIPT_NAME1][ORF1]
    will simply yeild a list of matches (where each match is a list from a
    line in the diamond file).

    diamond_dict.keys()s gives all unique TRANSCRIPT_NAMES, and
    diamond_dict[TRANSCRIPT_NAME1].keys() gives all unique ORFs for TRANSCRIPT_NAME1.
    """
    diamond_dict = dict()

    with open(diamond_infile_path) as infile:

        for line in infile:

            line = line.rstrip('\n').split("\t")

            # Parse read name, orf, and pident
            read_name = line[0].split("xxx")[0]
            orf = line[0].split("xxx")[1]
            pident = float(line[diamond_fields.index("pident")])

            # Sort based no pident here
            if pident < min_pident:
                continue

            # Initiate and load dictionary and subdictionary
            if read_name not in diamond_dict:
                diamond_dict[read_name] = dict()
            if not orf in diamond_dict[read_name]:
                diamond_dict[read_name][orf] = []

            diamond_dict[read_name][orf].append(line)

    return diamond_dict

def count_most_5prime(diamond_dict, diamond_fields):
    """
    This function iterates through each transcript in the diamond_dict,
    and finds the subject/match of the most 5' ORF of each transcript.
    These matches are then counted, and a Counter object is returned
    with the counts for each subject. This function does NOT
    check for cannonical vs variant. If the 5' ORF has more than one
    match, the match with the longest alignment length will be counted.
    """
    counts = Counter()

    for transcript, transcript_subdict in diamond_dict.items():

        # Determine which ORF is most 5'
        orfs = transcript_subdict.keys()
        fiveprime_orf =  get_most_fiveprime(orfs)

        # If multiple matches for that 5'-most ORF, need to
        # take the one with the longest alignment. If there is
        # only one match, this function will return it.
        matches = transcript_subdict[fiveprime_orf]
        longest_match = get_longest_match(matches, diamond_fields)

        # Add to the count for the longest_match's subject
        subject = longest_match[diamond_fields.index("sseqid")]
        counts[subject] += 1

    return counts

def split_cannonical_and_variant(diamond_dict, diamond_fields, pr_length_dict, min_ok_nterm_extension):
    """
    This function doesn't look at position of the transcript,
    but iterates through each ORF's matches to determine if the
    ORF is cannonical or variant. Variant ORFs are considered
    ORFs where qstart != sstart, or qend != send. IF an ORF has
    multiple matches, it will output the matches as a single tuple
    that will contain each match list, and the tubple will be saved
    to the subject with the longest length.

    The output to this function is therefore a dictoinary of structure
    subject:[[match1], [match2]...]

    As mentioned, if there is a case where there are multiple matches
    of the same ORF to a subject, there will be an instance of a
    tuple containing all matches. This tuple should still be counted
    as a single instance of that variant ORF. This would look like:
    subject:[[match1], [match2], (match3-1, match3-2),...]
    """
    pass


    out_dict = dict()

    for transcript, transcript_subdict in diamond_dict.items():

        # Iterate through every ORF.
        for orf, match_list in transcript_subdict.items():

            # If only one match, check if variant or cannonical
            if len(match_list) == 1:

                # Pull out the subject
                subject = match_list[0][diamond_fields.index("sseqid")]

                # Investigate the entire match to see if it is cannonical
                if is_cannonical(match_list[0], diamond_fields, pr_length_dict, min_ok_nterm_extension):
                    key = subject
                else:
                    key = subject + "_variant"

                if not key in out_dict:
                    out_dict[key] = []

                out_dict[key].append(match_list[0])

            # if more than one match, need to determine which match is longer
            # and assign accordingly
            else:
                longest_match = get_longest_match(match_list, diamond_fields)

                subject = longest_match[diamond_fields.index("sseqid")]

                if is_cannonical(longest_match, diamond_fields, pr_length_dict, min_ok_nterm_extension):
                    key = subject
                else:
                    key = subject + "_variant"

                if not key in out_dict:
                    out_dict[key] = []

                # Here, append all the matches as a single tuple
                out_dict[key].append(tuple(match_list))

    return out_dict

def cannonical_and_variant_dict_to_coord_dict(in_dict, diamond_fields):
    """
    This function takes in the cannonical_and_variant_split dict
    and outputs information as pandas dataframe.
    This input dictionary has keys which are subjects
    (for example "M_variant", "M", "N", and so on). The value is
    a list of matches, where each match is itself a list of the
    contents of a diamond line.

    Notably, if an ORF had two matches, one to the subject but one
    to another entry in the database, both matches would be here and
    that ORF would be present as a tuple of the matches.

    in_dict structure
    -----------------
    subject: [match1, (match2-1, match2-2), match3,...]
    * match1 and match3 are typical, but the ORF in match2 had two
    matches (one of which to this subject), and both matches are kept
    in that tuple.
    * each match is a list of diamond match information with the fields
    detailed in diamond_fields.


    out_df structure
    ----------------
    subject, transcript_name, start, end, strand
    """
    out_info = ''

    for subject, match_list in in_dict.items():

        # Iterate through each match for each subject
        for match in match_list:

            # Handle the tuple cases where there were multiple matches.
            # The match to this subject was set to be the longest one,
            # so need to extract that one here
            if type(match) == tuple:
                match = get_longest_match(list(match), diamond_fields)

            # Parse the required information
            transcript = match[diamond_fields.index("qseqid")].split("xxx")

            name = transcript[0]
            start = transcript[1].split("_")[0]
            end = transcript[1].split("_")[1]
            strand = transcript[1].split("_")[2]

            out_info += "{}\t{}\t{}\t{}\t{}\n".format(subject, name, start, end, strand)

    # Convert to a dataframe
    out_df = pd.read_csv(StringIO(out_info),
                         sep ="\t",
                        )
    out_df.columns = ["subject", "transcript_name", "start", "end", "strand"]

    return out_df

def get_variant_info(cannonical_and_variant_split, diamond_fields):
    """
    cannonical_and_variant_split: key:List of lists. But sometimes,
        the inner item will be a tuple of lists instead of a list
        if there was multiple matches.

    This function will skip the predicted ORFs for now.

    This funciton outputs a dataframe with one row for each match.
    Columns are:
    query, subject, sstart, ssend, orf_type
    """
    output = ""

    for variant, l in cannonical_and_variant_split.items():

        # Ignore predicted ORFs for now..
        if variant.startswith("NC"):
            continue

        # Label variant vs cannonical
        if len(variant.split("_")) == 2:
            orf_type = "variant"
        else:
            orf_type = "cannonical"

        # Iterate through each match
        for sublist in l:

            # Placeholders to deal with fusions
            fusion_subject = ""
            fusion_terminal = ""
            current_orf_type = orf_type

            # Deal with the multi-matched folks. Will extract the ones
            # that matched to the desired ORF.
            if type(sublist) == tuple:

                # IF we're here, this one is a fusion...
                current_orf_type = "fusion"

                desired_subject = variant.split("_variant")[0]

                # Get fusion subject
                fusion_subject = find_fusions_from_match_tuple(sublist, desired_subject, diamond_fields)

                # Get fusion terminal
                fusion_terminal = get_fusion_terminal(sublist, desired_subject, diamond_fields)

                if fusion_terminal == "":
                    print(sublist)
                    print(desired_subject)
                    print()

                # Set primary match
                sublist = get_primary_match(sublist, desired_subject, diamond_fields)


            # Extrand information...
            query = sublist[diamond_fields.index("qseqid")]
            subject = sublist[diamond_fields.index("sseqid")]
            start = sublist[diamond_fields.index("sstart")]
            end = sublist[diamond_fields.index("send")]

            # Add to output
            output += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query, subject, start, end, current_orf_type, fusion_subject, fusion_terminal)

    # Convert to a dataframe
    out_df = pd.read_csv(StringIO(output),
                         sep ="\t",
                        )
    out_df.columns = ["query", "subject", "sstart", "ssend", "orf_type", "fusion_subject", "fusion_terminal"]

    return out_df

#-------------------------------------------#
# Fusion subfunctions for get_variant_info
#-------------------------------------------#

def get_primary_match(match_tuple, subject, diamond_fields):
    """
    If given a tuple of multiple matches, will return the primary match
    that has the desired subject
    """
    for match in match_tuple:
        if match[diamond_fields.index("sseqid")] == subject:
            return match

def find_fusions_from_match_tuple(match_tuple, subject, diamond_fields):
    """
    Finds subject of fusions. Will return 'multiple' if there
    are multiple fusions. Otherwise, will return the fusion.
    """
    fusion_subject = []

    for match in match_tuple:
        if match[diamond_fields.index("sseqid")] == subject:
            continue

        else:
            fusion_subject.append(match[diamond_fields.index("sseqid")])

    if len(fusion_subject) > 1:
        return "multiple"
    elif fusion_subject == []:
        # If it got here, it's fused to itself...
        return "self"
    else:
        return fusion_subject[0]


def get_fusion_terminal(match_tuple, subject, diamond_fields):
    """
    First finds the query coordinates of the main match, and then finds
    query coordinats of the fusions. Reports if the fusion is N, C, or both
    """

    # First find primary match
    primary_match = ""
    for match in match_tuple:
        if match[diamond_fields.index("sseqid")] == subject:
            primary_match = match
            break

    # Find the primary match start and end coords
    primary_match_start = primary_match[diamond_fields.index("qstart")]
    primary_match_end = primary_match[diamond_fields.index("qend")]

    # Now determine N/C terminal status, assuming they are not overlapping
    terminal = ""
    for match in match_tuple:

        # Skip the primary match
        if match[diamond_fields.index("sseqid")] == subject:
            continue

        # Determine current qstart/qend
        start = match[diamond_fields.index("qstart")]
        end = match[diamond_fields.index("qend")]

        term = ""
        # If start is after primary_match_end, it is C-terminal
        if start > primary_match_end:
            term = "C"

        # If the end is before the primary_match_start, it is N-terminal
        if end < primary_match_start:
            term = "N"

        # If terminal is currently empty, add term.
        if terminal == "":
            terminal = term

        # Otherwise, check if we need to change it to both
        elif terminal == term:
            continue

        else:
            terminal = "both"

        # If terminal is still blank, they're overlapping. This may have to do with
        # noise on the edge of the diamond alignments. Therefore, will compare their
        # starts and ends.
        if terminal == "":

            # Check to see which is the most on the N/C side...
            if primary_match_start > start:
                terminal = "N"
            elif primary_match_end < end:
                terminal = "C"
            else:
                terminal = "Ambiguous"

    # If terminal is still "", then it aborted when checking for self.
    # Thus, this is a self fusion, and N/C is ambiguous/irrelevant.
    if terminal == "":
        terminal = "self"

    return terminal

#==============================================================================#
# Constants
#==============================================================================#
pr_length_dict = {
    "ORF1a": 4405,
    "ORF1b": 2595,
    "S": 1273,
    "ORF3a": 275,
    "E": 75,
    "M": 222,
    "ORF6": 61,
    "ORF7a": 121,
    "ORF7b": 43,
    "ORF8": 121,
    "N": 419,
    "ORF10": 38
}

# Diamond fields that are required.
required_fields = {'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send'}

#==============================================================================#
# Main
#==============================================================================#
def main():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to parse the DIAMOND output file
            generated by mapping transcript-predicted ORFs against SARS-CoV-2.
            This script makes several important output files:

            1) PREFIX_subject_starts_and_ends.tsv:
                This outputs the predicted start and stop positions of each
                query ORF based on their alignment starts and stops against
                cannonical proteins. So, this is  *protein level*
                information of all variant and canonical ORFs. This output file
                also gives information on any ORF fusions, if applicable.
            2) PREFIX_total_ORF_counts.tsv:
                This file contains total counts for a given ORF. This does NOT
                consider whether the ORF is canonical or variant. This DOES
                include all unannotated ORFs.
            3) PREFIX_fiveprime-most-orf_counts.tsv:
                This file lists the number of transcripts which contain each
                cannonical and variant ORF as the most 5' ORF. Does include
                unannotated ORFs.
            4) PREFIX_cannonical_vs_variant_counts.tsv:
                This file contains total counts for variant vs canonical ORFs.
                This does NOT include unannotated ORFs.
            5) PREFIX_cannonical_vs_variant_startsites.tsv:
                This file contains transcript ORF start and end sites for
                variant vs canonical ORFs. So, this is the nucleotide start and
                end points for each ORF on each transcript. This does not
                include unannotated ORFs.
            """)

    # Required arguments
    parser.add_argument(
        '-d',
        '--diamond_infile_path',
        type=str,
        required=True,
        help='''
        Path to the input diamond file containing the alignments.
        '''
    )
    parser.add_argument(
        '-o',
        '--output_prefix',
        type=str,
        required=True,
        help='''
        Prefix of the output files. Any required directories will be created.
        Outputs will be of format {PREFIX}_output_file_name.tsv
        '''
    )

    # Optional arguments
    parser.add_argument(
        '-f',
        '--diamond_fields',
        type=str,
        required=False,
        default = 'qseqid qlen sseqid evalue bitscore pident length qstart qend sstart send',
        help='''
        Fields of the diamond file. Should be a space-delimited list.
        '''
    )
    parser.add_argument(
        '-p',
        '--min_pident',
        type=int,
        required=False,
        default = 95,
        help='''
        The minimum percent identity for an alignment to count. Percent should
        be a whole integer, not a decimal. Alignments below this threshold
        are not even parsed from the diamond infile.
        '''
    )
    parser.add_argument(
        '-n',
        '--min_ok_nterm_extension',
        type=int,
        required=False,
        default = 20,
        help='''
        Sometimes an alignment is *almost* perfectly canonical, but kind of
        looks varriant due to an elongated N-terminal side. This is caused
        because ORFs are predicted based on any nTG, while the canonical
        proteins are definted at a specific ATG. Thus, this switch addresses
        cases an ORF is essentially canonical but just was predicted from an
        earlier start codon. The integer value entered is how many amino acids
        extra on the N terminus are permissible. In any case, the stop position
        of the ORF in question must be the expected stop.
        '''
    )

    args = parser.parse_args()

    # Define input variables
    diamond_infile_path = args.diamond_infile_path
    output_prefix = args.output_prefix
    diamond_fields = args.diamond_fields.split(" ")
    min_pident = args.min_pident
    min_ok_nterm_extension = args.min_ok_nterm_extension

    # Make sure all required fields are in diamond_fields
    validate_fields(diamond_fields, required_fields)

    # Read in the DIAMOND file
    diamond_dict = parse_diamond(diamond_infile_path, diamond_fields, min_pident)

    # Get the counts of the most 5' ORF in each transcript
    fiveprime_counts = count_most_5prime(diamond_dict, diamond_fields)

    # Split the cannonical and variants
    cannonical_and_variant_split = split_cannonical_and_variant(diamond_dict, diamond_fields, pr_length_dict, min_ok_nterm_extension)
    cannonical_and_variant_split_counts = dict_value_counts(cannonical_and_variant_split)
    cannonical_and_variant_split_coord_df = cannonical_and_variant_dict_to_coord_dict(cannonical_and_variant_split, diamond_fields)


    # Extract information of the subject start and end sites for each match for use
    # making a histogram in R
    varaint_starts_and_ends = get_variant_info(cannonical_and_variant_split, diamond_fields)

    # Make output directory if necessary
    out_dir = os.path.dirname(output_prefix)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Write output files

    # 1)
    varaint_starts_and_ends.to_csv('{}_subject_starts_and_ends.tsv'.format(output_prefix), index=False, sep="\t")

    # 2)
    # Aggregate cannonical and variant counts in cannonical_and_variant_split_counts, and output them
    counts = aggregate_cannonical_and_varaint_split_counts(cannonical_and_variant_split_counts)
    counts = pd.DataFrame.from_dict(counts, orient="index")
    counts = counts.reset_index()
    counts.columns = ["ORF", 'Counts']
    counts.to_csv('{}_total_ORF_counts.tsv'.format(output_prefix), index=False, sep="\t")

    # 3)
    fiveprime_counts_df = pd.DataFrame.from_dict(fiveprime_counts, orient="index")
    fiveprime_counts_df = fiveprime_counts_df.reset_index()
    fiveprime_counts_df.columns = ["ORF", 'Counts']
    fiveprime_counts_df.to_csv('{}_fiveprime-most-orf_counts.tsv'.format(output_prefix), index=False, sep="\t")


    # 4)
    cannonical_and_variant_df_COUNT = pd.DataFrame.from_dict(cannonical_and_variant_split_counts, orient="index")
    cannonical_and_variant_df_COUNT = cannonical_and_variant_df_COUNT.reset_index()
    cannonical_and_variant_df_COUNT.columns = ["ORF", 'Counts']
    cannonical_and_variant_df_COUNT_subset = cannonical_and_variant_df_COUNT[~cannonical_and_variant_df_COUNT["ORF"].str.contains("NC_045512.2")]
    cannonical_and_variant_df_COUNT_subset.to_csv('{}_cannonical_vs_variant_counts.tsv'.format(output_prefix), index=False, sep="\t")


    # 5)
    cannonical_and_variant_split_coord_df_subset = cannonical_and_variant_split_coord_df[~cannonical_and_variant_split_coord_df["subject"].str.contains("NC_045512.2")]
    cannonical_and_variant_split_coord_df_subset.to_csv('{}_cannonical_vs_variant_startsites.tsv'.format(output_prefix), index=False, sep="\t")

if __name__ == '__main__':
    main()
