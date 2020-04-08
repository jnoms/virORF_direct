#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import pathlib
import os
import argparse
import gzip
import sys

#==============================================================================#
# Auxiliary functions
#==============================================================================#
def there_is_overlap_quick(pr1, pr2):
    """
    This just checks if pr1 overlaps pr2 or vice-versa. If there
    is overlap the function returns True, else it returns False.

    This is a quick search! But, it can miss some cases.
    """
    if pr1 in pr2:
        return True
    elif pr2 in pr1:
        return True
    else:
        return False

def ORF_to_sequence(ORF, ref_seq, translate=False):
    """
    ORF is a tuple of structure (start, end, strand).
    """
    # Extract, start, end, and strand
    (start, end, strand) = ORF
    start = int(start)
    end = int(end)

    # Find ORF sequence based on reference sequence
    ORF_dna_seq = ref_seq[start-1:end]

    # If ORF is on negative strand, need to translate
    if strand == "-":
        ORF_dna_seq = str(Seq(ORF_dna_seq).reverse_complement())

    # If not translating, return the nucleotide seq
    if translate == False:
        return ORF_dna_seq

    # Otherwise, return the protein sequence
    else:
        return str(Seq(ORF_dna_seq).translate())[:-1]

#==============================================================================#
# Main functions
#==============================================================================#
def read_fasta_to_memory(input_fasta):
    """
    Reads fasta into a memory as a dictionary with header:sequence.
    This function can handle .gzip files, but input_fasta needs to
    end with .gz
    """
    fasta_dict = dict()

    if not input_fasta.endswith(".gz"):
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            fasta_dict[seq_record.id] = str(seq_record.seq)

    elif input_fasta.endswith(".gz"):
        with gzip.open(input_fasta, "rt") as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                fasta_dict[seq_record.id] = str(seq_record.seq)

    return fasta_dict


def parse_prodigal(prodigal_path, pos_strand_only=True):
    """
    Prodigal_path is the path to the prodigal file. If it is
    gzipped, file path must end in .gz.

    This function returns a dictionary where every key is the
    accession of a sequence, and the value is a LIST of tuples,
    where each tuple is (begining, end, strand) of each ORF.dd
    """
    # Read in contents to memory
    prodigal_contents = ""
    if prodigal_path.endswith(".gz"):
        with gzip.open(prodigal_path, "rt") as infile:
            prodigal_contents = infile.read()
    else:
        with open(prodigal_path) as infile:
            prodigal_contents = infile.read()

    # Make output dictionary, which will hold accessions and
    # a list of ORFs
    accesion_to_orfs = dict()

    # Iterate through each sequence in the prodigal file
    for entry in prodigal_contents.split("# Sequence Data: "):

        # Currently - each entry is a prodigal output for a given input sequence.
        # There are lines (separated by \n) and each line has tab-delimited
        # columns. The exception is the first line, which starts with "seqnum".
        # This is semicolon delimited, and the accession can be parsed from it.

        # Remove first empty line, which isn't actually an entry
        if entry == "":
            continue

        # Make a placeholder to hold the accession
        current_accession = ""

        # Split entry into lines, and start parsing
        for line in entry.split('\n'):

            # Ignore any blank lines
            if line == '':
                continue

            # Parse the first line - get the accession
            if line.startswith("seqnum"):

                # Isolate the header
                header = line.split("seqhdr=")[1].rstrip('\n').strip('"')

                # Get the accession
                accession = header.split(" ")[0]

                # Save it as "current accession"
                current_accession = accession

                # Initiate the accession_to_orfs dictionary for that accession
                accesion_to_orfs[current_accession] = []

                continue

            # Ignore the "# Run Data:" line
            if line.startswith("# Run Data:"):
                continue

            # Ignore the header line
            if line.startswith("Beg"):
                continue

            # Find the begining, end, and strand
            begining = int(line.split("\t")[0])
            end = int(line.split("\t")[1])
            strand = line.split("\t")[2]

            if pos_strand_only == True:
                if strand == "-":
                    continue

            # Write to orf dictionary
            accesion_to_orfs[current_accession].append((begining, end, strand))

    return accesion_to_orfs

def remove_shorter_orf_with_same_start_or_stop(orf_list):
    """
    This function is used as a precursor in remove_overlapping_orfs.
    This function simply finds if two orfs in an orf_list have the
    same start of stop, and then only keeps the longer one.
    """

    reduced_orf_list = []

    for orf in orf_list:
        start,stop = orf[0], orf[1]
        length = stop - start

        # init the output list
        if reduced_orf_list == []:
            reduced_orf_list.append(orf)
            continue

        # find previous start/stop and length
        previous_orf = reduced_orf_list[-1]
        previous_start, previous_stop = previous_orf[0], previous_orf[1]
        previous_length = previous_stop - previous_start

        # if this orf has the same start of stop as the previous orf,
        # check if the length is longer. If so, overwrite. If not,
        # continue.
        if start == previous_start or stop == previous_stop:

            if length > previous_length:
                reduced_orf_list[-1] = orf
                continue

        else:
            reduced_orf_list.append(orf)

    return reduced_orf_list

def remove_overlapping_orfs(ORF_list, ref_seq):
    """
    ORF_list is a list of ORFs, while the ref_seq is
    the sequence the ORFs where derived from. This function
    translates each ORF into protein sequence, and then
    determins if any of the ORFs are completely contained
    within one another (a fast process). If so, it takes
    the longer one.
    """

    # Make dict for storage
    current_orfs_and_seqs = dict()


    # First do a quick function where, if multiple ORFs have the
    # same start of stop, takes the longer one.
    ORF_list = remove_shorter_orf_with_same_start_or_stop(ORF_list)

    # Process each ORF in the list
    for ORF in ORF_list:

        # translate
        pr = ORF_to_sequence(ORF, ref_seq, translate=True)

        #Handle first one
        if current_orfs_and_seqs == dict():
            current_orfs_and_seqs[ORF] = pr

        # Next, determine if the ORF has a protein seq overlapping one of the
        # current protein sequences. Will iterate through the contents of
        # current_orfs_and_seqs and, for each entry there, check if overlap.

        # make a status indicator incase there isn't overlap
        status = ''

        for collected_orf, collected_pr in current_orfs_and_seqs.items():

            # if this sequence didn't overlap, continue
            if not there_is_overlap_quick(pr, collected_pr):
                continue

            # If there is overlap...
            if there_is_overlap_quick(pr, collected_pr):

                # If current protein seq is longer, replace the old one with this one
                if len(pr) > len(collected_pr):

                    # Remove smaller one
                    del current_orfs_and_seqs[collected_orf]

                    # Add this one
                    current_orfs_and_seqs[ORF] = pr

                    # Update status
                    status = "Done"

                    # No need to continue with this loop.
                    break

                # If the stored one is bigger, don't overwrite it and break the llop
                if len(pr) <= len(collected_pr):
                    status="Done"
                    break

        # if the status is empty, there was no overlap anywhere so we're good to add
        if status == "":
            current_orfs_and_seqs[ORF] = pr

    # Can go ahead and return the ORFs - can generate the protein sequences again when
    # needed
    return list(current_orfs_and_seqs.keys())

def remove_overlapping_orfs_from_dict(header_to_orfs, header_to_sequence):
    """
    This function calls remove_overlapping_orfs a bunch of time to format
    each orf list in header_to_orfs. header_to_sequence is used to get
    the reference sequence for each header.
    """

    # Keeping track of progress
    total = len(header_to_orfs)
    progress = 0

    non_overlapping_orf_dict = dict()
    for header, orf_list in header_to_orfs.items():

        sequence = header_to_sequence[header]
        non_overlapping_orf_dict[header] = remove_overlapping_orfs(orf_list, sequence)

        # Mark progress at every 10,000 headers written out.
        progress += 1
        if progress%10000/100000 == 0:
            print("Progress: {}/{}".format(progress, total))

    return non_overlapping_orf_dict

def write_pr_and_nt_fastas(header_to_orfs,
                           header_to_sequence,
                           nt_fasta_path,
                           pr_fasta_path,
                           force_methionine_starts = False
                          ):

    # Determine if we're writing out nt and/or pr
    types = set()
    if nt_fasta_path != '':
        types.add("nt")
    if pr_fasta_path != '':
        types.add("pr")

    # Generate the output directories if necessary, and open files.
    if "nt" in types:
        out_dir = os.path.dirname(nt_fasta_path)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
        nt_outfile_handle = open(nt_fasta_path, "w")


    if "pr" in types:
        out_dir = os.path.dirname(pr_fasta_path)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
        pr_outfile_handle = open(pr_fasta_path, "w")

    # Make items for keeping progress
    total = len(header_to_orfs)
    progress = 0

    for header, ORF_list in header_to_orfs.items():

        # Get the ref_seq
        ref_seq = header_to_sequence[header]

        # Parse the ORF_list
        for ORF in ORF_list:

            (start, end, strand) = ORF

            # Define the header
            out_header = ">{}xxx{}_{}_{}".format(header,
                                             start,
                                             end,
                                             strand,
                                             )

            # Find nucleotide and protein sequences
            nt = ORF_to_sequence(ORF, ref_seq, translate=False)
            pr = ORF_to_sequence(ORF, ref_seq, translate=True)

            # If methioine starts are required, change the pr sequences
            if force_methionine_starts == True:
                pr_list = list(pr)
                pr_list[0] = "M"
                pr = "".join(pr_list)

            # Write to outfiles
            if "nt" in types:
                nt_outfile_handle.write(out_header + "\n" + nt + "\n")
            if "pr" in types:
                pr_outfile_handle.write(out_header + "\n" + pr + "\n")

        # Mark progress at every 10,000 headers written out.
        progress += 1
        if progress%10000/100000 == 0:
            print("Progress: {}/{}".format(progress, total))

    if "nt" in types:
        nt_outfile_handle.close()
    if "pr" in types:
        pr_outfile_handle.close()

def main():
    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is generate protein and/or nucleotide
            sequences from prodigal-computed orfs directly from the input fasta
            file. This function also checks if orfs predicted by prodigal are
            overlapping, and if so, only outputs the longest one.
            """)

    # Required arguments
    parser.add_argument(
        '-f',
        '--fasta',
        type=str,
        required=True,
        help='''
        Path to the input fasta containing the sequences input to prodigal.
        '''
    )

    parser.add_argument(
        '-p',
        '--prodigal_file',
        type=str,
        required=True,
        help='''
        Path to the prodigal file that details all possible ORFs. Should have
        been generated with prodigal using the -s switch. Can contain the
        results of multiple queries - i.e. you can plug a multifasta into
        prodigal. Node, that the accession needs to be the first thing in the
        input fasta to prodigal so the accession is reported in the prodigal
        output file.
        '''
    )

    # Optional
    parser.add_argument(
        '-N',
        '--nt_fasta_out_path',
        type=str,
        required=False,
        default="",
        help='''
        Path to the output file containing nucleotide sequences of the ORFs.
        nt_fasta_out_path and/or pr_fasta_out_path must be input.
        '''
    )

    parser.add_argument(
        '-P',
        '--pr_fasta_out_path',
        type=str,
        required=False,
        default="",
        help='''
        nt_fasta_out_path and/or pr_fasta_out_path must be input.
        '''
        )

    parser.add_argument(
        '-m',
        '--force_methionine_starts',
        type=bool,
        required=False,
        default=False,
        help='''
        If option is set to True, will make all protein sequences start with
        methionine even if there is a non-cannonical start codon. Default False.
        '''
        )

    parser.add_argument(
        '-s',
        '--pos_strand_only',
        type=bool,
        required=False,
        default=True,
        help='''
        If option is set to True, will only parse ORFs from prodigal file that
        are on the + strand.
        '''
        )

    args = parser.parse_args()

    # Define input variables
    fasta = args.fasta
    prodigal_file = args.prodigal_file
    nt_fasta_out_path = args.nt_fasta_out_path
    pr_fasta_out_path = args.pr_fasta_out_path
    force_methionine_starts = args.force_methionine_starts
    pos_strand_only = args.pos_strand_only

    # Parse fasta file into memory
    print("{}: Starting.".format(sys.argv[0]))
    print("{}: Reading in fasta to memory.".format(sys.argv[0]))
    header_to_sequence = read_fasta_to_memory(fasta)

    # Parse prodigal file
    print("{}: Reading in prodigal file to memory.".format(sys.argv[0]))
    header_to_orfs = parse_prodigal(prodigal_file, pos_strand_only)

    # Remove overlapping ORFs from each sequence in header_to_orfs
    print("{}: Removing overlapping ORFs.".format(sys.argv[0]))
    header_to_orfs = remove_overlapping_orfs_from_dict(header_to_orfs, header_to_sequence)

    # Write out the nucleotide and protein sequences for each ORF.
    print("{}: Writing output files.".format(sys.argv[0]))
    write_pr_and_nt_fastas(header_to_orfs,
                           header_to_sequence,
                           nt_fasta_out_path,
                           pr_fasta_out_path,
                           force_methionine_starts == force_methionine_starts
                          )
    print("{}: Finished!".format(sys.argv[0]))

if __name__ == '__main__':
    main()
