#!/bin/bash

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose is to use fmlrc to correct long nanopore reads with
        illumina reads.

        USAGE: $0

        Required:
        -l <LONG_READS>
            Input file containing dRNAseq reads. Can be in either fastq or
            fasta format, doesn't matter.
        -s <SHORT_READS>
            Path to the fasta containing the reference genome. Short reads
            should be gzipped!
        -o <OUTPUT_PATH>
            Path to the resultant output bam. Will contain only the mapped reads.

        Optional:
        -t <THREADS> [1]
            Number of computing threads available.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts l:s:o:t: option ; do
        case "${option}"
        in
                l) LONG_READS=${OPTARG};;
                s) SHORT_READS=${OPTARG};;
                o) OUTPUT_PATH=${OPTARG};;
                t) THREADS=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Set defaults
#------------------------------------------------------------------------------#
THREADS=${THREADS:-1}

#------------------------------------------------------------------------------#
# Run cd-hit
#------------------------------------------------------------------------------#
echo "
INPUTS -

LONG_READS: $LONG_READS
SHORT_READS: $SHORT_READS
OUTPUT_PATH: $OUTPUT_PATH
THREADS: $THREADS
"

echo "Starting script."

# Make directories if necessary
mkdir -p $(dirname $OUTPUT_PATH)

# Generate fmlrc index from the short reads
gunzip -c $SHORT_READS |\
  awk 'NR % 4 == 2' |\
  sort | tr NT TN |\
  ropebwt2 -LR | tr NT TN |\
  fmlrc-convert $(basename $SHORT_READS)_msbwt.npy

# Correct reads using fmlrc
fmlrc \
-p 5 \
$(basename $SHORT_READS)_msbwt.npy \
$LONG_READS $OUTPUT_PATH

echo "Finished."
