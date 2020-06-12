#!/bin/bash

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose is to use minimap to map input nanopore direct RNAseq
        (dRNAseq) reads to the SARS-CoV-2 reference genome fasta. This script
        sets specific parameters that are optimized for SARS-CoV-2/other ~30kB
        small viruses, so consider that prior to using this script to map
        against other things.

        USAGE: $0

        Required:
        -i <INPUT_FASTX>
            Input file containing dRNAseq reads. Can be in either fastq or
            fasta format, doesn't matter.
        -r <REF_FASTA>
            Path to the fasta containing the reference genome.
        -o <OUT_BAM>
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
while getopts i:r:o:t: option ; do
        case "${option}"
        in
                i) INPUT_FASTX=${OPTARG};;
                r) REF_FASTA=${OPTARG};;
                o) OUT_BAM=${OPTARG};;
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

INPUT_FASTX: $INPUT_FASTX
REF_FASTA: $REF_FASTA
OUT_BAM: $OUT_BAM
THREADS: $THREADS
"


echo "Starting script."

# Make directories if necessary
mkdir -p $(dirname $OUT_BAM)

# Run minimap
minimap2 \
-a \
-x splice \
-t $THREADS \
--sam-hit-only \
-k 8 \
-w 1 \
--splice \
-g 30000 \
-G 30000 \
-A1 \
-B2 \
-O2,24 \
-E1,0 \
--no-end-flt \
-F 40000 \
-N 32 \
--splice-flank=no \
--max-chain-skip=40 \
-p 0.7 \
$REF_FASTA $INPUT_FASTX | samtools view -bh -@ $THREADS | samtools sort -@ $THREADS > $OUT_BAM
samtools index $OUT_BAM

echo "Finished."
