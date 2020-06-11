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
        -m <MODE> [long]
            Options: 'long', 'short'
            - long:
              Maps input long (nanopore) reads using settings optimized
              for a small viral reference like SARS-Cov-2
            - short:
              Maps input *short reads* against a desired target. This simply
              engages the sr (genomic short-read mapping) preset. No other
              mapping settings are changed.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts i:r:o:t:m: option ; do
        case "${option}"
        in
                i) INPUT_FASTX=${OPTARG};;
                r) REF_FASTA=${OPTARG};;
                o) OUT_BAM=${OPTARG};;
                t) THREADS=${OPTARG};;
                m) MODE=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Set defaults
#------------------------------------------------------------------------------#
THREADS=${THREADS:-1}
MODE=${MODE:-long}

#------------------------------------------------------------------------------#
# Run cd-hit
#------------------------------------------------------------------------------#
echo "
INPUTS -

INPUT_FASTX: $INPUT_FASTX
REF_FASTA: $REF_FASTA
OUT_BAM: $OUT_BAM
THREADS: $THREADS
MODE: $MODE
"


echo "Starting script."

# Make directories if necessary
mkdir -p $(dirname $OUT_BAM)

# Run minimap
if [[ $MODE == "long" ]] ; then
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

elif [[ $MODE == "short" ]] ; then
  minimap2 \
  -t $THREADS \
  -a \
  -x sr \
  --sam-hit-only \
  $REF_FASTA \
  $INPUT_FASTX | samtools view -bh -@ $THREADS | samtools sort -@ $THREADS > $OUT_BAM
  samtools index $OUT_BAM
else
  echo "MODE must be set to 'long' or 'short'. You specified $MODE"
  exit 1
fi

echo "Finished."
