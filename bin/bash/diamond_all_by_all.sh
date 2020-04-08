#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2

module load gcc

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#

usage() {
        echo "
        This purpose of this file is to use DIAMOND blastp to map all sequences
        in the input fasta against a database made of those sequences. Input
        sequences are expected to be protein.

        USAGE: $0

        Required:
        -i <IN_FASTA>
            Path to the multifasta containing protein sequences.
        -o <OUTPUT_PATH>
            Path to the output DIAMOND file.

        Optional:
        -m <MEMORY> [10]
            Memory available in GB. DIAMOND blocksize is set to MEMORY/10 or, if
            MEMORY <10, to 1.
        -t <TEMP_DIR> [./diamond_temp]
            Path to the temporary directory
        -e <EVALUE> [0.001]
            Max evalue.
        -f <OUT_FORMAT>  ['6 qseqid sseqid evalue bitscore pident length']
            Output format.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts i:o:m:t:e:f: option ; do
        case "${option}"
        in
                i) IN_FASTA=${OPTARG};;
                o) OUTPUT_PATH=${OPTARG};;
                m) MEMORY=${OPTARG};;
                t) TEMP_DIR=${OPTARG};;
                e) EVALUE=${OPTARG};;
                f) OUT_FORMAT=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Setting defaults
#------------------------------------------------------------------------------#
MEMORY=${MEMORY:-10}
TEMP_DIR=${TEMP_DIR:-./diamond_temp}
EVALUE=${EVALUE:-0.001}
OUT_FORMAT=${OUT_FORMAT:-"6 qseqid sseqid evalue bitscore pident length"}

# Calculate DIAMOND block size
if [[ $MEMORY < 10 ]] ; then
  BLOCK=1
else
  BLOCK=$((MEMORY/10))
fi

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#

# Make necessary directories
#------------------------------------------------------------------------------#
mkdir -p $(dirname $OUTPUT_PATH)
mkdir -p ref
mkdir -p $TEMP_DIR

# Generate DIAMOND database
#------------------------------------------------------------------------------#
echo "Making DIAMOND database."
date

diamond makedb \
--in $IN_FASTA \
-d ref/$(basename ${IN_FASTA}).dmnd

# Running DIAMOND
#------------------------------------------------------------------------------#
echo "Running DIAMOND."
date

#Running DIAMOND
diamond blastp \
-d ref/$(basename ${IN_FASTA}).dmnd \
-q $IN_FASTA \
--sensitive \
-o ${OUTPUT_PATH}.tmp \
--tmpdir $TEMP_DIR \
--evalue $EVALUE \
--outfmt $OUT_FORMAT \
--index-chunks 1 \
--block-size $BLOCK \
--max-target-seqs 0

# Quit if DIAMOND failed
if [ ! -f ${OUTPUT_PATH}.tmp ] ; then
  echo "Cannot find the output file $OUTPUT_PATH."
  echo "Even if there are no alignments, there should be an empty output file..."
  echo "This suggests someone went wrong."
  exit 1
fi

# Keep only lines where the query didn't match to itself
awk '$1 != $2 { print $0 }' ${OUTPUT_PATH}.tmp > ${OUTPUT_PATH} && \
  rm ${OUTPUT_PATH}.tmp

echo "Script complete."
date
