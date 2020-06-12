#!/bin/bash

#SBATCH -t 0-2:00:0
#SBATCH -p short
#SBATCH --mem=30GB
#SBATCH -c 5

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose is to map input short reads to a reference index using
        the STAR aligner. Settings are adapted from Kim et al's settings and
        is desired to limit penalties for non-canonical read junctions.

        USAGE: $0

        Required:
        -i <INPUT_FASTQ>
            Input file containing reads in fastq format.
        -d <INDEX>
            Path to the star index directory.
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
while getopts i:d:o:t: option ; do
        case "${option}"
        in
                i) INPUT_FASTQ=${OPTARG};;
                d) INDEX=${OPTARG};;
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

INPUT_FASTQ: $INPUT_FASTQ
INDEX: $INDEX
OUT_BAM: $OUT_BAM
THREADS: $THREADS
"


echo "Starting script."

# Make directories if necessary
mkdir -p $(dirname $OUT_BAM)

SAMPLE=$(basename ${OUT_BAM%.bam})

STAR \
--runThreadN 5 \
--genomeDir $INDEX \
--readFilesIn $INPUT_FASTQ \
--outFileNamePrefix $(dirname $OUT_BAM)/${SAMPLE}_working/${SAMPLE}_ \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--outSJfilterOverhangMin 12 12 12 12 \
--outSJfilterCountUniqueMin 1 1 1 1 \
--outSJfilterCountTotalMin 1 1 1 1 \
--outSJfilterDistToOtherSJmin 0 0 0 0 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04\
 --scoreGapNoncan -4 \
--scoreGapATAC -4 \
--chimOutType WithinBAM HardClip \
--chimScoreJunctionNonGTAG 0 \
--alignSJstitchMismatchNmax -1 -1 -1 -1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000

# Extract just the bam
mv $(dirname $OUT_BAM)/${SAMPLE}_working/${SAMPLE}*bam $OUT_BAM

echo Finished
