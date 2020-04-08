#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2

set -e

#------------------------------------------------------------------------------#
#Defining usage and setting input
#------------------------------------------------------------------------------#

usage() {
        echo "
        This purpose of this file is to extract reads from a fasta that contain
        the leader sequence. The position and percent identity of the leader
        sequence can be specified. This script detects the leader via blastn.

        USAGE: $0

        Required:
        -q <QUERY>
            Path to the query file in fasta format. This is the file containing
            the reads that you want to detect the leader in.
        -s <SUBJECT_FA>
            Path to the fasta containing the leader sequence, in fasta format.
        -o <OUTPUT_FA>
            Path to the output fasta that contains only reads that have the
            leader.

        Optional:
        -b <BLAST_TYPE> [blastn]
            Options are 'blastn', 'megablast', 'dc-megablast'...
        -t <THREADS> [1]
            Number of computing threads
        -S <SAMPLE_ID> [\$(basename \$QUERY)]
            Name of the sample - used for naming some intermediate files.
        -a <MIN_ALIGN_LENGTH> [45]
            Minimum alignment length between leader and read.
        -b <MAX_START_POS> [50]
            Maximum position at which the leader can start on the read.
        -c <MIN_PIDENT> [0.85]
            Minimum percent identity between the leader and the aligning
            sequence of the read.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts q:s:o:B:t:S:a:b:c: option ; do
        case "${option}"
        in
                q) QUERY=${OPTARG};;
                s) SUBJECT_FA=${OPTARG};;
                o) OUTPUT_FA=${OPTARG};;
                B) BLAST_TYPE=${OPTARG};;
                t) THREADS=${OPTARG};;
                S) SAMPLE_ID=${OPTARG};;
                a) MIN_ALIGN_LENGTH=${OPTARG};;
                b) MAX_START_POS=${OPTARG};;
                c) MIN_PIDENT=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Setting defaults
#------------------------------------------------------------------------------#
BLAST_TYPE=${BLAST_TYPE:-blastn}
THREADS=${THREADS:-1}
SAMPLE_ID=${SAMPLE_ID:-$(basename $QUERY)}
MIN_ALIGN_LENGTH=${MIN_ALIGN_LENGTH:-45}
MAX_START_POS=${MAX_START_POS:-50}
MIN_PIDENT=${MIN_PIDENT:-0.85}

# Constants
OUT_FORMAT='6 qseqid sseqid evalue bitscore pident qstart qend length'
MAX_HSPS=1
MAX_TARGETS=1

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
echo "$0: Starting. Running blast."
blastn \
-task $BLAST_TYPE \
-query $QUERY \
-subject $SUBJECT_FA \
-outfmt "$OUT_FORMAT" \
-max_hsps $MAX_HSPS \
-max_target_seqs $MAX_TARGETS \
-out ${SAMPLE_ID}_blast.txt.tmp \
-num_threads $THREADS

# Now, look for reads that meet alignment length, start position, and pident
# thresholds.
awk \
  -v MIN_ALIGN_LENGTH=$MIN_ALIGN_LENGTH \
  -v MAX_START_POS=$MAX_START_POS \
  -v MIN_PIDENT=$MIN_PIDENT \
  '($8 >= MIN_ALIGN_LENGTH && $6 < MAX_START_POS && $5 >= MIN_PIDENT)' \
  ${SAMPLE_ID}_blast.txt.tmp > ${SAMPLE_ID}_blast.txt &&\
  rm ${SAMPLE_ID}_blast.txt.tmp || echo "AWK DIDNT WORK"

# Extract only the passing reads using seqtk
echo "$0: Extracting passing reads."
cut -f1 ${SAMPLE_ID}_blast.txt > ${SAMPLE_ID}_passing_reads.txt
seqtk subseq $QUERY ${SAMPLE_ID}_passing_reads.txt > $OUTPUT_FA

echo "$0: Finished."
