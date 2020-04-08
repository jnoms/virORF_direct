#!/bin/bash

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose of this script is to take an input fasta containing
        protein sequences and group them using cd-hit. Automatically uses all
        available threads.

        USAGE: $0

        Required:
        -i <INFILE>           Path to the input fasta.
        -o <CLSTR_FILE>       Path to the cluster information outfile.
        -r <REPRESENTATIVES>  Path to fasta outfile with cluster representatives.

        Optional:
        -n <WORD_SIZE> Word size.[2]
        -c <MIN_SEQ_ID> Minimum sequence identity for clustering. [0.5]
        -M <MEMORY> Memory in megabytes. [10000]
        -l <MIN_PROTEIN_LEN> Minimum protein length [5]
        -g <MODE>
           1=Cluster into most similar - more computationally intensive
           0=Cluster into first matching - faster [0]
        -s <SORT_STRATEGY>
           1=Sort cluster file such that groups with more members are at
           the top.
           0=Largest sequences at top. cd-hit param is -sc. [1]
        -a <MIN_OVERLAP_LEN> Minimum percent of the larger sequence the smaller
            sequnce needs to cover for it to be grouped. cd-hit param is -aL [0.7]
        -S <SHORTNESS> [0]
            How short the smaller sequence can be relative to the cluster
            representative. For example, a value of 0.9 means the short sequence
            must be at least 90% the length of the longer sequence. Default is 0
            which means no cutoff.
        "
}

#If less than 2 options are input, show usage and exit script.
if [ $# -le 2 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts i:o:r:n:c:M:l:g:s:a:S: option ; do
        case "${option}"
        in
                i) INFILE=${OPTARG};;
                o) CLSTR_FILE=${OPTARG};;
                r) REPRESENTATIVES=${OPTARG};;

                n) WORD_SIZE=${OPTARG};;
                c) MIN_SEQ_ID=${OPTARG};;
                M) MEMORY=${OPTARG};;
                l) MIN_PROTEIN_LEN=${OPTARG};;
                g) MODE=${OPTARG};;
                s) SORT_STRATEGY=${OPTARG};;
                a) MIN_OVERLAP_LEN=${OPTARG};;
                S) SHORTNESS=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Set default parameters
#------------------------------------------------------------------------------#
WORD_SIZE=${WORD_SIZE:-2}
MIN_SEQ_ID=${MIN_SEQ_ID:-0.5}
MEMORY=${MEMORY:-10000}
MIN_PROTEIN_LEN=${MIN_PROTEIN_LEN:-5}
MODE=${MODE:-0}
SORT_STRATEGY=${SORT_STRATEGY:-1}
MIN_OVERLAP_LEN=${MIN_OVERLAP_LEN:-0.7}
SHORTNESS=${SHORTNESS:-0}

#------------------------------------------------------------------------------#
# Run cd-hit
#------------------------------------------------------------------------------#

echo "Starting script."

# Make directories if necessary
mkdir -p $(dirname $CLSTR_FILE)
mkdir -p $(dirname $REPRESENTATIVES)

# Run cd-hit
cd-hit \
-T 0 \
-n $WORD_SIZE \
-p 1 \
-c $MIN_SEQ_ID \
-d 0 \
-M $MEMORY \
-l $MIN_PROTEIN_LEN \
-g $MODE \
-sc $SORT_STRATEGY \
-aL $MIN_OVERLAP_LEN \
-s $SHORTNESS \
-i $INFILE \
-o $(basename ${INFILE})_cd-hit

# Rename the output files...
mv $(basename ${INFILE})_cd-hit.clstr $CLSTR_FILE
mv $(basename ${INFILE})_cd-hit $REPRESENTATIVES

echo "Finished."
