#!/bin/bash

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose of this script is to take an input fasta and generate a
        prodigal output file.

        USAGE: $0

        Required:
        -i <INFILE>
            Path to the input fasta.
        -o <OUTFILE>
            Path the the tabular prodigal output file.

        Optional:
        -p <MODE> [single]
            Either 'single' or 'meta'. single is default because, if meta, it
            sometimes uses organisms that have a genetic code that isn't 1 and
            then overrides my genetic code input... Single makes the assumption
            all sequences are from the same organism. It also might not work if
            input size is too small.
        -g <GENETIC_CODE> [11]
            Defaults to 11, which is similar to standard (1) but also include
            nTG start codons EXCEPT CTG (which is a huge drag).
        -c <COMPRESS> [FALSE]
            If set to 'TRUE', will gzip output.
            This is helpful because prodigal output files can be huge when
            processing direct RNAseq reads. Note that the output file name
            will still be OUTFILE, so recommend adding a .gz to the end of it.
        -e <CLOSED_ENDS> [0]
            If set to 1, will not consider ORFs that start off the edge. If set
            to 0, an ORF could 'run in' from before the sequence started.
        "
}

#If less than 2 options are input, show usage and exit script.
if [ $# -le 2 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts i:o:p:g:c:e: option ; do
        case "${option}"
        in
                i) INFILE=${OPTARG};;
                o) OUTFILE=${OPTARG};;

                p) MODE=${OPTARG};;
                g) GENETIC_CODE=${OPTARG};;
                c) COMPRESS=${OPTARG};;
                e) CLOSED_ENDS=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Set default parameters
#------------------------------------------------------------------------------#
MODE=${MODE:-single}
GENETIC_CODE=${GENETIC_CODE:-11}
COMPRESS=${COMPRESS:-FALSE}
CLOSED_ENDS=${CLOSED_ENDS:-0}

#------------------------------------------------------------------------------#
# Run prodigal
#------------------------------------------------------------------------------#
echo "Starting prodigal script. Processing the infile $INFILE."

# Make outdir if needed
mkdir -p $(dirname $OUTFILE)

# Run prodigal
if [[ $CLOSED_ENDS == 0 ]] ; then
  prodigal \
  -p $MODE \
  -i $INFILE \
  -s ${OUTFILE}.tmp \
  -g $GENETIC_CODE
else
  prodigal \
  -p $MODE \
  -i $INFILE \
  -s ${OUTFILE}.tmp \
  -g $GENETIC_CODE \
  -c
fi

# Compress if desired
if [[ $COMPRESS == "TRUE" ]] ; then
  gzip ${OUTFILE}.tmp && mv ${OUTFILE}.tmp.gz ${OUTFILE}
else
  mv ${OUTFILE}.tmp ${OUTFILE}
fi

echo "Finished."
