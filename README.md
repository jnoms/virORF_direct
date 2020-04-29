# virORF_direct
virORF_direct: Protein-level analysis of direct RNA sequencing data.

## Description
This repository and pipeline is sufficient to reproduce the work reported in the manuscript "Noncanonical junctions in subgenomic RNAs of SARS-CoV-2 lead to variant open reading frames."

The purpose of this pipeline is to process direct RNA sequencing (dRNAseq) reads from SARS-CoV-2-infected cells. This script maps the reads against the SARS-CoV-2 genome and generates coordinate-derived transcripts and junction coordinates from these alignments. 

Coordinate-derived transcripts are used to characterize the viral open reading frames (ORFs) present in the input dRNAseq transcriptome. ORFs are predicted from each coordinate-derived transcript using prodigal, prodigal ORFs are processed to generate protein sequences, and these protein sequences are mapped against annotated and predicted SARS-CoV-2 proteins using the DIAMOND aligner. This pipeline then incorporates a script to process the resultant DIAMOND alignments into a series of output files which can be visualized as desired. In addition, this pipeline predicts common TRS sequences from the junction coordinates generated at the coordinate-derived transcripts generation step.

R-markdown files used to generate the figures present in the manuscript "Noncanonical junctions in subgenomic RNAs of SARS-CoV-2 lead to variant open reading frames" are present in bin/R and are not run automatically.

All python and bash scripts can be run independently. Execute any python script (with -h) or bash script (without arguments) for detailed instructions on how to run them.

This pipeline automatically uses the Conda package manager to handle all software depencies besides Nextflow and Conda itself. Any process that can be run in parallel will automatically be run in parallel by this pipeline. As long as the "executor" is configured in the nextflow.config file, Nextflow will also automatically handle submitting jobs to the executor (i.e. slurm). 

## Quickstart
1. Download Nextflow version 19.10 or higher. Nextflow can be downloaded by running ```curl -s https://get.nextflow.io | bash``` , which makes the nextflow executible in your working directory. See https://www.nextflow.io/ for more information.
2. If necessary, download Conda. Miniconda can be downloaded [here](https://docs.conda.io/en/latest/miniconda.html). **This pipeline handles all software dependencies besides Nextflow and Conda itself automatically using a conda enviornment.**
2. Enter essential parameters to nextflow.config:
    - dRNAseq_reads: Path to the direct RNAseq reads, in fasta or fastq format.
    - You may want to alter other paremeters (described below), but there are sensible defaults in nextflow.config.
3. Select or configure profile in nextflow.config:
    - The profile o2_slurm is set up to run virORF_direct in a linux environment using the slurm scheduler.
    - To change the scheduler, alter the "executor" to the desired scheduler.
    - You may also need to change the `beforeScript` depending on your cluster requirements. It is currently set to load conda at my institution's cluster. 
    - The profile o2_local runs this script in a linux environment locally.
    - The profile mac runs this script locally on a mac.
    - To run virORF_direct with the selected profile, will use the `-profile` switch upon running this pipeline.
    - Finally, you also need to change the Conda `cacheDir` in the profile you're using to be a path to the directory you want the Conda environment to be saved to. This way, the Conda environment will only need to be downloaded once. You can also exclude this switch, in which case the Conda environment will be downloaded to your working directory.  
4. Run this pipeline. Switches can be entered on the command line and/or set in nextflow config. To run:
```
nextflow run virORF_direct.nf \
--dRNAseq_reads "path_to_reads/*fastq" \
--out_dir "output" \
-profile o2_slurm \
-resume
```
The `-resume` switch lets you resume this pipeline from the last finished step. The `-profile` switch should be set to whatever profile you decided was appropriate based on your machine or cluster. The SARS-CoV-2 reference genome, leader sequence, and diamond database are provided in this repository and coded into the nextflow.config (or you can use your own), as are conda yml files.

## Description of Parameters present in nextflow.config.

### General parameters  
**dRNAseq_reads:** Path to the direct RNAseq reads in **fasta** or **fastq** format. Even if you only have one input file, this should include at least one glob! `input/*fastq`  
**out_dir:** The directory that will contain the output files.  
**ref_fasta:** Path to the reference fasta containing the SARS-CoV-2 reference genome. This is included with this repository at resources/sars_cov2_NC_045512.2_genome.fasta and is programmed in to the config file. `$baseDir/resources/sars_cov2_NC_045512.2_genome.fasta`  

### Synthetic read generation  
**min_junction_length:** The minimum length of a minimap2-called deletion or intron for the region to be considered a junction. `1000`  

### Checking for presence of the leader sequence  
**leader_fa:** Path to the leader sequence fasta. This is provided in this repository at resources/sars_cov2_NC_045512.2_genome.fasta and is coded into nextflow.config. `$baseDir/resources/SARS2_leader.fasta`  
**blast_leader_type:** The method of blast to use. `blastn`  
**blast_leader_min_alignment_length:** Minimum alignment length between the query read and the leader sequence. `45`  
**blast_leader_max_start_position:** Maximum location on the read the leader alignment can start. `50`  
**blast_leader_min_pident:** The minimum percent identity between the aligned segment of the read and the leader sequence `0.85`  

### ORF calling with Prodigal  
**prodigal_mode:** The mode of prodigal to use. Recommend sticking with default. `single`  
**prodigal_genetic_code:** The NCBI genetic code to use. Recommend sticking with 11, as this is like the standard code (1) but includes alternative nTG start codons. `11`  
**prodigal_compress:** Automatically gzips the prodigal output because it can be quite large. Downstream processes can handle it either way. `TRUE`  
**prodigal_closed_edges:** If set to 0, allows ORFs to "run in" off the edge - i.e. the ORF doesn't need a start codon if it's close to the start of the sequence. If 1, it requires a start codon. `1`  

### DIAMOND alignment  
**diamond_database:** Path to the diamond database. A diamond database containing canonical and predicted SARS-CoV-2 proteins is included with this repository. `$baseDir/resources/sars_cov2_cannonical_and_virORFsense_ORF1_split.dmnd`  
**diamond_outfmt:** The output format to use when calling diamond. You can add output fields on, but probably removing them might break things. `"6 qseqid qlen sseqid evalue bitscore pident length qstart qend sstart send"`  
**diamond_temp_dir:** Temporary directory for diamond to use. `temp`  
**diamond_evalue:** The minimum evalue for an alignment for diamond to report out. I set this at 10 because downstream scripts filter at percent identity instead. `10`   

### Parsing the DIAMOND output file  
**parse_orf_assignments_diamond_fields:** Should be the same as diamond_outfmt minus the preceeding 6. `"qseqid qlen sseqid evalue bitscore pident length qstart qend sstart send"`  
**parse_orf_assignments_min_pident:** Minimum percent identity to keep an alingment. Should be entered as a full interger `95`  
**parse_orf_assignments_min_ok_nterm_extension:** The allowable amount of additional sequence allowed at the N terminus of an alignment for it to still be called canonical. This addresses situations where prodigal calls a truly canonical ORF from an upstream nTG, while the DIAMOND database has it form it's downstream ATG, and we don't want to call it variant by mistake. A setting of 20 allows 20 additional amino acids to the N-terminus of the query to still consider the query canonical. `20`  

### Predicting TRS sequences  
**TRS_proximal_sequence:** Amount of sequence to consider on either side of each junction point. For example, an entry of 30 means it will look for 15 bases on either side of the 5' and 3' junctions, and compare those two 30 base sequences to find the homologous sequence. This must be even. `30`  
**TRS_top:** The top N TRS' for each category to output. `2`  

## Note on Conda
There is a linux and mac conda yml programmed into the config file and provided in this repository. The linux yml is used for o2_local and o2_slurm profiles, while the mac yml is used for the mac profile. If for some reason the enviornment specified by either of these ymls does not work on your machine, you can generate a machine-specific conda yml by running
```
conda create --name virORF_direct jupyter biopython prodigal pathlib pandas blast diamond minimap2 samtools ete3 pysam seqtk
conda activate virORF_direct
conda env export > virORF_direct.yml
conda deactivate
```
And then specify the path to virORF_direct.yml in the "conda" process setting of nextflow.config.

## Description of tab-delimited output files  
I'll only highlight the tab-delimited files which can be plotted. The outputs after each step are by default transfered to the output directory - this includes the minimap2 bam, the coordinate-derived transcripts, and so on.  

### generate_synthetic_transcripts.py  
**${sampleID}_junctions.tsv:** Each line is a junction. Has the read_name, the 5' junction location, and the 3' junction location. In the rare cases that there are multiple junctions for a single transcript, they'll end up on different lines.  

### find_trs_sites.py  
**${sampleID}_TRS_sites.tsv:**  Has the top N homology (TRS) sequences between the 5' and 3' regions of each junction by category. More detail is present in the manuscript methods, and in the find_trs_sites.py -h information.  

### parse_orf_assignments.py  
**${sampleID}_subject_starts_and_ends.tsv:** This contains the predicted start and stop positions of each query ORF based on their alignment starts and stops against cannonical proteins. So, this is  *protein level* information of all variant and canonical ORFs. This output file also gives information on any ORF fusions, if applicable.  

**${sampleID}_total_ORF_counts.tsv:** This file contains total counts for a given ORF. This does NOT consider whether the ORF is canonical or variant. This DOES include all unannotated ORFs.  

**${sampleID}_fiveprime-most-orf_counts.tsv:** This file lists the number of transcripts which contain each cannonical and variant ORF as the most 5' ORF. Does include unannotated ORFs.  

**${sampleID}_cannonical_vs_variant_counts.tsv:** This file contains total counts for variant vs canonical ORFs. This does NOT include unannotated ORFs.  

**${sampleID}_cannonical_vs_variant_startsites.tsv:** This file contains transcript ORF start and end sites for variant vs canonical ORFs. So, this is the nucleotide start and end points for each ORF on each transcript. This does not include unannotated ORFs.  

