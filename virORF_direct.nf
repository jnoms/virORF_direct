#!/usr/bin/env nextflow

nextflow.preview.dsl=2

// Log file must be set here to access laucnhDir information
params.log_file = "${workflow.launchDir}/${params.out_dir}/reports/virID.log"

// Require nextflow version 19.10 or higher
// if( !nextflow.version.matches('20.00+') ) {
//     println "This workflow requires Nextflow version 19.10 or greater -- You \
//     are running version $nextflow.version"
//     exit 1
// }

//============================================================================//
// Set up modules
//============================================================================//
include './bin/modules/minimap_sars2'
include './bin/modules/blast_leader'
include './bin/modules/fmlrc'
include './bin/modules/prodigal'
include './bin/modules/prodigal_to_orfs_direct'
include './bin/modules/cdhit'
include './bin/modules/annotate_cdhit_representatives'
include './bin/modules/write_cdhit_clusters'
include './bin/modules/diamond'
include './bin/modules/generate_synthetic_transcripts'

//============================================================================//
// Defining functions
//============================================================================//
def sampleID_set_from_infile(input) {
  // The purpose of this function is to take an input list of file paths and
  // to return a channel of sets of structure basename:file_path.
  sample_set = []
  for (String item : file(input)) {
    file = file(item)
    name = file.baseName
    sample_set.add([name, file])
  }
  ch = Channel.from(tuple(sample_set))
  return ch
}

//============================================================================//
// Define workflows
//============================================================================//
workflow virORF_direct_synthetic {

  get: input_ch
  main:

  // Mapping reads to reference fasta
  minimap_sars2(input_ch)
    .filter{ it[2].size() > 0 }
    .set{ minimap_sars2_result }

  // Generate synthetic reads
  generate_synthetic_transcripts(minimap_sars2_result)

  // Only keep reads that contain the leader sequence
  blast_leader(generate_synthetic_transcripts.out)
    .filter{ it[1].size() > 0 }
    .set{ blast_leader_result }

  // Predict ORFs with prodigal
  prodigal(blast_leader.out) \

  // Generate ORFs
  blast_leader.out
    .join(prodigal.out) \
    | prodigal_to_orfs_direct

  // Run CD-hit and label the representatives
  cdhit(prodigal_to_orfs_direct.out) \
    | annotate_cdhit_representatives

  // Write out the clusters into separate fastas
  prodigal_to_orfs_direct.out
    .join(cdhit.out) \
    | write_cdhit_clusters

  // Run cluster representatives through diamond
  diamond(annotate_cdhit_representatives.out)
}

workflow virORF_direct_correction {

  get: input_ch
  main:

  // Mapping reads to reference fasta
  minimap_sars2(input_ch)
    .filter{ it[2].size() > 0 }
    .set{ minimap_sars2_result }

  // Only keep reads that contain the leader sequence
  blast_leader(minimap_sars2_result)
    .filter{ it[1].size() > 0 }
    .set{ blast_leader_result }

  // Read in short reads as a channel, and mix them with blast_leader_result
  short_reads = Channel.fromPath(params.short_reads)
  blast_leader_result
    .combine(short_reads)
    .set{ long_and_short_reads }

  // Correct reads with fmlrc
  fmlrc(long_and_short_reads)

  // Predict ORFs with prodigal
  prodigal(fmlrc.out) \

  // Generate ORFs
  fmlrc.out
    .join(prodigal.out) \
    | prodigal_to_orfs_direct

  // Run CD-hit and label the representatives
  cdhit(prodigal_to_orfs_direct.out) \
    | annotate_cdhit_representatives

  // Write out the clusters into separate fastas
  prodigal_to_orfs_direct.out
    .join(cdhit.out) \
    | write_cdhit_clusters

  // Run cluster representatives through diamond
  diamond(annotate_cdhit_representatives.out)
}

//============================================================================//
// Validate inputs
//============================================================================//
if( (params.mode != "synthetic") && (params.mode != "correction") ) {
  error "params.mode must be set to either 'synthetic' or 'correction'.\
  Current set to ${params.mode}."
}

if( params.diamond_database == "" ) {
  error "params.diamond_database is a required parameter, please enter a path."
}

//============================================================================//
// Define main workflow
//============================================================================//
workflow {

  main:
    if ( params.mode == "synthetic" ) {
      input_ch = sampleID_set_from_infile(params.dRNAseq_reads)
      virORF_direct_synthetic(input_ch)
    }
    else if ( params.mode == "direct" ) {
     input_ch = sampleID_set_from_infile(params.dRNAseq_reads)
     virORF_direct_correction(input_ch)
    }
}
