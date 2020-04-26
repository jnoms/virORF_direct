#!/usr/bin/env nextflow
nextflow.preview.dsl=2

//============================================================================//
// Set up modules
//============================================================================//
include './bin/modules/minimap_sars2' params(params)
include './bin/modules/blast_leader' params(params)
include './bin/modules/prodigal' params(params)
include './bin/modules/prodigal_to_orfs_direct' params(params)
include './bin/modules/diamond' params(params)
include './bin/modules/generate_synthetic_transcripts' params(params)
include './bin/modules/parse_orf_assignments' params(params)
include './bin/modules/find_trs_sites' params(params)

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

  // Run protein ORFs through diamond
  diamond(prodigal_to_orfs_direct.out)

  // Parse the diamond output file
  parse_orf_assignments(diamond.out)

  // Find TRS sequences from the junctions
  find_trs_sites(generate_synthetic_transcripts.out)
}

//============================================================================//
// Validate inputs
//============================================================================//
if( params.diamond_database == "" ) {
  error "params.diamond_database is a required parameter, please enter a path."
}

//============================================================================//
// Define main workflow
//============================================================================//
workflow {

  main:
    input_ch = sampleID_set_from_infile(params.dRNAseq_reads)
    virORF_direct_synthetic(input_ch)
}
