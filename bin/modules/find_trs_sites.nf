//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.ref_fasta = "$baseDir/resources/sars_cov2_NC_045512.2_genome.fasta"
params.TRS_proximal_sequence = "30"
params.TRS_top = "2"

//============================================================================//
// Define process
//============================================================================//
process find_trs_sites {
  tag "$sampleID"
  publishDir "$params.out_dir/find_trs_sites", mode: "copy"

  input:
  tuple sampleID, junctions, synthetic_transcripts

  output:
  tuple sampleID, file("${sampleID}_TRS_sites.tsv")

  script:
  """
  python $workflow.projectDir/bin/python/find_trs_sites.py \
  -j ${junctions} \
  -o ${sampleID}_TRS_sites.tsv \
  -r ${params.ref_fasta} \
  -p ${params.TRS_proximal_sequence} \
  -N ${params.TRS_top}
  """
}
