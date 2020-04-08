//============================================================================//
// Default params
//============================================================================//
params.out_dir = 'output'

//============================================================================//
// Define process
//============================================================================//
process annotate_cdhit_representatives {
  tag "$sampleID"
  publishDir "$params.out_dir/cdhit", mode: "copy"

  input:
  tuple sampleID, cluster_file, representative_fasta

  output:
  tuple sampleID, file("${sampleID}_pr_reps_annotated.fasta")

  script:
  """
  python $workflow.projectDir/bin/python/order_and_quant_cdhit_reps.py \
  -c ${cluster_file} \
  -f ${representative_fasta}  \
  -o ${sampleID}_pr_reps_annotated.fasta
  """
}
