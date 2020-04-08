//============================================================================//
// Default params
//============================================================================//
params.out_dir = 'output'

//============================================================================//
// Define process
//============================================================================//
process fmlrc {
  tag "$sampleID"
  publishDir "$params.out_dir/corrected_reads", mode: "copy"

  input:
  tuple sampleID, long_reads_fasta, short_reads

  output:
  tuple sampleID, file("${sampleID}_corrected.fasta")

  script:
  """
  $workflow.projectDir/bin/bash/fmlrc.sh \
  -l ${long_reads_fasta} \
  -s ${short_reads} \
  -t ${task.cpu} \
  -o ${sampleID}_corrected.fasta
  """
}
