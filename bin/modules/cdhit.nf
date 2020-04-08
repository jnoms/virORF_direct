//============================================================================//
// Default params
//============================================================================//
params.out_dir = 'output'
params.cdhit_pident = 0.9
params.cdhit_min_overlap = 0.97

//============================================================================//
// Define process
//============================================================================//
process cdhit {
  tag "$sampleID"
  publishDir "$params.out_dir/cdhit", mode: "copy"

  // default memory
  memory "10G"

  input:
  tuple sampleID, pr_orfs, nt_orfs

  output:
  tuple sampleID, file("${sampleID}_pr_clusters.txt"), file("${sampleID}_pr_reps.fasta")

  script:
  """
  $workflow.projectDir/bin/bash/cd-hit.sh \
  -i ${pr_orfs} \
  -o ${sampleID}_pr_clusters.txt \
  -r ${sampleID}_pr_reps.fasta \
  -c ${params.cdhit_pident} \
  -a ${params.cdhit_min_overlap} \
  -M ${task.memory.toMega()}
  """
}
