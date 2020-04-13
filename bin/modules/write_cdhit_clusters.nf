//============================================================================//
// Default params
//============================================================================//
params.out_dir = 'output'

//============================================================================//
// Define process
//============================================================================//
process write_cdhit_clusters {
  tag "$sampleID"
  publishDir "$params.out_dir/cdhit/${sampleID}_clusters", mode: "copy"

  input:
  tuple sampleID, pr_orfs, nt_orfs, cluster_file, representative_fasta

  output:
  tuple sampleID, file("pr_clusters/${sampleID}_pr*"), file("nt_clusters/${sampleID}_nt*")

  script:
  """
  python $workflow.projectDir/bin/python/write_cdhit_clusters.py \
  -c ${cluster_file}  \
  -f ${pr_orfs} \
  -o pr_clusters/${sampleID}_pr

  python $workflow.projectDir/bin/python/write_cdhit_clusters.py \
  -c ${cluster_file}  \
  -f ${nt_orfs} \
  -o nt_clusters/${sampleID}_nt
  """
}
