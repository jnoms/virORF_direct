//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process prodigal_to_orfs_direct {
  tag "$sampleID"
  publishDir "$params.out_dir/prodigal", mode: "copy"

  input:
  tuple sampleID, in_fasta, prodigal_file

  output:
  tuple sampleID, file("${sampleID}_pr.fasta"), file("${sampleID}_nt.fasta")

  script:
  """
  python $workflow.projectDir/bin/python/prodigal_to_orfs_direct.py \
  -f ${in_fasta} \
  -p ${prodigal_file} \
  -N ${sampleID}_nt.fasta \
  -P ${sampleID}_pr.fasta
  """
}
