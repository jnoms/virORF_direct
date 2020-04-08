//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.prodigal_mode = "single"
params.prodigal_genetic_code = 11 //standard with nTG except no CTG
params.prodigal_compress = "TRUE"
params.prodigal_closed_edges = 0 //not closed edges

//============================================================================//
// Define process
//============================================================================//
process prodigal {
  tag "$sampleID"
  publishDir "$params.out_dir/prodigal", mode: "copy"

  input:
  tuple sampleID, in_fasta

  output:
  tuple sampleID, file("${sampleID}_prodigal.txt*")

  script:
  """
  # Determine if prefix should be .gz
  if [[ ${params.prodigal_compress} == "TRUE" ]] ; then
      PREFIX=prodigal.txt.gz
  else
      PREFIX=prodigal.txt
  fi

  $workflow.projectDir/bin/bash/prodigal.sh \
  -i ${in_fasta} \
  -o ${sampleID}_\$PREFIX \
  -p ${params.prodigal_mode} \
  -g ${params.prodigal_genetic_code} \
  -c ${params.prodigal_compress} \
  -e ${params.prodigal_closed_edges}
  """
}
