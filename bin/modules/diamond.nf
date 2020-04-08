//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.diamond_database = "$baseDir/resources/sars_cov2_cannonical_and_virORFsense.dmnd"
params.diamond_outfmt = "6 qseqid qlen sseqid evalue bitscore pident length qstart qend sstart send"
params.diamond_temp_dir = "temp"
params.diamond_evalue = "10"

//============================================================================//
// Define process
//============================================================================//
process diamond {
  tag "$sampleID"
  publishDir "$params.out_dir/diamond", mode: "copy"

  // default memory
  memory "10G"

  input:
  tuple sampleID, sequences

  output:
  tuple sampleID, file("*_diamond.out")

  script:
  """
  $workflow.projectDir/bin/bash/diamond.sh \
  -d ${params.diamond_database} \
  -q ${sequences} \
  -o ${sampleID}_diamond.out \
  -m ${task.memory.toGiga()} \
  -t ${params.diamond_temp_dir} \
  -e ${params.diamond_evalue} \
  -f "${params.diamond_outfmt}" \
  -s ${sampleID}
  """
}
