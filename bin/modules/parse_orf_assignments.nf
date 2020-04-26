//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.parse_orf_assignments_diamond_fields = "qseqid qlen sseqid evalue bitscore pident length qstart qend sstart send"
params.parse_orf_assignments_min_pident = "95"
params.parse_orf_assignments_min_ok_nterm_extension = "20"

//============================================================================//
// Define process
//============================================================================//
process parse_orf_assignments {
  tag "$sampleID"
  publishDir "$params.out_dir/parse_orf_assignments", mode: "copy"

  input:
  tuple sampleID, diamond_file

  output:
  tuple sampleID, file("${sampleID}*")

  script:
  """
  python $workflow.projectDir/bin/python/parse_orf_assignments.py \
  -d ${diamond_file} \
  -o ${sampleID} \
  -f "${params.parse_orf_assignments_diamond_fields}" \
  -p ${params.parse_orf_assignments_min_pident} \
  -n ${params.parse_orf_assignments_min_ok_nterm_extension}
  """
}
