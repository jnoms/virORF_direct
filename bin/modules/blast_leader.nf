//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.leader_fa = "$baseDir/resources/SARS2_leader.fasta"
params.blast_leader_type = 'blastn'
params.blast_leader_min_alignment_length = 45
params.blast_leader_max_start_position = 50
params.blast_leader_min_pident = 0.85

//============================================================================//
// Define process
//============================================================================//
process blast_leader {
  tag "$sampleID"
  publishDir "$params.out_dir/blast_leader", mode: "copy"

  input:
  tuple sampleID, in_bam_or_junctions, in_fasta

  output:
  tuple sampleID, file("${sampleID}_leader.fasta")

  script:
  """
  $workflow.projectDir/bin/bash/blast_leader.sh \
  -q ${in_fasta} \
  -s ${params.leader_fa} \
  -o ${sampleID}_leader.fasta \
  -B ${params.blast_leader_type} \
  -t ${task.cpus} \
  -S ${sampleID} \
  -a ${params.blast_leader_min_alignment_length} \
  -b ${params.blast_leader_max_start_position} \
  -c ${params.blast_leader_min_pident}
  """
}
