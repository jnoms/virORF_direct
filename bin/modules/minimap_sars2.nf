//============================================================================//
// Default params
//============================================================================//
params.out_dir = 'output'
params.ref_fasta = "$baseDir/resources/sars_cov2_NC_045512.2_genome.fasta"

//============================================================================//
// Define process
//============================================================================//
process minimap_sars2 {
  tag "$sampleID"
  publishDir "$params.out_dir/minimap", mode: "copy"

  input:
  tuple sampleID, dRNAseq_reads

  output:
  tuple sampleID, file("${sampleID}_mapped.bam"), file("${sampleID}_mapped.fasta")

  script:
  """
  $workflow.projectDir/bin/bash/minimap_sars2.sh \
  -i ${dRNAseq_reads} \
  -r ${params.ref_fasta} \
  -t ${task.cpus} \
  -o ${sampleID}_mapped.bam

  samtools fasta ${sampleID}_mapped.bam > ${sampleID}_mapped.fasta
  """
}
