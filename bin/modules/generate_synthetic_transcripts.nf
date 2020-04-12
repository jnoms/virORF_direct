//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.ref_fasta = "$baseDir/resources/sars_cov2_NC_045512.2_genome.fasta"
params.min_junction_length = 1000

//============================================================================//
// Define process
//============================================================================//
process generate_synthetic_transcripts {
  tag "$sampleID"
  publishDir "$params.out_dir/synthetic_transcripts", mode: "copy"

  input:
  tuple sampleID, bam, fasta

  output:
  tuple sampleID, file("${sampleID}_junctions.tsv"), file("${sampleID}_synthetic.fasta")

  script:
  """
  python $workflow.projectDir/bin/python/generate_synthetic_transcripts.py \
  -g ${params.ref_fasta} \
  -b ${bam} \
  -o ${sampleID}_synthetic.fasta \
  -j ${sampleID}_junctions.tsv \
  -l ${params.min_junction_length}
  """
}
