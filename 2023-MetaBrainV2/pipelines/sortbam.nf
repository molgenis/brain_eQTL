process sortBAM {
  publishDir "$projectDir/results/sorted_bam", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}.sorted.bam"
  
  script:
  """
  samtools sort $sample -o ${sample.baseName}.sorted.bam
  """
}

workflow {
    def samples = Channel.fromPath('../samples/*.bam')
    sortBAM(samples)
}