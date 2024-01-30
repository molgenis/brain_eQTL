params.referenceGenome = "../static_files/GRCh38_latest_genom.fna"

process convertBAMToCRAM {
  publishDir "$projectDir/results/cram", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}.cram"
  
  script:
  """
  samtools view -T ${params.referenceGenome} -C -o ${sample.baseName}.cram ${sample}
  """
}

workflow {
  def bamFiles = Channel.fromPath('../sorted_bam/*.bam')
  convertBAMToCRAM(bamFiles)
}