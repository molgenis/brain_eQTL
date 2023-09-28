params.referenceGenome = "../static_files/GRCh38_latest_genomic.fna"


process QCwithMultipleMetrics {
  publishDir "$projectDir/results/multiple_metrics", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}/*"
  
  script:
  """
  mkdir ${sample.baseName}
  
  PicardCommandLine CollectMultipleMetrics I=${sample} \
  O=${sample.baseName}/multiple_metrics \
  R=${params.referenceGenome}
  """
}

workflow {
    def bamFiles = Channel.fromPath('../sorted_bam/*.bam')
    QCwithMultipleMetrics(bamFiles)
}