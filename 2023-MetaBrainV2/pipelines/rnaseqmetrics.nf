params.refFlat = "../static_files/refFlat.txt"

process QCwithRNASeqMetrics {
  publishDir "$projectDir/results/rna_seq_metrics", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}_rnaseqmetrics"
  
  script:
  """
  PicardCommandLine CollectRnaSeqMetrics \
  I=${sample} \
  O=${sample.baseName}_rnaseqmetrics \
  REF_FLAT=${params.refFlat} \
  STRAND=NONE
  """
}

workflow {
    def bamFiles = Channel.fromPath('../sorted_bam/*.bam')
    QCwithRNASeqMetrics(bamFiles)
}