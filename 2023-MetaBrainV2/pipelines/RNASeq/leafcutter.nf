process IdentifyAlternativeSplicingSitesLeafCutter {
  publishDir "$projectDir/results/leafcutter", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}/*.junc"
  
  shell:
  '''
  # 1. Index BAM file
  samtools index !{sample}

  # 2. Run regtools command
  mkdir !{sample.baseName}
  regtools junctions extract -s XS -a 8 -m 50 -M 500000 !{sample} -o !{sample.baseName}/!{sample.baseName}.junc 
  '''
}

workflow {
    def bamFiles = Channel.fromPath('../sorted_bam/*.bam')
    identifyAlternativeSplicingSitesLeafCutter(bamFiles)
}