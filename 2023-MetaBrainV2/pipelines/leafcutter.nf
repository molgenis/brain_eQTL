process identifyAlternativeSplicingSitesLeafCutter {
  publishDir "$projectDir/results/leafcutter", mode: 'copy'

  input:
  path sample
  
  output:
  path "leafcutter/*.junc"
  
  shell:
  '''
  # 1. Index BAM file
  samtools index !{sample}

  # 2. Run regtools command
  regtools junctions extract -s XS -a 8 -m 50 -M 500000 !{sample} \
  -o leafcutter/!{sample.baseName}.junc
  '''
}

workflow {
    def bamFiles = Channel.fromPath('../sorted_bam/*.bam')
    identifyAlternativeSplicingSitesLeafCutter(bamFiles)
}