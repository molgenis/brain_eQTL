process fastqcQualityControl {
  publishDir "$projectDir/results/fastqc/", mode: 'copy'

  input:
  val sample

  output:
  file "*_fastqc.zip"
  

  script:
  """
  mkdir fastqc
  fastqc ${sample} -o .
  """
}

workflow {
  inputFiles = Channel.fromPath("samples_fq/*.fastq")
  .concat(Channel.fromPath("samples_fq/*fq"))
  .concat(Channel.fromPath("samples_fq/*.fastq.gz"))
  .concat(Channel.fromPath("samples_fq/*.fq.gz"))

  inputFiles.view()  
  fastqcQualityControl(inputFiles)
}