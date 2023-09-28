params.in = "$projectDir/samples_fq/sample_1.fq"
out_dir = "$projectDir/results/fastqc/"


process fastqcQualityControl {
  publishDir "$projectDir/results/fastqc/", mode: 'copy'

  input:
  val sample

  output:
  

  script:
  """
  fastqc -o $projectDir/results/fastqc/$sample
  """
}

workflow {
    fastqcQualityControl(params.in)
}