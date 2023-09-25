#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process  align {
  input:
  publishDir "${params.out_dir}"
    path sample_dir
  
  when:
  (saples.d== 10)
  script:
  """
  STAR --runThreadN ${params.n_thread} --runMode genomeGenerate --genomeDir ${params.ref_dir} --genomeFastaFiles ${sample_dir}/${sample_dir}.fasta
  """
}

process  align_paired {
  input:
  publishDir "${params.out_dir}"
    path sample_dir
  script:
  """
  STAR --runThreadN ${params.n_thread} --runMode genomeGenerate --genomeDir ${params.ref_dir} --genomeFastaFiles ${sample_dir}/${sample_dir}.fasta
  """
}

workflow {
  Channel.fromPath("${params.samples_dir}/*", type: 'dir') |
  align
}