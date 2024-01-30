#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Standard process without external user defined output paths.
process align {
  // SLURM Params.
  time '2h'
  memory '16 GB'
  cpus 1

  // Take sample directory as input.
  input:
  path sample_dir

  // All output files which STAR produces.
  output:
  path "${sample_dir}/${sample_dir}_Aligned.out.bam"
  path "${sample_dir}/${sample_dir}_Log.final.out"
  path "${sample_dir}/${sample_dir}_Log.progress.out"
  path "${sample_dir}/${sample_dir}_SJ.out.tab"

  // Differing scripts dependong on whether input is single/paired.
  script:
  """
  module load STAR
  if [[ -f ${sample_dir}/${sample_dir}_2.fastq ]]; then
    STAR --runThreadN 1 --outFileNamePrefix ${sample_dir}/${sample_dir}_ --outSAMtype BAM Unsorted \
    --genomeDir ${params.ref_dir} --readFilesIn ${params.samples_dir}/${sample_dir}/${sample_dir}_1.fastq ${params.samples_dir}/${sample_dir}/${sample_dir}_2.fastq
  else
    STAR --runThreadN 1 --outFileNamePrefix ${sample_dir}/${sample_dir}_ --outSAMtype BAM Unsorted \
    --genomeDir ${params.ref_dir} --readFilesIn ${params.samples_dir}/${sample_dir}/${sample_dir}.fastq
  fi  
  """
}

// Standard process with user defined output paths.
process align_paths {
  time '2h'
  memory '16 GB'
  cpus 1
  input:
  publishDir "${params.out_dir}"
  path sample_dir
  output:
  path "${params.out_dir}/${sample_dir}/${sample_dir}"

  script:
  """
  mkdir ${params.out_dir}/${sample_dir}
  module load STAR
  if [[ -f ${sample_dir}/${sample_dir}_2.fastq ]]; then
    STAR --runThreadN 1 --outFileNamePrefix ${params.out_dir}/${sample_dir}/${sample_dir} --outSAMtype BAM Unsorted \
    --genomeDir ${params.ref_dir} --readFilesIn ${params.samples_dir}/${sample_dir}/${sample_dir}_1.fastq ${params.samples_dir}/${sample_dir}/${sample_dir}_2.fastq
  else
    STAR --runThreadN 1 --outFileNamePrefix ${params.out_dir}/${sample_dir}/${sample_dir} --outSAMtype BAM Unsorted \
    --genomeDir ${params.ref_dir} --readFilesIn ${params.samples_dir}/${sample_dir}/${sample_dir}.fastq
  fi  
  """
}

// Create process channel for every sample in a user defined directory.
workflow {
  Channel.fromPath("${params.samples_dir}/*", type: 'dir') |
  align
}