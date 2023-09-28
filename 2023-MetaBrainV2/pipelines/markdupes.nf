#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Standard process without external user defined output paths.
process markduplicates {
  // SLURM Params.
  time '2h'
  memory '16 GB'
  cpus 1

  // Take an unsorted .BAM as input.
  input:
  path bam_file

  // All output files which sorting produces.
  output:
  path "${bam_file.getBaseName()}.bam"
  path "${bam_file.getBaseName()}_dupes.txt"

  // Differing scripts dependong on whether input is single/paired.
  script:
  """
  module load picard
  java -jar \${EBROOTPICARD}/picard.jar MarkDuplicates\
      I=${bam_file} \
      O=${bam_file.getBaseName()}.bam \
      M=${bam_file.getBaseName()}_dupes.txt
  """
}
