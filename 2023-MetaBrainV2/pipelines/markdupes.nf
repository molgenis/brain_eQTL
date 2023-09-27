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

  // All output files which STAR produces.
  output:

  // Differing scripts dependong on whether input is single/paired.
  script:
  """

  """
}

// Create process channel for every sample in a user defined directory.
workflow {
  Channel.fromPath("${params.samples_dir}/*", type: 'dir') |
  align
}