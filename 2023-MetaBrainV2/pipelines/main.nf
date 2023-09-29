params.refFlat = ""
params.referenceGenome = ""
params.gtfAnnotationFile = ""
params.sampleDir = ""

process convertBAMToFASTQ {
  publishDir "$projectDir/results/fastq", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}"
  
  shell:
  '''
  # 1. Determine if BAM file is single or paired end
  numFASTQFiles=$(samtools view -H !{sample} | \
  grep "@PG" | \
  tr ' ' '\n' | \
  grep -oE '(.fq.gz|.fastq.gz|.fq|.fastq)($)' | \
  wc -l)

  # 2. Create command
  if [ "${numFASTQFiles}" -eq 1 ]; 
  then
    arguments="fastq -1 !{sample.baseName}/!{sample.baseName}.fastq !{sample}"
  elif [ "${numFASTQFiles}" -gt 1 ];
  then
    arguments="fastq -1 !{sample.baseName}/!{sample.baseName}_1.fastq -2 !{sample.baseName}/!{sample.baseName}_2.fastq !{sample}"
  else
    exit 1
  fi

  # 3. Make sample directory
  mkdir !{sample.baseName}

  # 4. Run command
  samtools ${arguments} > /dev/null 2>&1
  '''
}

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

workflow {
    def bamFiles = Channel.fromPath(params.sampleDir)
    convertBAMToFASTQ(bamFiles)
    align(convertBAMToFASTQ.out)
}