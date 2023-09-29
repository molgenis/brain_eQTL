process convertBAMToFASTQ {
  publishDir "$projectDir/results/fastq", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}/*.fastq"
  
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

workflow {
    def bamFiles = Channel.fromPath('../sorted_bam/*.bam')
    convertBAMToFASTQ(bamFiles)
}