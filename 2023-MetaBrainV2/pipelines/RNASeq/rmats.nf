params.gtfAnnotationFile = "../static_files/gencode.v44.annotation.gtf"

process identifyAlternativeSplicingSitesrMATS {
  publishDir "$projectDir/results/rmats", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}/*.txt.gz"
  
  shell:
  '''
  # 1. Check if the BAM file is derived from single or paired end
  numFASTQInputFiles=$(samtools view -H !{sample} | \
  grep "user command line" | \
  awk -F"user command line:" '{ print $2}' | \
  grep -o ".fastq.gz\\s\\|.fq.gz\\s\\|.fastq\\s\\|.fqs\\s\\ " | wc -l)

  if [ "${numFASTQInputFiles}" -eq 1 ];
  then
      end="single"
  elif [ "${numFASTQInputFiles}" -eq 2 ];
  then
      end="paired"
  else
      exit 1;
  fi

  # 2. Check the read length
  readLength=$(samtools view !{sample} | awk '{print length($10)}' | head -1000 | sort -u)
  
  # 3. Create config file
  echo !sample > config.txt

  # 4. Run rMATS command
  python ~/rmats-turbo/rmats.py --b1 config.txt \
  --gtf !{params.gtfAnnotationFile} \
  --readLength ${readLength} \
  --od !{sample.baseName} \
  --tmp rmats_tmp \
  --task both \
  -t ${end}  \
  --statoff

  # 5. Gzip all output files
  gzip !{sample.baseName}/*.txt
  '''
}

workflow {
    def bamFiles = Channel.fromPath('../sorted_bam/*.bam')
    identifyAlternativeSplicingSitesrMATS(bamFiles)
}