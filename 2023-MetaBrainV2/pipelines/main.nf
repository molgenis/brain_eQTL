params.refFlat = ""
params.referenceGenome = ""
params.gtfAnnotationFile = ""
params.sampleDir = ""

process convertBAMToFASTQ {
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "$projectDir/results/fastq", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}"
  
  shell:
  '''
  # 1. Determine if BAM file is single or paired end reads
  numFASTQFiles=$(samtools view -H !{sample} | \
  grep "@PG" | \
  tr ' ' '\n' | \
  grep -oE '(.fq.gz|.fastq.gz|.fq|.fastq)($)' | \
  wc -l)

  # 2. Create command arguments
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

process align {
  // SLURM Params.
  time '6h'
  memory '16 GB'
  cpus 1

  // Take sample directory as input.
  input:
  path sample_dir

  // All output files which STAR produces.
  output:
  path "${sample_dir}/${sample_dir}_Aligned.out.bam", emit: bamFile
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

process sortBAM {
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "$projectDir/results/sorted_bam", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}.sorted.bam"
  
  script:
  """
  samtools sort $sample -o ${sample.baseName}.sorted.bam
  """
}

// Standard process without external user defined output paths.
process markDuplicates {
  time '6h'
  memory '8 GB'
  cpus 1

  // Take an sorted .BAM as input.
  input:
  path bam_file

  // All output files which sorting produces.
  output:
  path "${bam_file.getBaseName()}.bam"
  // path "${bam_file.getBaseName()}_dupes.txt"

  script:
  """
  PicardCommandLine MarkDuplicates \
      I=${bam_file} \
      O=${bam_file.getBaseName()}.bam \
      M=${bam_file.getBaseName()}_dupes.txt
  """
}

process QCwithRNASeqMetrics {
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "$projectDir/results/rna_seq_metrics", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}_rnaseqmetrics"
  
  script:
  """
  PicardCommandLine CollectRnaSeqMetrics \
  I=${sample} \
  O=${sample.baseName}_rnaseqmetrics \
  REF_FLAT=${params.refFlat} \
  STRAND=NONE
  """
}

process QCwithMultipleMetrics {
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "$projectDir/results/multiple_metrics", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}/*"
  
  script:
  """
  mkdir ${sample.baseName}
  
  PicardCommandLine CollectMultipleMetrics I=${sample} \
  O=${sample.baseName}/multiple_metrics \
  R=${params.referenceGenome} \
  PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=QualityScoreDistribution \
  PROGRAM=MeanQualityByCycle \
  PROGRAM=CollectInsertSizeMetrics
  """
}

process identifyAlternativeSplicingSitesrMATS {
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "$projectDir/results/rmats", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}/*.txt.gz"
  
  shell:
  '''
  # 1. Check if the BAM file is derived from single or paired end reads
  numFASTQFiles=$(samtools view -H !{sample} | \
  grep "@PG" | \
  tr ' ' '\n' | \
  grep -oE '(.fq.gz|.fastq.gz|.fq|.fastq)($)' | \
  wc -l)

  if [ "${numFASTQFiles}" -eq 1 ]; 
  then
    end="single"
  elif [ "${numFASTQFiles}" -gt 1 ];
  then
    end="paired"
  else
    exit 1
  fi

  # 2. Check the read length
  readLength=$(samtools view !{sample} |head -n1 |awk '{print $10}'|tr -d "\\n" |wc -m)

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

process identifyAlternativeSplicingSitesLeafCutter {
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "$projectDir/results/leafcutter", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}/*.junc"
  
  shell:
  '''
  # 1. Index BAM file
  samtools index !{sample}

  # 2. Run regtools command
  mkdir !{sample.baseName}
  regtools junctions extract -s XS -a 8 -m 50 -M 500000 !{sample} -o !{sample.baseName}/!{sample.baseName}.junc 
  '''
}

process convertBAMToCRAM {
  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "$projectDir/results/cram", mode: 'copy'

  input:
  path sample
  
  output:
  path "${sample.baseName}.cram"
  
  script:
  """
  samtools view -T ${params.referenceGenome} -C -o ${sample.baseName}.cram ${sample}
  """
}

workflow {
    def bamFiles = Channel.fromPath(params.sampleDir)
    convertBAMToFASTQ(bamFiles)
    fastqcQualityControl(convertBAMToFASTQ.out)
    align(convertBAMToFASTQ.out)
    sortBAM(align.out.bamFile)
    markDuplicates(sortBAM.out)
    QCwithRNASeqMetrics(markDuplicates.out)
    QCwithMultipleMetrics(markDuplicates.out)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out)
    convertBAMToCRAM(markDuplicates.out)
}