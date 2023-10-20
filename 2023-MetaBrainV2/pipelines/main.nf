nextflow.enable.dsl=2

process convertBAMToFASTQ {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1
  
  input:
  val sample
  
  output:
  path "fastq_output"
  
  shell:
  '''
  # Split the input string by ';' to separate file paths
  IFS=';' read -r -a sampleFiles <<< "!{sample}"

  # Get the extension of the sample
  extension=$(awk -F '.' '{print $NF}' <<< "${sampleFiles[0]}")

  # Split the file path by '/' to and get the last part (this is the sample name)
  IFS='/' read -r -a directories <<< "${sampleFiles[0]}"
  fileName="${directories[-1]}"
  sampleName=${fileName%".$extension"}

  # Clean the path (remove invisble characters)
  cleaned_path=$(echo "${sampleFiles[0]}" | tr -d '[[:cntrl:]]')
  
  # Check if the sample is a .bam file
  if [[ "$extension" == *"bam" ]]; then
    sample=$cleaned_path

    # 1. Determine if BAM file is single or paired end reads
    numFASTQFiles=$(samtools view -H $cleaned_path | \
    grep "@PG" | \
    tr ' ' '\n' | \
    grep -oE '(.fq.gz|.fastq.gz|.fq|.fastq)($)' | \
    wc -l)

    # 2. Sort the BAM file by name
    samtools sort -n ${sample} -o "sorted_${sampleName}.bam"

    # 3. Create command arguments
    if [[ "${numFASTQFiles}" -eq 1 ]]; 
    then
      arguments="fastq -0 fastq_output/${sampleName}.fastq.gz -n sorted_${sampleName}.bam"
    elif [ "${numFASTQFiles}" -gt 1 ];
    then
      arguments="fastq -1 fastq_output/${sampleName}_1.fastq.gz -2 fastq_output/${sampleName}_2.fastq.gz -n sorted_${sampleName}.bam"
    else
      exit 1
    fi

    # 4. Make sample directory
    mkdir fastq_output

    # 5. Run command
    samtools ${arguments}

  # If file is not .bam (so fastq)
  else
    # 1. Make output directory
    mkdir fastq_output

    # 2. Copy the fastq file(s) to the output directory
    for path in "${sampleFiles[@]}"; do
      strippedPath=$(echo "$path" | tr -d '\\r')
      cp $strippedPath fastq_output
    done
  fi
  '''
}

process fastqcQualityControl {
  containerOptions "--bind ${params.bindFolder}"
  publishDir "${params.outDir}/fastqc/", mode: 'copy'

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  val fastqDir

  output:
  file "*_fastqc.zip"
  

  shell:
  '''
  for file in !{fastqDir}/*; do
    fastqc ${file} -o . --noextract
  done
  '''
}

process alignWithSTAR {
  containerOptions "--bind ${params.bindFolder}"
  publishDir "${params.outDir}/star/", mode: 'copy', pattern: "*/*.{gz}"

  time '6h'
  memory '50 GB'
  cpus 4

  input:
  path sample_dir

  output:
  path "*/*_Aligned.out.bam", emit: bamFile
  path "*/*.gz"

  shell:
  '''
  # Get path of the first file in the input directory
  firstFile=$(ls -1 "!{sample_dir}" | sort | head -n 1)

  # Extract sample name from the file path
  IFS='/' read -r -a directories <<< "${firstFile}"
  fileName="${directories[-1]}"

  # Remove .gz extension from the file name
  sampleName=${fileName%".gz"}

  # Remove the fastq extension from the file name
  extension=$(awk -F '.' '{print $NF}' <<< "$sampleName")
  sampleName=${sampleName%".$extension"}

  # For paired end, remove _1 and _2 suffix from the sample name
  if [[ $sampleName == *_1 || $sampleName == *_2 ]]; then
    sampleName="${sampleName%??}"
  fi

  # Create output directory
  mkdir ${sampleName}

  # Determine allowed number of mismatches based on read length
  readLength=$(samtools view !{sample_dir}/${firstFile} |head -n1 |awk '{print $10}'|tr -d "\\n" |wc -m)

  if [ $readLength -ge 90 ]; then
    numMism=4
  elif [ $readLength -ge 60 ]; then
    numMism=3
  else
    numMism=2
  fi

  # Set different --readFilesIn values for paired and single end
  if [[ -f !{sample_dir}/${sampleName}_2.fastq ]]; then
     let numMism=$numMism*2
     readFilesInArgument="--readFilesIn !{sample_dir}/${sampleName}_1.fastq.gz !{sample_dir}/${sampleName}_2.fastq.gz"
  else
     readFilesInArgument="--readFilesIn !{sample_dir}/${sampleName}.fastq.gz"
  fi

  # Run the STAR command
  STAR --runThreadN 8 \
  --outFileNamePrefix ${sampleName}/${sampleName}_ \
  --outSAMtype BAM Unsorted \
  --genomeDir !{params.refDir} \
  --genomeLoad NoSharedMemory \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNmax ${numMism} \
  --twopassMode Basic \
  --quantMode GeneCounts \
  --outSAMunmapped Within \
  --readFilesCommand zcat \
  ${readFilesInArgument}

  # Gzip all output files
  gzip ${sampleName}/*.tab
  gzip ${sampleName}/*.out
'''
}

process sortBAM {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path sample

  output:
  path "${sample.baseName}.sorted.bam"
  
  script:
  """
  samtools sort $sample -o ${sample.baseName}.sorted.bam
  """
}

process markDuplicates {
  containerOptions "--bind ${params.bindFolder}"
  publishDir "${params.outDir}/mark_duplicates/", mode: 'copy', pattern: "*.{gz}"

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path bam_file

  output:
  path "${bam_file.getBaseName()}.dupes.bam", emit: bamFile
  path "${bam_file.getBaseName()}_dupes.txt.gz"

  script:
  """
  java -Xmx10g -jar /usr/bin/picard.jar MarkDuplicates \
      I=${bam_file} \
      O=${bam_file.getBaseName()}.dupes.bam \
      M=${bam_file.getBaseName()}_dupes.txt

  gzip ${bam_file.getBaseName()}_dupes.txt
  """
}

process QCwithRNASeqMetrics {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '12 GB'
  cpus 1

  publishDir "${params.outDir}/rna_seq_metrics", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}_rnaseqmetrics"
  
  script:
  """
  java -Xmx10g -jar /usr/bin/picard.jar CollectRnaSeqMetrics \
  I=${sample} \
  O=${sample.baseName}_rnaseqmetrics \
  REF_FLAT=${params.refFlat} \
  STRAND=NONE \
  RIBOSOMAL_INTERVALS=${params.ribosomalIntervalList}
  """
}

process QCwithMultipleMetrics {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '12 GB'
  cpus 1

  publishDir "${params.outDir}/multiple_metrics", mode: 'copy'

  input:
  path sample

  output:
  path "${sample.baseName}/*"
  
  script:
  """
  mkdir ${sample.baseName}
  
  java -Xmx10g -jar /usr/bin/picard.jar CollectMultipleMetrics I=${sample} \
  O=${sample.baseName}/multiple_metrics \
  R=${params.referenceGenome} \
  PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=QualityScoreDistribution \
  PROGRAM=MeanQualityByCycle \
  PROGRAM=CollectInsertSizeMetrics \
  """
}

process identifyAlternativeSplicingSitesrMATS {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "${params.outDir}/rmats", mode: 'copy'

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
  python /usr/bin/rmats_turbo_v4_1_2/rmats.py --b1 config.txt \
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
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "${params.outDir}/leafcutter", mode: 'copy'

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
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '10 GB'
  cpus 1

  publishDir "${params.outDir}/cram", mode: 'copy'

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
    // Load sample list from text file
    String samplePaths = new File(params.sampleFile).text
    String[] samplePathsList = samplePaths.split('\n')
    samples = Channel.from(samplePathsList)

    // Run pipeline
    convertBAMToFASTQ(samples)
    fastqcQualityControl(convertBAMToFASTQ.out)
    alignWithSTAR(convertBAMToFASTQ.out)
    sortBAM(alignWithSTAR.out.bamFile)
    markDuplicates(sortBAM.out)
    QCwithRNASeqMetrics(markDuplicates.out.bamFile)
    QCwithMultipleMetrics(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out.bamFile)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out.bamFile)
    convertBAMToCRAM(markDuplicates.out.bamFile)
}