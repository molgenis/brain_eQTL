nextflow.enable.dsl=2

process convertBAMToFASTQ {
  containerOptions "--bind ${params.bindFolder}"

  afterScript """
  mkdir -p ${params.outDir}/bam_to_fastq/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/bam_to_fastq/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/bam_to_fastq/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/bam_to_fastq/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/bam_to_fastq/${sampleName}/command.log.gz
  """
  
  errorStrategy 'retry'
  maxRetries 999
  
  input:
  tuple val(sampleName), val(samplePath)

  time '6h'
  memory '8 GB'
  cpus 1
  
  output:
  path "fastq_output", emit: fastqPath
  val sampleName, emit: sampleName
  val task.workDir, emit: workDir
  
  shell:
  '''
  # Split the input string by ';' to separate file paths
  IFS=';' read -r -a sampleFiles <<< "!{samplePath}"

  # Get the extension of the sample
  extension=$(awk -F '.' '{print $NF}' <<< "${sampleFiles[0]}")

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
    samtools sort -n ${sample} -o "sorted_!{sampleName}.bam"

    # 3. Create command arguments
    if [[ "${numFASTQFiles}" -eq 1 ]]; 
    then
      arguments="fastq -0 fastq_output/!{sampleName}.fastq.gz -n sorted_!{sampleName}.bam"
    elif [ "${numFASTQFiles}" -gt 1 ];
    then
      arguments="fastq -1 fastq_output/!{sampleName}_1.fastq.gz -2 fastq_output/!{sampleName}_2.fastq.gz -n sorted_!{sampleName}.bam"
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

    # 2. Create symlink from the fastq file(s) to the output directory
    for path in "${sampleFiles[@]}"; do
      strippedPath=$(echo "$path" | tr -d '\\r')
      IFS='/' read -r -a directories <<< "${path}"
      fileName="${directories[-1]}"
      finalOut=fastq_output/"$fileName"
      ln -s $strippedPath $finalOut
    done
  fi
  '''
}

process fastqcQualityControl {
  containerOptions "--bind ${params.bindFolder}"
  publishDir "${params.outDir}/fastqc/", mode: 'copy'
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/fastqc/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/fastqc/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/fastqc/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/fastqc/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/fastqc/${sampleName}/command.log.gz
  """

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  val fastqDir
  val sampleName

  output:
  file "${sampleName}/*_fastqc.zip"
  path ".command.*"
  val task.workDir, emit: workDir

  shell:
  '''
  mkdir !{sampleName}

  for file in !{fastqDir}/*; do
    fastqc ${file} -o !{sampleName} --noextract
  done
  '''
}

process alignWithSTAR {
  containerOptions "--bind ${params.bindFolder}"
  publishDir "${params.outDir}/star/", mode: 'copy', pattern: "${sampleName}/*.{gz}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/star/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/star/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/star/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/star/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/star/${sampleName}/command.log.gz
  """

  time '6h'
  memory '50 GB'
  cpus 4

  input:
  path sampleDir
  val sampleName

  output:
  path "*/*_Aligned.out.bam", emit: bamFile
  path "*/*.gz"
  val sampleName, emit: sampleName
  val task.workDir, emit: workDir

  shell:
  '''
  # Get path of the first file in the input directory
  firstFile=$(ls -1 "!{sampleDir}" | sort | head -n 1)

  # Extract sample name from the file path
  IFS='/' read -r -a directories <<< "${firstFile}"
  fileName="${directories[-1]}"

  # Create output directory
  mkdir !{sampleName}

  # Determine allowed number of mismatches based on read length
  readLength=$(samtools view !{sampleDir}/${firstFile} |head -n1 |awk '{print $10}'|tr -d "\\n" |wc -m)

  if [ $readLength -ge 90 ]; then
    numMism=4
  elif [ $readLength -ge 60 ]; then
    numMism=3
  else
    numMism=2
  fi

  # Count fastq files
  fastqCount=$(ls -1 "!{sampleDir}" | wc -l)

  # If data is paired end (more than one fastq file), double nuMism parameter
  if [[ "$fastqCount" -gt 1 ]]; then
     let numMism=$numMism*2
  fi

  # Set readFilesIn argument
  readFilesInArgument="--readFilesIn"
  for file in !{sampleDir}/*; do
    readFilesInArgument+=" $file "
  done

  # Run the STAR command
  STAR --runThreadN 8 \
  --outFileNamePrefix !{sampleName}/!{sampleName}_ \
  --outSAMtype BAM Unsorted \
  --genomeDir !{params.refDir} \
  --genomeLoad NoSharedMemory \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNmax ${numMism} \
  --twopassMode Basic \
  --quantMode GeneCounts \
  --readFilesCommand zcat \
  --outSAMunmapped Within \
  ${readFilesInArgument}

  # Gzip all output files
  gzip !{sampleName}/*.tab
  gzip !{sampleName}/*.out
'''
}

process sortBAM {
  containerOptions "--bind ${params.bindFolder}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/sort_bam/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/sort_bam/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/sort_bam/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/sort_bam/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/sort_bam/${sampleName}/command.log.gz
  """

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path samplePath
  val sampleName

  output:
  path "${sampleName}.sorted.bam", emit: bamFile
  val sampleName, emit: sampleName
  val task.workDir, emit: workDir
  
  script:
  """
  samtools sort ${samplePath} -o ${sampleName}.sorted.bam
  """
}

process markDuplicates {
  containerOptions "--bind ${params.bindFolder}"
  publishDir "${params.outDir}/mark_duplicates/", mode: 'copy', pattern: "*/*.{gz}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/mark_duplicates/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/mark_duplicates/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/mark_duplicates/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/mark_duplicates/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/mark_duplicates/${sampleName}/command.log.gz
  """

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path bam_file
  val sampleName

  output:
  path "${sampleName}/${sampleName}.duplicates.bam", emit: bamFile
  path "${sampleName}/${sampleName}_duplicates.txt.gz"
  val sampleName, emit: sampleName
  val task.workDir, emit: workDir

  script:
  """
  mkdir ${sampleName}

  java -Xmx10g -jar /usr/bin/picard.jar MarkDuplicates \
      I=${bam_file} \
      O=${sampleName}/${sampleName}.duplicates.bam \
      M=${sampleName}/${sampleName}_duplicates.txt

  gzip ${sampleName}/${sampleName}_duplicates.txt
  """
}

process QCwithRNASeqMetrics {
  containerOptions "--bind ${params.bindFolder}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/rna_seq_metrics/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/rna_seq_metrics/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/rna_seq_metrics/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/rna_seq_metrics/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/rna_seq_metrics/${sampleName}/command.log.gz
  """

  time '6h'
  memory '12 GB'
  cpus 1

  publishDir "${params.outDir}/rna_seq_metrics", mode: 'copy'

  input:
  path samplePath
  val sampleName

  output:
  path "${sampleName}/${sampleName}_rnaseqmetrics.gz"
  path "${sampleName}/${sampleName}.chart.pdf.gz"
  val sampleName
  val task.workDir, emit: workDir
  
  script:
  """
  mkdir ${sampleName}

  java -Xmx10g -jar /usr/bin/picard.jar CollectRnaSeqMetrics \
  I=${samplePath} \
  O=${sampleName}/${sampleName}_rnaseqmetrics \
  CHART_OUTPUT=${sampleName}/${sampleName}.chart.pdf \
  REF_FLAT=${params.refFlat} \
  STRAND=NONE \
  RIBOSOMAL_INTERVALS=${params.ribosomalIntervalList}

  gzip ${sampleName}/*
  """
}

process QCwithMultipleMetrics {
  containerOptions "--bind ${params.bindFolder}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/multiple_metrics/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/multiple_metrics/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/multiple_metrics/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/multiple_metrics/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/multiple_metrics/${sampleName}/command.log.gz
  """

  time '6h'
  memory '12 GB'
  cpus 1

  publishDir "${params.outDir}/multiple_metrics", mode: 'copy'

  input:
  path samplePath
  val sampleName

  output:
  path "${sampleName}/*"
  val task.workDir, emit: workDir
  
  script:
  """
  mkdir ${sampleName}
  
  java -Xmx10g -jar /usr/bin/picard.jar CollectMultipleMetrics I=${samplePath} \
  O=${sampleName}/multiple_metrics \
  R=${params.referenceGenome} \
  PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=QualityScoreDistribution \
  PROGRAM=MeanQualityByCycle \
  PROGRAM=CollectInsertSizeMetrics 

  gzip ${sampleName}/*
  """
}

process identifyAlternativeSplicingSitesrMATS {
  containerOptions "--bind ${params.bindFolder}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/rmats/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/rmats/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/rmats/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/rmats/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/rmats/${sampleName}/command.log.gz
  """

  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "${params.outDir}/rmats", mode: 'copy'

  input:
  path samplePath
  val sampleName
  
  output:
  path "${sampleName}/*.txt.gz"
  val task.workDir, emit: workDir
  
  shell:
  '''
  # 1. Check if the BAM file is derived from single or paired end reads
  numFASTQFiles=$(samtools view -H !{samplePath} | \
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
  readLength=$(samtools view !{samplePath} |head -n1 |awk '{print $10}'|tr -d "\\n" |wc -m)

  # 3. Create config file
  echo !{samplePath} > config.txt

  # 4. Run rMATS command
  python /usr/bin/rmats_turbo_v4_1_2/rmats.py --b1 config.txt \
  --gtf !{params.gtfAnnotationFile} \
  --readLength ${readLength} \
  --od !{sampleName} \
  --tmp rmats_tmp \
  --task both \
  -t ${end}  \
  --statoff

  # 5. Gzip all output files
  gzip !{sampleName}/*.txt
  '''
}

process identifyAlternativeSplicingSitesLeafCutter {
  containerOptions "--bind ${params.bindFolder}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/leafcutter/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/leafcutter/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/leafcutter/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/leafcutter/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/leafcutter/${sampleName}/command.log.gz
  """

  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "${params.outDir}/leafcutter", mode: 'copy'

  input:
  path samplePath
  val sampleName
  
  output:
  path "${sampleName}/*.junc.gz"
  val task.workDir, emit: workDir
  
  shell:
  '''
  # 1. Index BAM file
  samtools index !{samplePath}

  # 2. Run regtools command
  mkdir !{sampleName}
  regtools junctions extract -s XS -a 8 -m 50 -M 500000 !{samplePath} -o !{sampleName}/!{sampleName}.junc 

  # 3. Gzip the resulting junctions file
  gzip !{sampleName}/!{sampleName}.junc
  '''
}

process convertBAMToCRAM {
  containerOptions "--bind ${params.bindFolder}"
  errorStrategy 'retry'
  maxRetries 999

  afterScript """
  mkdir -p ${params.outDir}/cram/${sampleName}
  gzip ./.command.*
  cp ./.command.sh.gz ${params.outDir}/cram/${sampleName}/command.sh.gz
  cp ./.command.err.gz ${params.outDir}/cram/${sampleName}/command.err.gz
  cp ./.command.out.gz ${params.outDir}/cram/${sampleName}/command.out.gz
  cp ./.command.log.gz ${params.outDir}/cram/${sampleName}/command.log.gz
  """

  time '6h'
  memory '10 GB'
  cpus 1

  publishDir "${params.outDir}/cram", mode: 'copy'

  input:
  path samplePath
  val sampleName
  
  output:
  path "${sampleName}/${sampleName}.cram"
  val task.workDir, emit: workDir
  
  script:
  """
  mkdir ${sampleName}
  samtools view -T ${params.referenceGenome} -C -o ${sampleName}/${sampleName}.cram ${samplePath}
  """
}

def splitSampleNamesandPaths(String input) {
    return new Tuple2(input.split(',')[0], input.split(',')[1])
}

def checkIfSampleIsProcessed(String folderName, String sampleName) {
    
    // An array containing all the folders that are expected for a succesful pipeline run for a sample
    def expectedFolders = [
        folderName + '/fastqc/' + sampleName,
        folderName + '/star/' + sampleName,
        folderName + '/multiple_metrics/' + sampleName,
        folderName + '/rna_seq_metrics/' + sampleName,
        folderName + '/rmats/' + sampleName,
        folderName + '/mark_duplicates/' + sampleName,
        folderName + '/cram/' + sampleName,
    ];

    // An array containing the expected number of files in order of the folders above
    def expectedNumberOfFiles = [5, 9, 12, 6, 40, 5, 5];

    // Loop through expected folders and number of expected files
    for (int i = 0; i < expectedNumberOfFiles.size; i++) {
//        println "Looking for folder: "+expectedFolders[i]    
        // Return false if the expected folder does not exist
        if (!new File(expectedFolders[i]).exists()) {
           return false;
        }

        // Return flase if the number of items in folder does not match expected number of items
        if (new File(expectedFolders[i]).list().length < expectedNumberOfFiles[i]){
           return false;
        }
      }

      return true;
}

process removeWorkDirs {
  containerOptions "--bind ${params.bindFolder}"
  errorStrategy 'retry'
  maxRetries 5

  time '1h'
  memory '1 GB'
  cpus 1

  input:
  val bamToFastqWorkDir
  val fastQCWorkDir
  val alignWithStarWorkDir
  val sortBamWorkDir
  val markDuplicatesWorkDir
  val QCwithRNASeqMetricsWorkDir
  val QCwithMultipleMetricsWorkDir
  val identifyAlternativeSplicingSitesrMATSWorkDir
  val identifyAlternativeSplicingSitesLeafCutterWorkDir
  val convertBAMToCRAMWorkDir
  
  script:
  """
  sleep 5

  rm -r ${bamToFastqWorkDir} || echo 'Failed to remove work directory'
  rm -r ${fastQCWorkDir} || echo 'Failed to remove work directory'
  rm -r ${alignWithStarWorkDir} || echo 'Failed to remove work directory'
  rm -r ${sortBamWorkDir} || echo 'Failed to remove work directory'
  rm -r ${markDuplicatesWorkDir} || echo 'Failed to remove work directory'
  rm -r ${QCwithRNASeqMetricsWorkDir} || echo 'Failed to remove work directory'
  rm -r ${QCwithMultipleMetricsWorkDir} || echo 'Failed to remove work directory'
  rm -r ${identifyAlternativeSplicingSitesrMATSWorkDir} || echo 'Failed to remove work directory'
  rm -r ${identifyAlternativeSplicingSitesLeafCutterWorkDir} || echo 'Failed to remove work directory'
  rm -r ${convertBAMToCRAMWorkDir} || echo 'Failed to remove work directory'
  """
}

workflow {
    // Load list with sample paths from the input text file
    String samplePaths = new File(params.sampleFile).text
    String[] splittedSamplePaths = samplePaths.split('\n')
    println "Nr sample paths: "+splittedSamplePaths.length
    
    List<Tuple2<String, String>> tupleList = splittedSamplePaths.collect { splitSampleNamesandPaths(it) }
    println "Nr sample tuples: "+tupleList.size()

    Tuple2<String, String>[] tupleArray = tupleList as Tuple2<String, String>[]
    println "Nr sample tuples in arr: "+tupleArray.length

    def channel = Channel.of(tupleArray)

    // Remove samples from the channels that are already in the output folder
    filteredChannel = channel.filter { !checkIfSampleIsProcessed(params.outDir, it[0]) }

    // Run pipeline
    convertBAMToFASTQ(filteredChannel)
    fastqcQualityControl(convertBAMToFASTQ.out.fastqPath, convertBAMToFASTQ.out.sampleName)
    alignWithSTAR(convertBAMToFASTQ.out.fastqPath, convertBAMToFASTQ.out.sampleName)
    sortBAM(alignWithSTAR.out.bamFile, alignWithSTAR.out.sampleName)
    markDuplicates(sortBAM.out.bamFile, sortBAM.out.sampleName)
    QCwithRNASeqMetrics(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    QCwithMultipleMetrics(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    identifyAlternativeSplicingSitesrMATS(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    identifyAlternativeSplicingSitesLeafCutter(markDuplicates.out.bamFile, markDuplicates.out.sampleName)
    convertBAMToCRAM(markDuplicates.out.bamFile, markDuplicates.out.sampleName)

    // Remove all work directories for sample
    removeWorkDirs(
      convertBAMToFASTQ.out.workDir,
      fastqcQualityControl.out.workDir,
      alignWithSTAR.out.workDir,
      sortBAM.out.workDir,
      markDuplicates.out.workDir,
      QCwithRNASeqMetrics.out.workDir,
      QCwithMultipleMetrics.out.workDir,
      identifyAlternativeSplicingSitesrMATS.out.workDir,
      identifyAlternativeSplicingSitesLeafCutter.out.workDir,
      convertBAMToCRAM.out.workDir,
    )
}