nextflow.enable.dsl=2

process convertBAMToFASTQ {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1
  
  input:
  val sample
  
  output:
  path "fastq_output", emit: fastqPath
  val "${sample[1]}", emit: sampleName
  val task.workDir, emit: workDir
  
  shell:
  '''
  # Split the input string by ';' to separate file paths
  IFS=';' read -r -a sampleFiles <<< "!{sample[0]}"

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
      arguments="fastq -0 fastq_output/!{sample[1]}.fastq.gz -n sorted_${sampleName}.bam"
    elif [ "${numFASTQFiles}" -gt 1 ];
    then
      arguments="fastq -1 fastq_output/!{sample[1]}_1.fastq.gz -2 fastq_output/!{sample[1]}_2.fastq.gz -n sorted_${sampleName}.bam"
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
  val sampleName

  output:
  file "${sampleName}/*_fastqc.zip"
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
  publishDir "${params.outDir}/star/", mode: 'copy', pattern: "*/*.{gz}"

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

  # Set different --readFilesIn values for paired and single end
  if [[ -f !{sampleDir}/${sampleName}_2.fastq.gz ]]; then
     let numMism=$numMism*2
     readFilesInArgument="--readFilesIn !{sampleDir}/${sampleName}_1.fastq.gz !{sampleDir}/${sampleName}_2.fastq.gz"
  else
     readFilesInArgument="--readFilesIn !{sampleDir}/${sampleName}.fastq.gz"
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

  time '6h'
  memory '12 GB'
  cpus 1

  publishDir "${params.outDir}/rna_seq_metrics", mode: 'copy'

  input:
  path samplePath
  val sampleName

  output:
  path "${sampleName}/${sampleName}_rnaseqmetrics"
  path "${sampleName}/${sampleName}.chart.pdf"
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
  """
}

process QCwithMultipleMetrics {
  containerOptions "--bind ${params.bindFolder}"

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

  time '6h'
  memory '8 GB'
  cpus 1

  publishDir "${params.outDir}/leafcutter", mode: 'copy'

  input:
  path samplePath
  val sampleName
  
  output:
  path "${sampleName}/*.junc"
  val task.workDir, emit: workDir
  
  shell:
  '''
  # 1. Index BAM file
  samtools index !{samplePath}

  # 2. Run regtools command
  mkdir !{sampleName}
  regtools junctions extract -s XS -a 8 -m 50 -M 500000 !{samplePath} -o !{sampleName}/!{sampleName}.junc 
  '''
}

process convertBAMToCRAM {
  containerOptions "--bind ${params.bindFolder}"

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

def extractSampleName(String path) {
    // Split the path by ';' to handle multiple paths
    def paths = path.split(';')

    // Initialize a set to store unique sample names
    Set<String> sampleNameSet = new HashSet<>()
    
    // Extract and process the sample names for each path
    paths.each { pathPart ->
        // Extract the sample name and remove _1 or _2 if present
        def sampleNamePart = pathPart.tokenize("/").last().replaceAll(/\.(fastq|fq|bam|gz)/, '').replaceAll(/_(1|2)$/, '')
        
        // Add the cleaned sample name to the set
        sampleNameSet << sampleNamePart
    }
    
    // Join the unique sample names with '_'
    def sampleName = sampleNameSet.join("_")
    
    return sampleName
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
    def expectedNumberOfFiles = [1, 5, 8, 2, 36, 1, 1];

    // Loop through expected folders and number of expected files
    for (int i = 0; i < expectedNumberOfFiles.size; i++) {
        
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
    String[] samplePathsArray = samplePaths.split('\n')

    // Create a list of sample names extracted from the sample paths
    List<String> sampleNamesList = samplePathsArray.collect { extractSampleName(it) }
    String[] sampleNamesArray = sampleNamesList as String[]
    
    // Create channels of the sample paths and sample names
    samplePathsChannel = Channel.of(samplePathsArray)
    sampleNamesChannel = Channel.of(sampleNamesArray)

    // Merge the channels
    mergedChannel = samplePathsChannel.merge(sampleNamesChannel)

    // Remove samples from the channel that are already in the output folder
    filteredChannel = mergedChannel.filter { !checkIfSampleIsProcessed(params.outDir, it[1]) }

    // Run pipeline
    convertBAMToFASTQ(mergedChannel)
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