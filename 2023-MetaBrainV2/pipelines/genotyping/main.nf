nextflow.enable.dsl=2

process indexBam {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  
  output:
  path "indexBam/${bam_file.SimpleName}.bai"

  input:
  path bam_file

  script:
  """
  mkdir indexBam
  gatk --java-options "-Xmx64g" BuildBamIndex \
  I=${bam_file} \
  O=indexBam/${bam_file.SimpleName}.bai
  """
}

process encodeConvert {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "fixed/${sample_bam.SimpleName}.bam"
  
  script:
  """
  mkdir fastqc unzipped fixed
  fastqc ${sample_bam} --outdir fastqc
  unzip fastqc/*.zip -d unzipped
  num=\$(sed '6!d' unzipped/${sample_bam.SimpleName}_fastqc/fastqc_data.txt | cut -d " " -f 4)
  if [ echo "\$num < 1.8" | bc ]  
  then
      gatk --java-options "-Xmx16g" FixMisencodedBaseQualityReads \
      -I ${sample_bam} \
      -O fixed/${sample_bam.SimpleName}.bam
  else
      mv ${sample_bam} fixed
  fi
  """
}

process splitNCigarReads {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "split/${sample_bam.SimpleName}-splitreads.bam"
  
  script:
  """
  mkdir split
  gatk --java-options "-Xmx16g" SplitNCigarReads \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O split/${sample_bam.SimpleName}-splitreads.bam \
  """
}

process baseRecalibrator {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path split_bam
  path vcf_index
  path bam_index
  
  output:
  path "recal/${split_bam.SimpleName}.table"
  
  script:
  """
  mkdir recal
  gatk --java-options "-Xmx16g" BaseRecalibrator \
  -I ${split_bam} \
  -R ${params.referenceGenome} \
  --known-sites ${params.knownSites} \
  -O recal/${split_bam.SimpleName}.table
  """
}

process applyBQSR {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  path table_file

  output:
  path "bqsr/${sample_bam.SimpleName}.bqsr.bam"
  
  script:
  """
  mkdir bqsr
  gatk --java-options "-Xmx16g" ApplyBQSR \
  -I ${sample_bam} \
  -R ${params.referenceGenome} --bqsr-recal-file ${table_file} \
  -O bqsr/${sample_bam.SimpleName}.bqsr.bam 
  """
}

process haplotypeCaller {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 4

  input:
  path sample_bam
  path bam_index
  
  output:
  path "gvcf/${sample_bam.SimpleName}.gvcf.gz"
  
  script:
  """
  mkdir gvcf
  gatk --java-options "-Xmx16g" HaplotypeCaller \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O gvcf/${sample_bam.SimpleName}.gvcf.gz \
  -ERC GVCF
  """
}

process indexGvcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path gvcf_file
  
  output:
  path "index/${gvcf_file.SimpleName}.gvcf.gz.tbi"

  script:
  """
  mkdir index
  gatk --java-options "-Xmx16g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O index/${gvcf_file.SimpleName}.gvcf.gz.tbi
  """
}


process indexJointGvcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path gvcf_file
  
  output:
  path "index/${gvcf_file.SimpleName}.gvcf.gz.tbi"

  script:
  """
  mkdir index
  gatk --java-options "-Xmx16g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O index/${gvcf_file.SimpleName}.gvcf.gz.tbi
  """
}

process combineGvcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path gvcf_samples
  path vcf_indexes
  
  output:
  path "combined/combined.gvcf.gz"
  
  script:
  """
  mkdir combined
  sample_string=\$(echo ${gvcf_samples} | sed "s/ / --variant /g")
  index_string=\$(echo ${vcf_indexes} | sed "s/ / --read-index /g")
  gatk --java-options "-Xmx16g" CombineGVCFs \
  -R ${params.referenceGenome}\
  --variant \$sample_string \
  --read-index \$index_string \
  -O combined/combined.gvcf.gz
  """
}

process jointGenotype {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_gvcf
  path gvcf_index
  
  output:
  path "${sample_gvcf.SimpleName}.gvcf.gz"
  
  script:
  """
  gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R ${params.referenceGenome}\
  -V ${sample_gvcf} \
  -O ${sample_gvcf.SimpleName}.gvcf.gz
  """
}


workflow {
    bam_files = Channel.fromPath("${params.bamDir}/*.bam")
    indexBam(bam_files)
    encodeConvert(bam_files)
    splitNCigarReads(encodeConvert.output) 
    baseRecalibrator(splitNCigarReads.output, params.knownSitesIndex, indexBam.output)
    applyBQSR(splitNCigarReads.output, baseRecalibrator.output)
    gvcf_files = haplotypeCaller(applyBQSR.output, indexBam.output).collect()
    gvcf_indexes = indexGvcf(haplotypeCaller.output).collect()
    combineGvcf(gvcf_files, gvcf_indexes)
    indexJointGvcf(combineGvcf.output)
    jointGenotype(combineGvcf.output, indexJointGvcf.output)
}