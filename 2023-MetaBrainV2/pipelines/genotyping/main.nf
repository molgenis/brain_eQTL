nextflow.enable.dsl=2

params.referenceGenome = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa"
params.outDir = "/groups/umcg-biogen/tmp01/umcg-jbakker/results"
params.bamDir = "/groups/umcg-fg/tmp01/projects/downstreamer/2023-09-PublicRNASeqGenotypeCalls/metabrain/2023-MetaBrainV2/pipelines/genotyping/bam"
params.knownSites = "/groups/umcg-fg/tmp01/projects/downstreamer/2023-09-PublicRNASeqGenotypeCalls/variants/GCF_000001405.40.new.gz"
params.knownSitesIndex = "/groups/umcg-fg/tmp01/projects/downstreamer/2023-09-PublicRNASeqGenotypeCalls/variants/GCF_000001405.40.new.gz.tbi"

def getEncoding(bam_file) {
  """
  fastqc ${bam_file} --outdir fastqc
  unzip fastqc/*.zip -d unzipped
  unzipped/*/fastqc_data.txt
  sed '6!d' unzipped/*/fastqc_data.txt | cut -d " " -f 4
  """
}

process indexBam {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1
  
  output:
  path "indexBam/${bam_file.baseName}.bai"

  input:
  path bam_file

  script:
  """
  mkdir indexBam
  gatk --java-options "-Xmx64g" BuildBamIndex \
  I=${bam_file} \
  O=indexBam/${bam_file.baseName}.bai
  """
}

process indexVcf {
  containerOptions '--bind /groups/'
  time '6h'
  memory '64 GB'
  cpus 1
  
  output:
  path "indexVcf/*"

  input:
  path vcf_file

  script:
  """
  mkdir indexVcf
  gatk --java-options "-Xmx64g" IndexFeatureFile \
  --arguments_file ${vcf_file}
  """
}

process encodeConvert {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "fixed/${sample_bam.baseName}.bam"
  
  script:
  """
  mkdir fastqc unzipped fixed
  fastqc ${sample_bam} --outdir fastqc
  unzip fastqc/*.zip -d unzipped
  num=\$(sed '6!d' unzipped/${sample_bam.baseName}_fastqc/fastqc_data.txt | cut -d " " -f 4)
  if [ echo "\$num < 1.8" | bc ]  
  then
      gatk --java-options "-Xmx16g" FixMisencodedBaseQualityReads \
      -I ${sample_bam} \
      -O fixed/${sample_bam.baseName}.bam
  else
      mv ${sample_bam} fixed
  fi
  """
}

process createSeqDict {
  containerOptions '--bind /groups/'
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path ref_file
  
  output:
  path "${params.referenceGenome.baseName}.dict"
  
  script:
  """
  gatk --java-options "-Xmx8g" CreateSequenceDictionary \
  -R ${params.referenceGenome}\
  -O ${sample_bam.baseName.baseName}.dict
  """
}

process splitNCigarReads {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "split/${sample_bam.baseName}-splitreads.bam"
  
  script:
  """
  mkdir split
  gatk --java-options "-Xmx16g" SplitNCigarReads \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O split/${sample_bam.baseName}-splitreads.bam \
  """
}

process baseRecalibrator {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path split_bam
  path vcf_index
  path bam_index
  
  output:
  path "recal/${split_bam.baseName}.table"
  
  script:
  """
  mkdir recal
  gatk --java-options "-Xmx16g" BaseRecalibrator \
  -I ${split_bam} \
  -R ${params.referenceGenome} \
  --known-sites ${params.knownSites} \
  -O recal/${split_bam.baseName}.table
  """
}

process printReads {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  path vcf_index

  output:
  path "preads/${sample_bam.baseName}.bam"
  
  script:
  """
  mkdir preads
  gatk --java-options "-Xmx16g" PrintReads \
  -R ${params.referenceGenome}  \
  -I ${sample_bam} \
  -o preads/${sample_bam.baseName} \
  -BQSR ${table_file} \
  -nct 2
  """
}

process haplotypeCaller {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "gvcf/${sample_bam.baseName}.gvcf.gz"
  
  script:
  """
  mkdir gvcf
  gatk --java-options "-Xmx16g" HaplotypeCaller \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O gvcf/${sample_bam.baseName}.gvcf.gz \
  -ERC GVCF
  """
}

process indexGvcf {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path gvcf_file
  
  output:
  path "index/${gvcf_file.baseName}.gvcf.gz.tbi"
  path gvcf_file, emit: gvcf_file

  script:
  """
  mkdir index
  gatk --java-options "-Xmx16g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O index/${gvcf_file.baseName}.gvcf.gz.tbi
  """
}

process combineGvcf {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path gvcf_samples
  
  output:
  path "combined/combined.gvcf.gz"
  
  script:
  """
  mkdir combined
  gatk --java-options "-Xmx16g" CombineGVCFs \
  -R ${params.referenceGenome}\
  -V ${gvcf_samples} \
  -O combined/combined.gvcf.gz
  """
}

process jointGenotype {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_gvcf
  
  output:
  path "${sample_bam.baseName}.gvcf.gz"
  
  script:
  """
  gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O ${sample_bam.baseName}.gvcf.gz \
  -ERC GVCF
  """
}
;

workflow {
    bam_files = Channel.fromPath("${params.bamDir}/*.bam")
    //indexVcf(params.knownSites)
    indexBam(bam_files)
    encodeConvert(bam_files)
    splitNCigarReads(encodeConvert.output) 
    baseRecalibrator(splitNCigarReads.output, params.knownSitesIndex, indexBam.output)
    printReads(splitNCigarReads.output, baseRecalibrator.output)
    haplotypeCaller(printReads.output)
    indexGvcf(haplotypeCaller.output)
    combineGvcf(haplotypeCaller.output)
    //indexGvcf(gvcf_files)
}