nextflow.enable.dsl=2

process cramToBam {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '4 GB'
  cpus 1
  maxRetries 16
  
  output:
  path "${cram_file.SimpleName}.bam"

  input:
  path cram_file

  script:
  """
  samtools view -b  -T ${params.referenceGenome} -o ${cram_file.SimpleName}.bam ${cram_file}
  """
}

process indexBam {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '6 GB'
  cpus 1
  maxRetries 16
  
  output:
  path "indexBam/${bam_file.SimpleName}.bai"

  input:
  path bam_file

  script:
  """
  mkdir indexBam
  gatk --java-options "-Xmx4g" BuildBamIndex \
  I=${bam_file} \
  O=indexBam/${bam_file.SimpleName}.bai
  """
}

process encodeConvert {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '6 GB'
  cpus 1
  maxRetries 16
  
  input:
  path sample_bam
  path index_file

  output:
  path "fixed/${sample_bam.SimpleName}.bam", emit: bam_file
  path index_file, emit: index_file
  
  script:
  """
  mkdir fastqc unzipped fixed
  fastqc ${sample_bam} --outdir fastqc
  unzip fastqc/*.zip -d unzipped
  num=\$(sed '6!d' unzipped/${sample_bam.SimpleName}_fastqc/fastqc_data.txt | cut -d " " -f 4)
  if [ echo "\$num < 1.8" | bc ]  
  then
      gatk --java-options "-Xmx4g" FixMisencodedBaseQualityReads \
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
  memory '8 GB'
  cpus 1
  maxRetries 16

  input:
  path sample_bam
  path index_file

  output:
  path "split/${sample_bam.SimpleName}-splitreads.bam", emit: bam_file
  path index_file, emit: index_file
  
  script:
  """
  mkdir split
  gatk --java-options "-Xmx6g" SplitNCigarReads \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O split/${sample_bam.SimpleName}-splitreads.bam \
  """
}

process AddOrReplaceReadGroups {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '6 GB'
  cpus 1
  maxRetries 16

  input:
  path split_bam
  path index_file
  
  output:
  path "readgroup/${split_bam.SimpleName}.bam", emit: bam_file
  path index_file, emit: index_file
  
  script:
  """
  mkdir readgroup
  gatk --java-options "-Xmx4g" AddOrReplaceReadGroups \
  I=${split_bam} \
  O=readgroup/${split_bam.SimpleName}.bam \
  RGID=4 \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=${split_bam.SimpleName}
  """
}


process baseRecalibrator {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '16 GB'
  cpus 1
  maxRetries 16

  input:
  path split_bam
  path index_file
  path vcf_index
  
  output:
  path "recal/${split_bam.SimpleName}.table", emit: table_file
  path index_file, emit: index_file
  
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
  memory '6 GB'
  cpus 1
  maxRetries 16

  input:
  path sample_bam
  path index_file
  path table_file

  output:
  path "bqsr/${sample_bam.SimpleName}.bqsr.bam", emit: bam_file
  path index_file, emit: index_file
  
  script:
  """
  mkdir bqsr
  gatk --java-options "-Xmx4g" ApplyBQSR \
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
  gatk --java-options "-Xmx12g" HaplotypeCaller \
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
  memory '8 GB'
  cpus 1
  maxRetries 16

  input:
  path gvcf_file
  
  output:
  path "index/${gvcf_file.SimpleName}.gvcf.gz.tbi"

  script:
  """
  mkdir index
  gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O index/${gvcf_file.SimpleName}.gvcf.gz.tbi
  """
}


process indexJointGvcf {
  containerOptions '--bind /groups/'
  publishDir "/groups/umcg-biogen/tmp01/umcg-ogkourlias/combined_vcf", mode: 'move'
  errorStrategy 'retry'
  time '6h'
  memory '8 GB'
  cpus 1
  maxRetries 12

  input:
  path gvcf_file
  
  output:
  path "${gvcf_file.SimpleName}.vcf.gz.tbi", emit: tbi_file
  path "${gvcf_file}", emit: gvcf_file

  script:
  """
  gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O ${gvcf_file.SimpleName}.vcf.gz.tbi
  """
}

process combineGvcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '72h'
  memory '12 GB'
  cpus 1
  maxRetries 12

  input:
  path gvcf_files
  path gvcf_indexes
  
  output:
  path "combined/combined.gvcf.gz"
  
  script:
  """
  mkdir combined
  for GVCF in ${gvcf_files}
  do
    echo \$GVCF >> gvcf.list
  done
  gatk --java-options "-Xmx10g" CombineGVCFs \
  -R ${params.referenceGenome} \
  -V gvcf.list \
  -O combined/combined.gvcf.gz \
  --call-genotypes true
  """
}

process jointGenotype {
  publishDir "${params.outDir}", mode: 'copy'
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '18 GB'
  cpus 1
  maxRetries 18

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
  -O ${sample_gvcf.SimpleName}.gvcf.gz \
  --read-index ${gvcf_index}
  """
}

process chrSplit {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '8 GB'
  cpus 1
  maxRetries 18

  input:
  path gvcf_file
  path gvcf_index
  val i
  
  output:
  path "${gvcf_file.SimpleName}/chr${i}.vcf.gz"
  
  script:
  """
  mkdir ${gvcf_file.SimpleName}
  bcftools view ${gvcf_file} --regions ${i} -o ${gvcf_file.SimpleName}/chr${i}.vcf.gz -Oz
  """
}

process combineChrGvcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  publishDir "${chr_dir}", mode: 'copy'
  time '6h'
  memory '12 GB'
  cpus 1
  maxRetries 16

  input:
  path chr_dir
  
  output:
  path "combined/combined.vcf.gz"
  
  script:
  """
  mkdir combined
  for VCF in ${chr_dir}/*.vcf.gz
  do
    echo \$VCF >> vcf.list
  done
  gatk --java-options "-Xmx10g" CombineGVCFs \
  -R ${params.referenceGenome} \
  -V vcf.list \
  -O combined/combined_${chr_dir}.vcf.gz
  """
}


workflow {
    cram_files = Channel.fromPath("${params.cram_dir}/*/*.cram")
    cramToBam(cram_files)
    indexBam(cramToBam.output)
    encodeConvert(cramToBam.output, indexBam.output)
    splitNCigarReads(encodeConvert.output.bam_file, encodeConvert.output.index_file) 
    AddOrReplaceReadGroups(splitNCigarReads.output.bam_file, splitNCigarReads.output.index_file)
    baseRecalibrator(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, params.knownSitesIndex)
    applyBQSR(AddOrReplaceReadGroups.output.bam_file, AddOrReplaceReadGroups.output.index_file, baseRecalibrator.output.table_file)
    gvcf_files = haplotypeCaller(applyBQSR.output.bam_file, applyBQSR.output.index_file).collect()
    gvcf_indexes = indexGvcf(haplotypeCaller.output).collect()
    chrs = Channel.from( 1..22 )
    chrSplit(gvcf_files, gvcf_indexes, chrs)
    indexJointGvcf(chrSplit.output)
    jointGenotype(indexJointGvcf.output.gvcf_file, indexJointGvcf.output.tbi_file)
}