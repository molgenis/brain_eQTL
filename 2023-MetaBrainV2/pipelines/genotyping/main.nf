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
  memory '12 GB'
  cpus 4

  input:
  path sample_bam
  path bam_index
  
  output:
  path "gvcf/${sample_bam.SimpleName}.gvcf.gz"
  
  script:
  """
  mkdir gvcf
  gatk --java-options "-Xmx10g" HaplotypeCaller \
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
  path "${gvcf_file.SimpleName}.gvcf.gz.tbi"
  
  script:
  """
  gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O ${gvcf_file.SimpleName}.gvcf.gz.tbi
  """
}


process indexJointGvcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '8 GB'
  cpus 1
  maxRetries 12

  input:
  path gvcf_file
  
  output:
  path "${gvcf_file.SimpleName}.gvcf.gz.tbi", emit: tbi_file
  path "${gvcf_file}", emit: gvcf_file

  script:
  """
  gatk --java-options "-Xmx6g" IndexFeatureFile \
  -I ${gvcf_file} \
  -O ${gvcf_file.SimpleName}.gvcf.gz.tbi
  """
}

process combineGvcf {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '12h'
  memory '12 GB'
  cpus 1
  maxRetries 12

  input:
  path gvcfs
  path gvcf_idx
  val i
  
  output:
  path "chr${i}.gvcf.gz"
  path "chr${i}.gvcf.gz.tbi"
  
  script:
  """
  for GVCF in ${gvcfs}
  do
    bcftools view \$GVCF --regions chr${i} -o \$GVCF.chr${i}.gvcf.gz -Oz
    tabix \$GVCF.chr${i}.gvcf.gz
    echo \$GVCF.chr${i}.gvcf.gz >> gvcf.list
  done
  gatk --java-options "-Xmx10g" CombineGVCFs \
  -R ${params.referenceGenome} \
  -V gvcf.list \
  -O chr${i}.gvcf.gz
  """
}

process jointGenotype {
  publishDir "${params.outDir}", mode: 'move'
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '12 GB'
  cpus 1
  maxRetries 18

  input:
  path gvcf
  
  output:
  path "${gvcf.SimpleName}.vcf.gz"
  
  script:
  """
  gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R ${params.referenceGenome}\
  -V gendb://${gvcf} \
  -O ${gvcf.SimpleName}.vcf.gz
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
  path gvcf_dir
  val i
  
  output:
  path "chr${i}.gvcf.gz"
  
  script:
  """
  bcftools view ${gvcf_file} --regions chr${i} -o chr${i}.gvcf.gz -Oz
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

process genomicsDBImport {
  containerOptions '--bind /groups/'
  errorStrategy 'retry'
  time '6h'
  memory '12 GB'
  cpus 4
  maxRetries 1

  input:
  path gvcfs
  val i

  output:
  path "chr${i}.gdb"
  
  script:
  """
  for GVCF in ${gvcfs}
  do
    tabix \$GVCF
    bcftools view \$GVCF --regions chr${i} -o \$GVCF.chr${i}.gvcf.gz -Oz
    echo "\${GVCF%%.*}_chr${i}\t\$GVCF" >> cohort.sample_map
  done
  gatk --java-options "-Xmx10g" GenomicsDBImport \
  --genomicsdb-workspace-path chr${i}.gdb \
  --batch-size 200 \
  --sample-name-map cohort.sample_map \
  -L chr${i} \
  --reader-threads 4
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
    haplotypeCaller(applyBQSR.output.bam_file, applyBQSR.output.index_file)
    gvcfs = haplotypeCaller.output.collect()
    chrs = Channel.from( 1..22 )
    genomicsDBImport(gvcfs, chrs)
    jointGenotype(genomicsDBImport.output)

    // indexJointGvcf(chrSplit.output)
    // jointGenotype(indexJointGvcf.output.gvcf_file, indexJointGvcf.output.tbi_file)
    // haplotypeCaller(applyBQSR.output.bam_file, applyBQSR.output.index_file)
    // indexGvcf(haplotypeCaller.output)
    // jointGenotype(indexGvcf.output.gvcf_file, indexGvcf.output.tbi_file)
}
