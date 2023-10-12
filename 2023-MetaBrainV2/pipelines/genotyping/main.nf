nextflow.enable.dsl=2

params.refFlat = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.collapsedGenes.refflat"
params.referenceGenome = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa"
params.gtfAnnotationFile = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz"
params.sampleFile = "/groups/umcg-biogen/tmp01/umcg-jbakker/samples.txt"
params.outDir = "/groups/umcg-biogen/tmp01/umcg-jbakker/results"

process splitNCigarReads {
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "${sample_bam.baseName}-splitreads.bam"
  
  script:
  """
  gatk --java-options "-Xmx8g" SplitNCigarReads \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O ${sample_bam.baseName}-splitreads.bam
  """
}

process haplotypeCaller {
  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "${sample_bam.baseName}.g.vcf.gz"
  
  script:
  """
  gatk --java-options "-Xmx8g" HaplotypeCaller \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O ${sample_bam.baseName}.g.vcf.gz
  -ERC GVCF
  """
}

workflow {
    def bamFiles = Channel.fromPath(params.bamDir)
    convertBAMToFASTQ(bamFiles)
}