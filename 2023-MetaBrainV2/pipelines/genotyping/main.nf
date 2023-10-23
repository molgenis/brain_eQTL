nextflow.enable.dsl=2

params.refFlat = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.collapsedGenes.refflat"
params.referenceGenome = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa"
params.gtfAnnotationFile = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz"
params.sampleFile = "/groups/umcg-biogen/tmp01/umcg-jbakker/samples.txt"
params.outDir = "/groups/umcg-biogen/tmp01/umcg-jbakker/results"
params.bamDir1 = "/groups/umcg-biogen/tmp01/input/rawdata/2023-TargetALS/RNA_BAM/CGND_12616/Project_CGND_12616_B03_EXS_Lane.bam.2018-09-10/Sample_CGND-HRA-00076-v2-2-grch38/analysis/CGND-HRA-00076-v2-2-grch38.final.bam"
params.bamDir = "/groups/umcg-fg/tmp01/projects/downstreamer/2023-09-PublicRNASeqGenotypeCalls/metabrain/2023-MetaBrainV2/pipelines/genotyping/bam"

def getEncoding(bam_file) {
  """
  fastqc ${bam_file} --outdir fastqc
  unzip fastqc/*.zip -d unzipped
  unzipped/*/fastqc_data.txt
  sed '6!d' unzipped/*/fastqc_data.txt | cut -d " " -f 4
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

process fixMisencodedBaseQualityReads {
  containerOptions '--bind /groups/'
  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path sample_bam
  
  output:
  path "fixed/${sample_bam.baseName}-fixed.bam"
  
  shell:
  """
  mkdir fixed
  gatk --java-options "-Xmx16g" FixMisencodedBaseQualityReads \
  -I ${sample_bam} \
  -O fixed/${sample_bam.baseName}-fixed.bam \
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
  
  shell:
  """

  samtools view ${markDuplicatesBam} | \
    head -1000000 | \
    awk '{gsub(/./,"&\n",$11);print $11}'| \
    sort -u| \
    perl -wne '
    $_=ord($_);
    print $_."\n"if(not($_=~/10/));' | \
    sort -n | \
    perl -wne '
    use strict;
    use List::Util qw/max min/;
    my @ords=<STDIN>;
    if(min(@ords) >= 59 && max(@ords) <=104 ){
      print " --fix_misencoded_quality_scores ";
      warn "Illumina <= 1.7 scores detected using:--fix_misencoded_quality_scores.\n";
    }elsif(min(@ords) >= 33 && max(@ords) <= 74){
	  print " ";
	  warn "quals > illumina 1.8 detected no action to take.\n";
    }elsif(min(@ords) >= 33 && max(@ords) <= 80){
	  print " --allow_potentially_misencoded_quality_scores ";
	  warn "Strange illumina like quals detected using:--allow_potentially_misencoded_quality_scores."
    }else{
	  die "Cannot estimate quality scores here is the list:".join(",",@ords)."\n";
    }
  ')


  mkdir split
  gatk --java-options "-Xmx16g" SplitNCigarReads \
  -R ${params.referenceGenome}\
  -I ${sample_bam} \
  -O split/${sample_bam.baseName}-splitreads.bam \
  --fix_misencoded_quality_scores
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


workflow {
    bam_files = Channel.fromPath("${params.bamDir}/*.bam")
    splitNCigarReads(bam_files) 
    haplotypeCaller(splitNCigarReads.output)
    indexGvcf(haplotypeCaller.output)
    combineGvcf(haplotypeCaller.output)
    
    //indexGvcf(gvcf_files)
}