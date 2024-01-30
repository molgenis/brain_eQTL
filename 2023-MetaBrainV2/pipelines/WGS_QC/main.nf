nextflow.enable.dsl=2

process mergeVCFFiles {
  containerOptions "--bind ${params.bindFolder}"

  time '12h'
  memory '32 GB'
  cpus 1

  input:
  path sampleFile

  output:
  path "merged.vcf.gz"

  shell:
  '''
  bcftools merge -0 -l !{sampleFile} -Oz -o merged.vcf.gz
  '''
}

process createGteFile {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '1 GB'
  cpus 1
  
  input:
  path vcfFile
  
  output:
  path "gte.txt"
  
  shell:
  '''
  samples=$(bcftools query -l !{vcfFile})

  echo "$samples" | while read sample; do
      echo -e "$sample\t$sample"
  done > "gte.txt"
  '''
}

process splitByChromosome {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path vcfFile
  
  output:
  path "*.vcf.gz"
  
  shell:
  '''
  # 1. Create index file of merged VCF
  tabix -p vcf !{vcfFile}

  # 2. Define chromosomes that should be included
  chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")

  # 3. Split the merged VCF file by chromosome
  bcftools index -s !{vcfFile} | cut -f 1 | while read C; do
    if [[ " ${chromosomes[*]} " =~ ${C#"chr"} ]]; then
      bcftools view -Oz -o split.${C}.vcf.gz !{vcfFile} ${C}
    fi
  done
  '''
}

process removeMultiAllelicVariants {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path vcfFile
  
  output:
  path "no_multi_allelic.*.vcf.gz"
  
  shell:
  '''
  # 1. Extract chromosome number from file name
  IFS='.' read -r -a fileNameParts <<< "!{vcfFile}"
  chromosome="${fileNameParts[1]}"

  # 2. Split multi-allelic variants
  #bcftools norm -m -any !{vcfFile} -Oz -o "norm.${chromosome}.vcf.gz"

  # 2. Remove multi allelic variants
  bcftools view --max-alleles 2 !{vcfFile} -Oz -o "no_multi_allelic.${chromosome}.vcf.gz"
  '''
}

process wgsQC {
  containerOptions "--bind ${params.bindFolder}"
  publishDir "${params.outDir}/0_filter_vcf/", mode: 'copy', pattern: "*.{log,txt}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path vcfFile
  
  output:
  path "*-filtered.vcf.gz", emit: vcfFile
  path "*.log"
  path "*.txt"
  
  shell:
  '''
  # 1. Extract chromosome number from file name
  IFS='.' read -r -a fileNameParts <<< "!{vcfFile}"
  chromosome="${fileNameParts[1]}"

  # 2. Define command arguments
  commandArguments="--input !{vcfFile} \
  --output no_multi_allelic.${chromosome} \
  --sex !{params.sexFile}
  --ignore_homref_stats"

  # 3. If the chromosome does not contain any numbers (so it's a sex chromosome), add hardy weinberg argument
  if ! [[ ${chromosome} =~ [0-9] ]]; then
    commandArguments+=" --hardy_weinberg_equilibrium 0" 
  fi

  # 4. Run the custom VCF filter script
  python3 !{baseDir}/scripts/custom_vcf_filter.py ${commandArguments} \
  | tee custom_vcf_filter.log

  # 5. Create filter log
  python3 !{baseDir}/scripts/print_WGS_VCF_filter_overview.py \
    --workdir . --chr ${chromosome} \
    --vcf_file_format "no_multi_allelic.${chromosome}.vcf.gz"
  '''
}

process convertVCFToPlink {
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path vcfFile
  
  output:
  path "plink_out_*/*"
  
  shell:
  '''
  # 1. Extract chromosome number from file name
  IFS='.' read -r -a fileNameParts <<< "!{vcfFile}"
  chromosome="${fileNameParts[1]}"

  # 2. Make output dir
  mkdir plink_out_${chromosome}
  
  # 3. Convert to pgen (because of --sort-vars)
  ~/plink2 --vcf !{vcfFile}  \
  --const-fid \
  --make-pgen \
  --output-chr 26 \
  --set-missing-var-ids @:#[b38]\\$r,\\$a \
  --new-id-max-allele-len 1000 \
  --sort-vars \
  --out "plink_out_${chromosome}/${chromosome}_converted_vcf"

  # 4. Convert pgen to bed/bim/fam
  ~/plink2 --pfile "plink_out_${chromosome}/${chromosome}_converted_vcf" --make-bed --out "plink_out_${chromosome}/${chromosome}_converted_vcf"
  '''
}

process mergePlinkFiles {
  publishDir "${params.outDir}/1_merge_plink/", mode: 'copy'
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path plinkDir
  
  output:
  path 'chrAll.bed', emit: bedFile
  path 'chrAll.bim', emit: bimFile
  path 'chrAll.fam', emit: famFile
  
  shell:
  '''
  # 1. Split input filepath string into separate paths
  IFS='\s' read -r -a plinkFiles <<< "!{plinkDir}"
  
  # 2. Create empty arrays for the different file types
  bed_files=()
  bim_files=()
  fam_files=()

  # 3. Loop through the split file paths and add them to the corresponding array
  for file in "${plinkFiles[@]}"; do
      if [[ $file == *".bed" ]]; then
          bed_files+=("$file")
      elif [[ $file == *".bim" ]]; then
          bim_files+=("$file")
      elif [[ $file == *".fam" ]]; then
          fam_files+=("$file")
      fi
  done

  # 4. Loop through the arrays and add bed/bim/bam paths to merge file
  for i in "${!bed_files[@]}"; do
      if [[ -n "${bed_files[$i]}" && -n "${bim_files[$i]}" && -n "${fam_files[$i]}" ]]; then
        if [ -s "${bim_files[$i]}" ]; then
          echo "${bed_files[$i]} ${bim_files[$i]} ${fam_files[$i]}" >> merge_list.txt
        fi
      fi
  done

  # 5. Run plink merge command to merge files
  ~/plink2 --pmerge-list merge_list.txt --merge-par --make-bed --out "chrAll" 
  '''
}

process calculateAlleleFrequencies {
  publishDir "${params.outDir}/2_calculate_allele_frequencies/", mode: 'copy', pattern: "*.{log, gz}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1
  
  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path 'target.afreq.gz'
  
  script:
  """
  # 1. Get allele frequencies
  ~/plink2 --bed ${bedFile} \
    --bim ${bimFile} \
    --fam ${famFile} \
    --threads 4 \
    --freq 'cols=+pos' \
    --out target

  # 2. Gzip output
  gzip target.afreq --force
  """
}

process plotTargetAndRefAlleleFrequencies {
  publishDir "${params.outDir}/3_plot_allele_frequencies/", mode: 'copy', pattern: "*.{png}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path targetAFs
  
  output:
  path 'allele_frequencies.png'
  
  script:
  """
  Rscript ${baseDir}/scripts/plot_allele_frequencies.R --ref_afs ${params.refAFs} --target_afs ${targetAFs}
  """
}

process snpAndGenotypeQC {
  publishDir "${params.outDir}/4_snp_and_genotype_qc/", mode: 'copy', pattern: "*.{log}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path 'chrAll_QC.bed', emit: bedFile
  path 'chrAll_QC.bim', emit: bimFile
  path 'chrAll_QC.fam', emit: famFile

  shell:
  '''
  ~/plink2 --bed !{bedFile} \
    --bim !{bimFile} \
    --fam !{famFile} \
    --maf !{params.maf} \
    --geno !{params.geno} \
    --mind !{params.mind} \
    --hwe !{params.hwe} \
    --make-bed \
    --out chrAll_QC \
    --output-chr 26 \
    --not-chr 0 25-26 \
    --set-all-var-ids @:#[b38]\\$r,\\$a \
    --new-id-max-allele-len 10 truncate \
    --threads 4
  '''
}

process sexCheck {
  publishDir "${params.outDir}/5_sex_check/", mode: 'copy', pattern: "*.{txt,pdf,png,log}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '12 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path 'SexCheck.bed', emit: bedFile
  path 'SexCheck.bim', emit: bimFile
  path 'SexCheck.fam', emit: famFile
  path '*.png'
  path '*.txt'
  path '*.pdf'
  path '*.log'
  
  script:
  """
  # 1. Split sex chromosome
  ~/plink2 --bed ${bedFile} \
      --bim ${bimFile} \
      --fam ${famFile} \
      --chr X \
      --maf ${params.maf} \
      --split-par b38 \
      --make-bed \
      --out chrAll_split

  # 2. Pruning
  ~/plink2 --bfile chrAll_split \
      --rm-dup 'exclude-mismatch' \
      --indep-pairwise 20000 200 0.2 \
      --out check_sex_x \
      --threads 4

  # 3. Sex check
  plink2 --bfile 'chrAll_split' \
      --extract check_sex_x.prune.in \
      --check-sex 0.4 0.6 \
      --threads 4

  # 5. Plot results
  Rscript ${baseDir}/scripts/plot_sex_check.R --sex_check_file plink.sexcheck

  # 6. Remove samples that failed sex check
  ~/plink2 --bed ${bedFile} \
    --bim ${bimFile} \
    --fam ${famFile} \
    --maf ${params.maf} \
    --geno ${params.geno} \
    --mind ${params.mind} \
    --hwe ${params.hwe} \
    --autosome \
    --make-bed \
    --out SexCheck \
    --output-chr 26 \
    --remove SexCheckFailed.txt \
    --threads 4
  """
}

process heterozygosityCheck {
  publishDir "${params.outDir}/6_het_check/", mode: 'copy', pattern: "*.{txt,pdf,png,log}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  path alleleFrequencies
  
  output:
  path 'HetCheck.bed', emit: bedFile
  path 'HetCheck.bim', emit: bimFile
  path 'HetCheck.fam', emit: famFile
  path '*.png'
  path '*.pdf'
  path '*.txt'
  path '*.log'
  
  script:
  """
  # 1. Pruning
  ~/plink2 --bed ${bedFile} \
      --bim ${bimFile} \
      --fam ${famFile} \
      --rm-dup 'exclude-mismatch' \
      --indep-pairwise 50 1 0.2 \
      --threads 4

  # 2. Heterozygosity check
  ~/plink2 --bed ${bedFile} \
      --bim ${bimFile} \
      --fam ${famFile} \
      --extract plink2.prune.in \
      --het \
      --threads 4 \
      --read-freq ${alleleFrequencies} \
      --out het_check_results

  # 3. Plot results
  Rscript ${baseDir}/scripts/plot_het_check.R --het_check_file het_check_results.het

  # 4. Remove failed samples
  ~/plink2 --bed ${bedFile} \
    --bim ${bimFile} \
    --fam ${famFile} \
    --remove HeterozygosityFailed.txt \
    --make-bed \
    --out HetCheck
  """
}

process projectSamplesToRefPanel {
  publishDir "${params.outDir}/7_project_samples_to_ref_panel/", mode: 'copy', pattern: "*.{txt,pdf,png,log}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path '*.png'
  path '*.pdf'
  path '1000G_PC_projections.txt'
  path 'PopAssignResults.txt'
  
  script:
  """
  Rscript ${baseDir}/scripts/project_samples_to_superpop.R --ref_bed ${params.refPath}.bed \
    --target_bed ${bedFile} --ref_pop ${params.refPop}
  """
}

process relatednessCheck {
  publishDir "${params.outDir}/8_relatedness_check/", mode: 'copy', pattern: "*.{log,txt,kin0}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path 'RelatednessCheck.bed', emit: bedFile
  path 'RelatednessCheck.bim', emit: bimFile
  path 'RelatednessCheck.fam', emit: famFile
  path 'related.kin0'
  path 'RelatednessPassedSamples.txt'
  
  script:
  """
  # 1. Do relatedness check
  ~/plink2 --bed ${bedFile} \
    --bim ${bimFile} \
    --fam ${famFile} \
    --make-king-table \
    --king-table-filter ${params.kingTableFilter} \
    --out related \
    --threads 4 \

  # 2. Create sample list of non-related samples
  Rscript ${baseDir}/scripts/find_related_samples.R --kin_file related.kin0 --target_bed ${bedFile}

  # 3. Remove samples that are not on the list created above
  ~/plink2 --bed ${bedFile} \
    --bim ${bimFile} \
    --fam ${famFile} \
    --keep RelatednessPassedSamples.txt \
    --make-bed \
    --out RelatednessCheck
  """
}

process targetPCA {
  publishDir "${params.outDir}/9_target_pca/", mode: 'copy', pattern: "*.{txt,pdf,png,log}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '16 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path 'toImputation.bed', emit: bedFile
  path 'toImputation.bim', emit: bimFile
  path 'toImputation.fam', emit: famFile
  path 'SamplesToInclude.txt', emit: sampleFile
  
  script:
  """
  # 1. Do PCA and find outliers
  Rscript ${baseDir}/scripts/target_pca.R --target_bed ${bedFile}

  # 2. Remove outlier samples
  ~/plink2 --bed ${bedFile} \
    --bim ${bimFile} \
    --fam ${famFile} \
    --output-chr 26 \
    --keep SamplesToInclude.txt \
    --make-bed \
    --threads 4 \
    --out toImputation
  """
}

process shuffleSampleOrder {
  publishDir "${params.outDir}/10_shuffle_sample_order/", mode: 'copy', pattern: "*.{txt,pdf,png,log}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '1 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  path sampleFile
  
  output:
  path 'shuffled.bed', emit: bedFile
  path 'shuffled.bim', emit: bimFile
  path 'shuffled.fam', emit: famFile
  
  script:
  """
  # 1. Create shuffled sample list
  Rscript ${baseDir}/scripts/create_shuffled_sample_list.R --sample_file ${sampleFile}

  # 2. Shuffle samples 
  ~/plink2 --bed ${bedFile} \
    --bim ${bimFile} \
    --fam ${famFile} \
    --indiv-sort f ShuffledSampleOrder.txt \
    --make-bed \
    --out shuffled \
    --threads 4
  """
}

process finalSNPandGenotypeQC {
  publishDir "${params.outDir}/11_final_snp_and_genotype_qc/", mode: 'copy', pattern: "*.{txt,pdf,png,log}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '8 GB'
  cpus 1

  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path 'chrAll_QC.bed'
  path 'chrAll_QC.bim'
  path 'chrAll_QC.fam'
  
  shell:
  '''
  ~/plink2 --bed !{bedFile} \
    --bim !{bimFile} \
    --fam !{famFile} \
    --maf !{params.maf} \
    --geno !{params.geno} \
    --mind !{params.mind} \
    --hwe !{params.hwe} \
    --autosome \
    --make-bed \
    --out chrAll_QC \
    --output-chr 26 \
    --not-chr 0 25-26 \
    --set-all-var-ids @:#[b38]\\$r,\\$a \
    --new-id-max-allele-len 10 truncate \
    --threads 4
  '''
}

process finalPCA {
  publishDir "${params.outDir}/12_final_pca/", mode: 'copy', pattern: "*/*.{txt,pdf,png}"
  containerOptions "--bind ${params.bindFolder}"

  time '6h'
  memory '16 GB'
  cpus 1

  publishDir "${params.outDir}/het_check/", mode: 'copy', pattern: "*/*.{txt,pdf,png}"

  input:
  path bedFile
  path bimFile
  path famFile
  
  output:
  path '*.png'
  path '*.pdf'
  path '*.txt'
  
  script:
  """
  Rscript ${baseDir}/scripts/final_pca.R --target_bed ${bedFile} 
  """
}


workflow {
    mergeVCFFiles(params.sampleFile)
    createGteFile(mergeVCFFiles.out)
    splitByChromosome(mergeVCFFiles.out) | flatten | removeMultiAllelicVariants
    wgsQC(removeMultiAllelicVariants.out)
    convertVCFToPlink(wgsQC.out.vcfFile) | collect | mergePlinkFiles
    calculateAlleleFrequencies(
      mergePlinkFiles.out.bedFile,
      mergePlinkFiles.out.bimFile,
      mergePlinkFiles.out.famFile)
    plotTargetAndRefAlleleFrequencies(calculateAlleleFrequencies.out)
    snpAndGenotypeQC(
      mergePlinkFiles.out.bedFile,
      mergePlinkFiles.out.bimFile,
      mergePlinkFiles.out.famFile)
    sexCheck(
      snpAndGenotypeQC.out.bedFile,
      snpAndGenotypeQC.out.bimFile,
      snpAndGenotypeQC.out.famFile
      )
    heterozygosityCheck(
      sexCheck.out.bedFile, 
      sexCheck.out.bimFile, 
      sexCheck.out.famFile, 
      calculateAlleleFrequencies.out
      )
    projectSamplesToRefPanel(
      heterozygosityCheck.out.bedFile,
      heterozygosityCheck.out.bimFile,
      heterozygosityCheck.out.famFile
      )
    relatednessCheck(
      heterozygosityCheck.out.bedFile,
      heterozygosityCheck.out.bimFile,
      heterozygosityCheck.out.famFile
      )
    targetPCA(
      relatednessCheck.out.bedFile,
      relatednessCheck.out.bimFile,
      relatednessCheck.out.famFile
      )
    shuffleSampleOrder(
      targetPCA.out.bedFile, 
      targetPCA.out.bimFile, 
      targetPCA.out.famFile, 
      targetPCA.out.sampleFile
      )
    finalSNPandGenotypeQC(
      shuffleSampleOrder.out.bedFile,
      shuffleSampleOrder.out.bimFile,
      shuffleSampleOrder.out.famFile
      )
    finalPCA(
      shuffleSampleOrder.out.bedFile,
      shuffleSampleOrder.out.bimFile,
      shuffleSampleOrder.out.famFile
      )  
}