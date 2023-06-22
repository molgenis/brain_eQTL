
library(data.table)
library(ggplot2)
library(plyr)
library("dplyr")
library(stringr)
library("optparse")
library(RColorBrewer)
library(tidyr)

# Get command line arguments 
option_list = list(
  make_option(c("-i", "--input"), type="character", default=getwd(), 
             help="path to input file base directory (will recursive search)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=getwd(),
                help="path to output dir", metavar="character"),
  make_option(c("-p", "--pcaOutliers"), type="character", default=NULL,
              help="path to file with PCA outliers", metavar="character"),
  make_option(c("-f", "--freeze"), type="character", default=NULL,
              help="freeze version, e.g. 2 or freeze2dot1", metavar="character"),
  make_option(c("-n", "--NABEC_pheno"), type="character", default=NULL,
              help="phenotype data for NABEC", metavar="character"),
  make_option(c("-a", "--plot_all"), type="character", default="TRUE",
              help="if set as false, don't make all the plots", metavar="character"),
  make_option(c("-s", "--expression_samples"), type="character",
              help="File with expression samples", metavar="character")
  
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

test <- function(){
  opt <- list()
  opt$input <- getwd()
  opt$output <- getwd()
  opt$expression_samples <- '2020-02-03-all-expression-samples.txt'
  opt$plot_all="FALSE"
  return(opt)
}
# comment out when not testing
#opt <- test()
dir.create(file.path(opt$output,'figures/' ), showWarnings = FALSE)
dir.create(file.path(paste0(opt$output,'/figures/'), "QC_figures_separate/" ), showWarnings = FALSE)

samples <- read.table(opt$expression_samples, header=F)

print(paste('Searching for input files in: ', opt$input))
if(!is.null(opt$pcaOutliers)){
  print(paste('Reading pca outlier file from: ',opt$pcaOutliers))
}
print(paste('Writing output files to: ',opt$output))


##### read in samples that are filtered out during PC #####
if(!is.null(opt$pcaOutliers)){
    pca_filtered_samples <- fread(opt$pcaOutliers,header=F)
}
#####
nabec_phenotypes <- NULL
#### Read in NABEC phenotype data ####
if(!is.null(opt$nabec_phenotypes)){
  nabec_phenotypes <- fread(opt$NABEC_pheno)
}
####

##### get cohort from path function that's being used a couple of times ####
get_cohort_from_path <- function(path){
  # order is important, goes from most to leas specific. 
  # e.g. CMC_HBCC needs to be before CMC, otherwise would return CMC due to grep
  if(grepl('ENA',path)) { return('ENA') }
  if(grepl('NABEC',path)) { return('NABEC') }
  if(grepl('GTEx',path)) { return('GTEx') }
  if(grepl('TargetALS',path)) { return('TargetALS') }
  if(grepl('Braineac',path)) { return('Braineac') }
  if(grepl('Brainseq',path)) { return('Brainseq') }
  if(grepl('CMC_HBCC',path)) { return('CMC_HBCC') }
  if(grepl('CMC',path)) { return('CMC') }
  if(grepl('BipSeq',path)) { return('BipSeq') }
  if(grepl('BrainGVEx',path)) { return('BrainGVEx') }
  if(grepl('UCLA_ASD',path)) { return('UCLA_ASD') }
  if(grepl('MayoTCX',path)) { return('MayoTCX') }
  if(grepl('MayoCBE',path)) { return('MayoCBE') }
  if(grepl('MSBB',path)) { return('MSBB') }
  if(grepl('ROSMAP',path)) { return('ROSMAP') }
  stop(paste('unknown cohort:',path))
}
#####


##### Read RnaMetrics alignment QC #######
print('Readin RNAseq QC metrics files')
RnaMetrics_multiQC_files <- list.files(path=opt$input, 
                                       pattern = "multiqc_picard_RnaSeqMetrics.txt$", 
                                       recursive = TRUE,
                                       full.names = T)
BPD <- RnaMetrics_multiQC_files[grepl("BPD", RnaMetrics_multiQC_files)]

RnaMetrics_multiQC_files <- RnaMetrics_multiQC_files[!grepl("BPD", RnaMetrics_multiQC_files)]
RnaMetric_qc_all <- data.frame()
for(f in RnaMetrics_multiQC_files){
  cohort <- get_cohort_from_path(f)
  if(cohort == 'ENA'){
    if(!grepl('multiQC_alignment',  f)) { next }
  }else{print(f)
    print(cohort)}
  
  RnaMetrics_qc <- fread(f)
  RnaMetrics_qc$cohort <- cohort
  RnaMetrics_qc$SAMPLE <- NULL
  RnaMetrics_qc$study <- ''
 if(cohort == 'ENA'){
    pattern <- ".*/multiqc_data/ENA/(.*)_batch[0-9]+.*"
    RnaMetrics_qc$study <- sub(pattern, "\\1", f)
  }
  RnaMetric_qc_all <- rbind(RnaMetric_qc_all, RnaMetrics_qc)
}
RnaMetric_qc_all <- data.frame(RnaMetric_qc_all)
RnaMetric_qc_all$Sample <- gsub('.cram','',RnaMetric_qc_all$Sample)

print(RnaMetric_qc_all[RnaMetric_qc_all$cohort!='BrainGVEx',])
RnaMetric_qc_all[RnaMetric_qc_all$cohort!='BrainGVEx',]$Sample <- gsub('-','.',RnaMetric_qc_all[RnaMetric_qc_all$cohort!='BrainGVEx',]$Sample)

if('TargetALS' %in% RnaMetric_qc_all$cohort){
  RnaMetric_qc_all[RnaMetric_qc_all$cohort=="TargetALS",]$Sample <- str_match(RnaMetric_qc_all[RnaMetric_qc_all$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
}
######

###### Read MultiMetrics QC #####
print('Readin MultipleMetrics files')
MultipleMetrics_multiQC_files <- list.files(path=opt$input, 
                                            pattern = "multiqc_picard_AlignmentSummaryMetrics.txt$", 
                                            recursive = TRUE,
                                            full.names = T)

MultipleMetrics_multiQC_files <- MultipleMetrics_multiQC_files[!grepl("BPD", MultipleMetrics_multiQC_files)]

MultipleMetric_qc_all <- data.frame()
for(f in MultipleMetrics_multiQC_files){
  cohort <- get_cohort_from_path(f)
  if(cohort == 'ENA'){
    if(!grepl('multiQC_alignment',  f)) { next }
  }
  MultipleMetrics_qc <- fread(f)
  MultipleMetrics_qc$cohort <- cohort
  MultipleMetrics_qc$SAMPLE <- NULL
  MultipleMetric_qc_all <- rbind(MultipleMetric_qc_all, MultipleMetrics_qc)
}
MultipleMetric_qc_all <- data.frame(MultipleMetric_qc_all)
MultipleMetric_qc_all$Sample <- gsub('.cram','',MultipleMetric_qc_all$Sample)
MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=='UCLA_ASD',]$Sample <- gsub('-','.',MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=='UCLA_ASD',]$Sample)

MultipleMetric_qc_all[MultipleMetric_qc_all$cohort!='BrainGVEx',]$Sample <- gsub('-','_',MultipleMetric_qc_all[MultipleMetric_qc_all$cohort!='BrainGVEx',]$Sample)

if('TargetALS' %in% MultipleMetric_qc_all$cohort){
  MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=="TargetALS",]$Sample <- str_match(MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
}
####### 


###### Read star QC #####
print('Reading star QC files')
STAR_multiQC_files <- list.files(path=opt$input, 
                                 pattern = "multiqc_star.txt$", 
                                 recursive = TRUE,
                                 full.names = T)
STAR_multiQC_files <- STAR_multiQC_files[!grepl("BPD", STAR_multiQC_files)]

STAR_qc_all <- data.frame()
for(f in STAR_multiQC_files){
  STAR_qc <- fread(f)
  STAR_qc$cohort <- get_cohort_from_path(f)
  STAR_qc_all <- rbind(STAR_qc_all, STAR_qc)
}
STAR_qc_all <- data.frame(STAR_qc_all)
STAR_qc_all$Sample <- gsub('.cram','',STAR_qc_all$Sample)
STAR_qc_all[RnaMetric_qc_all$cohort!='BrainGVEx',]$Sample <- gsub('-','.',STAR_qc_all[RnaMetric_qc_all$cohort!='BrainGVEx',]$Sample)

if('TargetALS' %in% STAR_qc_all$cohort){
  STAR_qc_all[STAR_qc_all$cohort=="TargetALS",]$Sample <- str_match(STAR_qc_all[STAR_qc_all$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
}
####### 

###### Read fastqc QC #####
print('Readin fastqc files')
FastQC_multiQC_files <- list.files(path=opt$input, 
                                   pattern = "multiqc_fastqc.txt$", 
                                   recursive = TRUE,
                                   full.names=T)
FastQC_multiQC_files <- FastQC_multiQC_files[!grepl("BPD", FastQC_multiQC_files)]
FastQC_qc_all <- data.frame()
for(f in FastQC_multiQC_files){
  FastQC_qc <- fread(f)
  if(grepl("MSBB",f)){
    # in MSBB there are resequenced samples and ssamples where the resequenced samples have been merged
    # with the original samples. We are using the merged ones, so remove the rest
    FastQC_qc <- FastQC_qc[!(grepl('esequen', FastQC_qc$Filename) & !grepl('merged', FastQC_qc$Filename)),]
  }
  FastQC_qc$cohort <- get_cohort_from_path(f)
  FastQC_qc_all <- rbind(FastQC_qc_all, FastQC_qc, fill=T)
}
# Because FastQC gives results per fastq file instead of per sample, adjust the names and merge them together
FastQC_qc_all <- data.frame(FastQC_qc_all)
FastQC_qc_all$Sample <- gsub('.cram','',FastQC_qc_all$Sample)


FastQC_qc_all_R1 <- FastQC_qc_all[grepl('*_R1$|*\\.R1$|*_1$|*\\.r1$', FastQC_qc_all$Sample),]
FastQC_qc_all_R1$Sample <- gsub('*_R1$|*\\.R1$|*_1$|*\\.r1$','',FastQC_qc_all_R1$Sample)
FastQC_qc_all_R2 <- FastQC_qc_all[grepl('*_R2$|*\\.R2$|*_2$|*\\.r2$', FastQC_qc_all$Sample),]
FastQC_qc_all_R2$Sample <- gsub('*_R2$|*\\.R2$|*_2$|*\\.r2$','',FastQC_qc_all_R2$Sample)

# everything that is leftover is single end, save as R1
single_end <- FastQC_qc_all[!grepl('*_R1$|*\\.R1$|*_1$|*\\.r1$', FastQC_qc_all$Sample),]
single_end <- single_end[!grepl('*_R2$|*\\.R2$|*_2$|*\\.r2', single_end$Sample),]
FastQC_qc_all_R1 <- rbind(FastQC_qc_all_R1, single_end)

# fill.x = True because single end sample are only in R1
FastQC_qc_all_merged <- merge(FastQC_qc_all_R1, FastQC_qc_all_R2, by=c('Sample','cohort'), suffix = c("_R1", "_R2"),
                              all.x=T)
FastQC_qc_all_merged$FastQC_original_sample <- FastQC_qc_all_merged$Sample

# For a few samples fastqc also got information on a 3rd file that contained orphan reads. Remove these (recognizable by lower total reads
FastQC_qc_all_merged <- FastQC_qc_all_merged[order(FastQC_qc_all_merged$Sample, -abs(FastQC_qc_all_merged$Total.Sequences_R1) ), ] #sort by id and reverse of abs(value)
FastQC_qc_all_merged <- FastQC_qc_all_merged[ !duplicated(FastQC_qc_all_merged$Sample), ]    

# TargetALS sample naming is different between FastQC and picard results, change them here to make them match later
if('TargetALS' %in% FastQC_qc_all_merged$cohort){
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- gsub('-','.',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample )
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- str_match(FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- gsub('-b38','',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample)
}
# Same for Brainseq
if('Brainseq' %in% FastQC_qc_all_merged$cohort){
    Brainseq_full_sample <- MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=="Brainseq",]$Sample
    Brainseq_fastqc_sample <- gsub('fastq','',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="Brainseq",]$Sample)

    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="Brainseq",]$Sample <- unlist(
                                                                            lapply(Brainseq_fastqc_sample, function(x) {
                                                                                Brainseq_full_sample[grepl(x, Brainseq_full_sample)]
                                                                            }))
}

####### 



###### merge FastQC and picard metrics together #####
RnaMetric_qc_all[grepl('PCT', colnames(RnaMetric_qc_all))] <- RnaMetric_qc_all[grepl('PCT', 
                                                                                colnames(RnaMetric_qc_all))]/100




CombinedMetrics <- merge(RnaMetric_qc_all, MultipleMetric_qc_all, by=c('Sample','cohort'),fill=T)
CombinedMetrics$READ_GROUP <- CombinedMetrics$READ_GROUP.x
CombinedMetrics$READ_GROUP.x <- NULL
CombinedMetrics$READ_GROUP.y <- NULL
CombinedMetrics$PF_ALIGNED_BASES <- CombinedMetrics$PF_ALIGNED_BASES.x
CombinedMetrics$PF_ALIGNED_BASES.x <- NULL
CombinedMetrics$PF_ALIGNED_BASES.y <- NULL
CombinedMetrics$library <- CombinedMetrics$library.x
CombinedMetrics$library.x <- NULL
CombinedMetrics$library.y <- NULL


CombinedMetrics <- merge(CombinedMetrics, STAR_qc_all, by=c('Sample','cohort'), all = TRUE)

if('Braineac' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="Braineac",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="Braineac",]$Sample, ".*(A653.*)")[, 2]
}
if('GTEx' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="GTEx",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="GTEx",]$Sample, "(.*)_.*")[, 2]
}
if('NABEC' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="NABEC",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="NABEC",]$Sample, ".*_(.*)")[, 2]
}
if('ENA' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="ENA",]$Sample <- str_match(
                                                                CombinedMetrics[CombinedMetrics$cohort=="ENA",]$Sample,
                                                                                                         ".*_(.*)")[, 2]
}
if('TargetALS' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample, ".*(HRA.*)")[, 2]
    CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample <- gsub('-2','',CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample)
}

if('BrainGVEx' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample <- gsub('_','-',CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample )
  CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample,"^([0-9]+-[0-9]+)-[0-9]+-[0-9]+")[, 2]
}
if('BipSeq' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="BipSeq",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="BipSeq",]$Sample,"Br[0-9]+_(R[0-9]+)")[, 2]
}

if('UCLA_ASD' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="UCLA_ASD",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="UCLA_ASD",]$Sample,"[aA-zZ]+[0-9]+_(.+)")[, 2]
  FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="UCLA_ASD" & grepl('-',FastQC_qc_all_merged$Sample),]$Sample <- gsub('-','.', FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="UCLA_ASD" & grepl('-',FastQC_qc_all_merged$Sample),]$Sample)
}


if('CMC_HBCC' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="CMC_HBCC",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="CMC_HBCC",]$Sample,"individualID.*_specimenID.(.*)")[, 2]
  
}
if('CMC' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="CMC",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="CMC",]$Sample,"CMC_[aA-zZ]+_[0-9]+_(.*)")[, 2]
}

if('MayoTCX' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="MayoTCX",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="MayoTCX",]$Sample,"^([0-9]+_TCX)_[0-9]+_TCX")[, 2]
}

if('MayoCBE' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="MayoCBE",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="MayoCBE",]$Sample,"^([0-9]+_CER)_[0-9]+_CER")[, 2]
  FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="MayoCBE",]$Sample <- str_match(FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="MayoCBE",]$Sample,"^([0-9]+_CER)*")[, 2]
}

if('MSBB' %in% CombinedMetrics$cohort){
  merged_sampeles <- gsub('_mergedWithResequenced','',CombinedMetrics[grepl('_mergedWithResequenced',CombinedMetrics$Sample),]$Sample)
  MSBB_to_remove <- CombinedMetrics[grepl(paste(merged_sampeles,collapse="|"), CombinedMetrics$Sample),]$Sample
  MSBB_to_remove <- MSBB_to_remove[!grepl('merged',MSBB_to_remove)]
  CombinedMetrics <- CombinedMetrics[!CombinedMetrics$Sample %in% MSBB_to_remove,]
  
  CombinedMetrics[CombinedMetrics$cohort=="MSBB",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="MSBB",]$Sample,"^AMPAD_MSSM_[0-9]+_(.*)")[, 2]
  CombinedMetrics[CombinedMetrics$cohort=="MSBB",]$Sample <- gsub('_mergedWithResequenced','',CombinedMetrics[CombinedMetrics$cohort=="MSBB",]$Sample)
  
  FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="MSBB",]$Sample <- gsub('.accepted_hits.sort.coord.combined','',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="MSBB",]$Sample)
  
}


if('ROSMAP' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="ROSMAP",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="ROSMAP",]$Sample,"^(.*_.*)_.*_.*")[, 2]
}
FastQC_qc_all_merged[grepl('hB_RNA_10492',FastQC_qc_all_merged$Sample),]


CombinedMetricsWithFastQC <- merge(CombinedMetrics, FastQC_qc_all_merged, by=c('Sample','cohort'), all = TRUE)
CombinedMetricsWithFastQC <- CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$Sample!='1271_TCX',] 

######



###### columns with % dont actually show %, have to multiply with 100
CombinedMetricsWithFastQC[which(grepl('PCT',colnames(CombinedMetricsWithFastQC)))] <- CombinedMetricsWithFastQC[which(grepl('PCT',colnames(CombinedMetricsWithFastQC)))]*100

###### set filter based on 3 measurements ######
CombinedMetricsWithFastQC$FILTER <- "NO"
print('set filter')


if(sum(!is.na(CombinedMetricsWithFastQC$PCT_CODING_BASES) & CombinedMetricsWithFastQC$PCT_CODING_BASES <= 10) > 0){
    CombinedMetricsWithFastQC[!is.na(CombinedMetricsWithFastQC$PCT_CODING_BASES) & CombinedMetricsWithFastQC$PCT_CODING_BASES <= 10,]$FILTER <- "YES"
}
if(sum(!is.na(CombinedMetricsWithFastQC$PCT_PF_READS_ALIGNED) & CombinedMetricsWithFastQC$PCT_PF_READS_ALIGNED <= 60) > 0){
    CombinedMetricsWithFastQC[!is.na(CombinedMetricsWithFastQC$PCT_PF_READS_ALIGNED) & CombinedMetricsWithFastQC$PCT_PF_READS_ALIGNED <= 60,]$FILTER <- "YES"
}


# sample should have been merged but resequenced sample didn't get transfered to us until the first one was already processed.
# only keep cov. data for teh correct one

CombinedMetricsWithFastQC <- CombinedMetricsWithFastQC[!duplicated(CombinedMetricsWithFastQC$Sample),]


filtered_samples_file <- paste0(opt$output,"/",Sys.Date(),"-samplesToFilter-",opt$freeze,".txt")
write.table(CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$FILTER=='YES',]$Sample, filtered_samples_file, quote=F, sep="\t", row.names=F)
print(paste0("Written samples to filter to ",filtered_samples_file))


nSamples <- CombinedMetricsWithFastQC %>% group_by(cohort) %>% dplyr::summarize(no = sum(FILTER=="NO"),yes = sum(FILTER=="YES"))
nSamples$total <- nSamples$no + nSamples$yes

#CombinedMetricsWithFastQC$SampleFull <- gsub('individualID.','', CombinedMetricsWithFastQC$SampleFull)
#CombinedMetricsWithFastQC$SampleFull <- gsub('specimenID.','', CombinedMetricsWithFastQC$SampleFull)

CombinedMetricsWithFastQC$pcaFiltered <- 'NO'

if(!is.null(opt$pcaOutliers)){
    # Add a separate column for those that got filtered with PCA
    CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$Sample %in% pca_filtered_samples$V1,]$pcaFiltered <- 'YES'
}
#####

##### Merge NABEC phenotype into CombinedMetricsWithFastQC
CombinedMetricsWithFastQC$lib_selection <- 'unknown'
if(!is.null(nabec_phenotypes)){
  if('NABEC' %in% CombinedMetricsWithFastQC$cohort){
    nabec_phenotypes$Sample <- paste0(nabec_phenotypes$BioSample, '_', nabec_phenotypes$Run)
    CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$cohort=="NABEC",]$lib_selection <- nabec_phenotypes[match(CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$cohort=="NABEC",]$Sample,
                                                                             nabec_phenotypes$Sample),]$LibrarySelection
  }
}

#####

##### plot function so that it can be more easily repeated
QC_plotter <- function(CombinedMetricsWithFastQC, column, plot_pca_outliers = F, single_cohort =F, lib_selection=F){
  # add the super small nuber in case y-axis needs to be log scaled, but only for those columns where the max value > 100000 (so that columns with e.g. percentages don't get log scaled)
  if(max(CombinedMetricsWithFastQC[!is.na(CombinedMetricsWithFastQC[column]),][column])> 10000){
    CombinedMetricsWithFastQC[!is.na(CombinedMetricsWithFastQC[column]),][column] <-CombinedMetricsWithFastQC[!is.na(CombinedMetricsWithFastQC[column]),][column] +0.00000000000000000000001
  }
  
  # Because we have more than 8 colours, have to extend the Dark2 colour palette
  colourCount = length(unique(CombinedMetricsWithFastQC$cohort))
  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
  dark2ExtendedPallete <- getPalette(colourCount)
  
  p <- ggplot()
  
  if(single_cohort){
    p <- p + geom_point(data=CombinedMetricsWithFastQC, aes_string('sampleID', column, colour='FILTER'))
  }else{
    # if PCA outliers are plotted, make all others have lower alpha settings
    if(plot_pca_outliers){
      if(lib_selection){
        p <- p + geom_point(alpha=0.1)+ 
          geom_point(data=CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$pcaFiltered=='YES',],  aes_string('sampleID', column, colour='cohort', shape='lib_selection'))
      }else{
        p <- p + geom_point(alpha=0.1)+ 
             geom_point(data=CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$pcaFiltered=='YES',],  aes_string('sampleID', column, colour='cohort', shape='FILTER'))
      }
    }else{
      if(lib_selection){
        p <- p + geom_point(data=CombinedMetricsWithFastQC, aes_string('sampleID', column, colour='cohort', shape='lib_selection'))
      }else{
        p <- p + geom_point(data=CombinedMetricsWithFastQC, aes_string('sampleID', column, colour='cohort', shape='FILTER'))
      }
    }
  }
  p <- p + theme_bw(base_size=30)+
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank())+  
    xlab('Samples')+ 
    guides(colour = guide_legend(override.aes = list(size=10)),
           shape = guide_legend(override.aes = list(size=10)))
  
  if(single_cohort){
    p <- p+scale_colour_brewer(palette="Dark2")
  }else{
    p <- p+scale_colour_manual(values=dark2ExtendedPallete)
  }

  # Add the log scaling for columns with high values
  if(max(CombinedMetricsWithFastQC[!is.na(CombinedMetricsWithFastQC[column]),][column])> 10000){
    p <- p + scale_y_continuous(trans='log10') +
      ylab(paste0('log10( ',column, ' )'))
  }
  
  # Plot a line where we put the threshold
  if(column == "PCT_CODING_BASES"){
    p <- p + geom_hline(yintercept=10, colour="red",  linetype="dashed", size=2)
  }
  if(column == "PCT_PF_READS_ALIGNED"){
    p <- p + geom_hline(yintercept=60, colour="red",  linetype="dashed", size=2)
  }
  if(column == "uniquely_mapped_percent"){
    p <- p + geom_hline(yintercept=60, colour="red",  linetype="dashed", size=2)
  }
  return(p)
}
#####

##### plot STAR + multimetrics together #####

CombinedMetricsWithFastQC <- CombinedMetricsWithFastQC[order(CombinedMetricsWithFastQC$cohort),]
CombinedMetricsWithFastQC$sampleID <- 1:nrow(CombinedMetricsWithFastQC)

# make sure that the columns used for filtering are plotted first so that they are easy to find
columns_to_plot <- colnames(select_if(CombinedMetricsWithFastQC, is.numeric))
filter_columns <-  c("PCT_CODING_BASES","PCT_PF_READS_ALIGNED")
columns_to_plot <- columns_to_plot[!columns_to_plot %in% filter_columns]
columns_to_plot <- c(filter_columns, columns_to_plot)
CombinedMetricsWithFastQC <- CombinedMetricsWithFastQC[!is.na(CombinedMetricsWithFastQC$Sample),]
if(opt$plot_all == "TRUE"){
    print('start plotting')
 
   pdf(paste0(opt$output,'/figures/all_MultiMetrics_QC_plots.pdf'), width=8, height=8)
    for(column in columns_to_plot){
      if(column == "sampleID"){
        next
      }
      print(column)
      p <- QC_plotter(CombinedMetricsWithFastQC, column, FALSE)
      ggsave(paste0(opt$output,'/figures/QC_figures_separate/',column,'.pdf'), width=12, height = 9, plot=p)
      print(p)
    }
    dev.off()
}

colourCount = length(unique(CombinedMetricsWithFastQC$cohort))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
dark2ExtendedPallete <- getPalette(colourCount)

ggplot(CombinedMetricsWithFastQC[CombinedMetricsWithFastQC$FILTER=="NO",], aes(cohort,fill=cohort,colour='yes'))+
  geom_bar()+ 
  theme_bw(base_size=30)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=dark2ExtendedPallete)+
  scale_colour_manual(values='black')+
  guides(colour=F)+ 
  coord_flip()+
  guides(fill=guide_legend(ncol=2)) + theme(legend.position = "none")
ggsave(paste0(opt$output,'/figures/n_RNAseq_samples.pdf'),width=8, height=12)
#####

##### plot lib select #####
plot_with_lib_selection <- function(){
  # used for debugging nabec, not necesarry to run
  print('plot library selection')
  pdf(paste0(opt$output,'/figures/all_MultiMetrics_QC_plots_library_selection.pdf'), width=8, height=8)
  for(column in columns_to_plot){
    if(column == "sampleID"){
      next
    }
    print(column)
    p <- QC_plotter(CombinedMetricsWithFastQC, column, FALSE, lib_selection=T)
    print(p)
  }
  dev.off()
}
# plot_with_lib_selection()
#####

##### plot STAR + multimetrics again, but highlight PCA filtered samples #####
plot_with_pca_outliers <- function(){
  if(!is.null(opt$pcaOutliers)){
      print('plotting with pca outliers')
      pdf(paste0(opt$output,'/figures/all_MultiMetrics_QC_plots.highlightPcaFilteredSamples.pdf'), width=8, height=8)
      for(column in columns_to_plot){
        if(column == "sampleID"){
          next
        }
        print(column)
        p <- QC_plotter(CombinedMetricsWithFastQC, column, TRUE)
        ggsave(paste0(opt$output,'/figures/QC_figures_separate/',column,'.highlightPcaFilteredSamples.pdf'), width=12, height = 8, plot=p)
        print(p)
      }
      dev.off()
  }
}
samples$V1 <- as.character(samples$V1)
missing <- samples$V1[!samples$V1 %in% CombinedMetricsWithFastQC$Sample]
missing <- missing[!grepl('_resequenced', missing)]
missing <- missing[!grepl('merged', missing)]
if(length(missing) > 0){
    cat("Missing samples:\n")
    print(missing)
    cat("\n")
}
write.table(CombinedMetricsWithFastQC,file=paste0(opt$output,'/',Sys.Date(),'-',opt$freeze,'.TMM.Covariates.txt'),quote=F, 
    row.names=FALSE, sep='\t')
