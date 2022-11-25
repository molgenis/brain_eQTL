library("optparse")
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(plyr)
library(data.table)
# Get command line arguments 
option_list = list(
  make_option(c("-m", "--MetaBrain"), type="character",
              help="bin file for MetaBrain", metavar="character"),
  make_option(c("-e", "--eqtlGen"), type="character",
              help="bin files for eqtlGen", metavar="character"),
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- list()
opt$MetaBrain <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/2019-06-05-metabrain-eQTLgenGenes-bins.txt"
opt$eqtlGen <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/BIOS_expression_summary_bins_20180515.txt"


metaBrain <- fread(opt$MetaBrain)
eqtlGen <- fread(opt$eqtlGen)

diff_per_bin <- data.frame()
diff_per_bin_combined <- data.frame()
for(bin_number in unique(metaBrain$bin_order)){
  eqtl_gen_bin_name <- paste0('bin',bin_number)
  
  metaBrain_current_bin <- metaBrain[metaBrain$bin_order==bin_number,]
  eqtlGen_current_bin <- eqtlGen[eqtlGen$eprs_exp_mean_bin==eqtl_gen_bin_name,]
  
  eqtlGen_current_bin_cis_gene  <- table(eqtlGen_current_bin$cis_eGene)
  metaBrain_current_bin_cis_gene  <- table(metaBrain_current_bin$is_cis_qtl)
  
  
  eqtlGen_ratio <- eqtlGen_current_bin_cis_gene[['yes']] / (eqtlGen_current_bin_cis_gene[['yes']] + eqtlGen_current_bin_cis_gene[['no']])
  metaBrain_ratio <- metaBrain_current_bin_cis_gene[['TRUE']] / (metaBrain_current_bin_cis_gene[['TRUE']] + metaBrain_current_bin_cis_gene[['FALSE']])
  df <- data.frame('bin'=bin_number,'ratio'=eqtlGen_ratio,dataset='eqtlGen')
  diff_per_bin <- rbind(diff_per_bin,df)
  
  df <- data.frame('bin'=bin_number,'ratio'=metaBrain_ratio,dataset='metaBrain')
  diff_per_bin <- rbind(diff_per_bin,df)

  df <- data.frame('bin'=bin_number,'metaBrain_ratio'=metaBrain_ratio, 'eQTLgen_ratio'=eqtlGen_ratio)
  diff_per_bin_combined <- rbind(diff_per_bin_combined, df)
}

diff_per_bin_combined$diff <- diff_per_bin_combined$metaBrain_ratio-diff_per_bin_combined$eQTLgen_ratio


ggplot(diff_per_bin, aes(x=bin, y=ratio, fill=dataset))+
  geom_bar(stat='identity', position = 'dodge')+
  theme_bw(base_size = 16)

diff_per_bin_combined
