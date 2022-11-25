library("optparse")
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(plyr)
library(data.table)
# Get command line arguments 
option_list = list(
  make_option(c("-c", "--cisQTLs"), type="character",
              help="eQTLProbesFDR0.05-ProbeLevel.txt.gz of cis-eQTLs", metavar="character"),
  make_option(c("-e", "--expression"), type="character",
              help="File with mean expression per gene", metavar="character"),
  make_option(c("-p", "--proteinCodingGenes"), type="character",
              help="File with protein coding genes", metavar="character"),
  make_option(c("-i", "--pLI"), type="character",
              help="File with pLI scores", metavar="character"),
  make_option(c("-o", "--out_path"), type="character",
              help="Out path for plots", metavar="character"),
  make_option(c("-z", "--ensemble_IDs"), type="character",
              help="Path to file with ensembl IDS", metavar="character"),

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- list()
opt$cisQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/eQTLProbesFDR0.05-ProbeLevel.CIS.txt.gz"
opt$expression <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.mean_and_sd.txt"
opt$proteinCodingGenes <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/protein_coding_genes_ensembl84.txt"
opt$pLI <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"
opt$out_path <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/figures/freeze2/"
opt$ensemble_IDs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/ensembl_transcript_gene_IDs.txt"


##### read in data ####
cis_qtl <- fread(cmd=paste0("gzip -dc ",opt$cisQTLs), header=T, sep='\t')

protein_coding_genes <- read.table(opt$proteinCodingGenes, sep='\t', header=T)
pLI <- read.table(opt$pLI, sep='\t', header=T)
ensemble_IDs <- fread(opt$ensemble_IDs)
expression <- read.table(opt$expression, header=T, sep='\t', row.names=1)
expression$signal_to_nois <- expression$SD/expression$mean_expression

# make sure all are FDR < 0.05
cis_qtl <- cis_qtl[cis_qtl$FDR < 0.05,]

# genes still have version numbers, remove
rownames(expression) <- sapply(strsplit(rownames(expression),"\\."), `[`, 1)
cis_qtl$ProbeName <- sapply(strsplit(as.character(cis_qtl$ProbeName),"\\."), `[`, 1)

# select only protein coding genes
expression_protCoding <- expression[rownames(expression) %in% protein_coding_genes$Ensembl.Gene.ID,]
cis_qtl_protCoding <- cis_qtl[cis_qtl$ProbeName %in% protein_coding_genes$Ensembl.Gene.ID,]
#####

##### make expression and sd bins bins #####
expression$expression_bin <- cut_number(expression$mean_expression, 10)
expression_protCoding$expression_bin <- cut_number(expression_protCoding$mean_expression, 10)
expression$sd_bin <- cut_number(expression$SD, 10)
expression_protCoding$sd_bin <- cut_number(expression_protCoding$SD, 10)
expression$signal_to_noise_bin <- cut_number(expression$signal_to_nois, 10)
expression_protCoding$signal_to_noise_bin <- cut_number(expression_protCoding$signal_to_nois, 10)



####
get_gene_per_bin <- function(expression, bin_type){
  # per expression bin calculate the proportion of cis/trans eQTLs and get info such as pLI and length per bin
  gene_per_bin <- data.frame()
  gene_lengths_per_bin <- data.frame()
  for(bin in sort(unique(expression[[bin_type]]))){
    expression_current_bin <- expression[expression[[bin_type]]==bin,]
    genes_current_bin <- rownames(expression_current_bin)
    n_genes <- length(genes_current_bin)
    n_cis_QTL <- sum(genes_current_bin %in% cis_qtl$ProbeName)

    mean_sd_cis <- round(mean(expression_current_bin[rownames(expression_current_bin) %in% cis_qtl$ProbeName,]$SD),2)
    mean_sd_not_cis <- round(mean(expression_current_bin[!rownames(expression_current_bin) %in% cis_qtl$ProbeName,]$SD),2)
  
    gene_lengths <- protein_coding_genes[match(rownames(expression_current_bin), protein_coding_genes$Ensembl.Gene.ID),]
    is_cis_qtl <- gene_lengths$Ensembl.Gene.ID %in% cis_qtl$ProbeName

    pLI$transcript <- gsub('\\.[0-9]+$','',pLI$transcript)
    pLI$gene <- ensemble_IDs[match(pLI$transcript, ensemble_IDs$`Transcript stable ID`),]$`Gene stable ID`
    df <- data.frame('bin'=bin, 'length'=gene_lengths$Transcript.length..including.UTRs.and.CDS.,
                   'is_cis_qtl'=is_cis_qtl,'gene'=gene_lengths$Ensembl.Gene.ID,
                   'expression'=expression_current_bin[as.character(gene_lengths$Ensembl.Gene.ID),]$mean_expression,
                   'sd'=expression_current_bin[as.character(gene_lengths$Ensembl.Gene.ID),]$SD,
                   'pLI'=pLI[match(gene_lengths$Ensembl.Gene.ID, pLI$gene),]$pLI)
    gene_lengths_per_bin <- rbind(gene_lengths_per_bin, df)
  
    df <- data.frame('n_genes'=c((n_cis_QTL/n_genes)*100, 
                               ((n_genes-n_cis_QTL)/n_genes)*100), 
                   'bin'=c(bin, bin),
                   'qtl_type'=c('cis','cis'),
                   'genes_are_QTL'=c('yes','no'))
    gene_per_bin <- rbind(gene_per_bin, df)
  }
  returnList <- list("gene_per_bin" = gene_per_bin, "gene_lengths_per_bin" = gene_lengths_per_bin)
  return(returnList)
}



bin_results_protCoding <- get_gene_per_bin(expression_protCoding,'expression_bin')
sd_bin_results_protCoding <- get_gene_per_bin(expression_protCoding,'sd_bin')
singal_to_noise_bin_results_protCoding<- get_gene_per_bin(expression_protCoding,'signal_to_noise_bin')

sd_bin_results <- get_gene_per_bin(expression,'sd_bin')
bin_results<- get_gene_per_bin(expression,'expression_bin')
singal_to_noise_bin_results <- get_gene_per_bin(expression,'signal_to_noise_bin')
#####
 
make_plots <- function(gene_per_bin, gene_lengths_per_bin, out_dir, xlab){
  #### make the histogram ####
  ggplot(gene_per_bin, aes(bin, n_genes, fill=genes_are_QTL))+
    geom_bar(stat='identity')+
    scale_fill_manual(values=c('grey95','lightblue'))+
    theme_pubr(base_size=18)+ 
    guides(fill=FALSE)+
    xlab(xlab)+
    ylab('Proportion of genes showing\ncis-eQTL effect')+
    scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                       labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))+
    theme(axis.text= element_text(colour='grey70'))+
    scale_x_discrete(breaks = c('[0.00778,0.774]','(10.9,19.3]'),
                     labels=c('low','high'))
  ggsave(paste0(out_dir,'proportion_of_QTL_per_bin_only_cis.png'),width=8, height=5)  
  
  ####
  #### same per sd bin ####
  ggplot(gene_per_bin, aes(bin, n_genes, fill=genes_are_QTL))+
    geom_bar(stat='identity')+
    scale_fill_manual(values=c('grey95','lightblue'))+
    theme_pubr(base_size=18)+ 
    guides(fill=FALSE)+
    xlab(xlab)+
    ylab('Proportion of genes showing\ncis-eQTL effect')+
    scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                       labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))+
    theme(axis.text= element_text(colour='grey70'))+
    scale_x_discrete(breaks = c('[0.00778,0.774]','(10.9,19.3]'),
                     labels=c('low','high'))
  ggsave(paste0(out_dir,'proportion_of_QTL_per_bin_only_cis.png'),width=8, height=5) 
  
  
  #### Make pLI plot ####
  gene_lengths_per_bin$is_cis_qtl <- factor(gene_lengths_per_bin$is_cis_qtl, levels=c(TRUE,FALSE))
  gene_lengths_per_bin_tmp <- gene_lengths_per_bin[!is.na(gene_lengths_per_bin$pLI),]
  ggplot(gene_lengths_per_bin_tmp, aes(is_cis_qtl,pLI, fill=is_cis_qtl))+
    geom_violin()+
    stat_summary(fun.y = median, geom = "point", shape=23, size = 3)+
    theme_bw(base_size=18)+
    facet_grid(~bin)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top')+
    scale_fill_manual(values=c('lightblue','grey95'))+
    xlab('')
  ggsave(paste0(out_dir,'pLI_per_bin.png'),height=5, width=20)
  
  ####
}

make_plots(bin_results_protCoding$gene_per_bin, bin_results_protCoding$gene_lengths_per_bin, paste0(opt$out_path,'proteinCoding/'),xlab='Average brain gene expression')
make_plots(sd_bin_results_protCoding$gene_per_bin, sd_bin_results_protCoding$gene_lengths_per_bin, paste0(opt$out_path,'proteinCoding_sd_bin/'),xlab='Brain gene expression standard deviation')
make_plots(sd_bin_results_protCoding$gene_per_bin, sd_bin_results_protCoding$gene_lengths_per_bin, paste0(opt$out_path,'proteinCoding_signalToNoise_bin/'),xlab='Signal to noise')

make_plots(bin_results$gene_per_bin, bin_results$gene_lengths_per_bin, paste0(opt$out_path,'all/'),xlab='Average brain gene expression')
make_plots(bin_results$gene_per_bin, bin_results$gene_lengths_per_bin, paste0(opt$out_path,'all_sd_bin/'),xlab='Brain gene expression standard deviation')
make_plots(bin_results$gene_per_bin, bin_results$gene_lengths_per_bin, paste0(opt$out_path,'all_signalToNoise_bin/'),xlab='Signal to noise')

