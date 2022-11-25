library("optparse")
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(plyr)
library(data.table)
library(stringr)
# Get command line arguments 
option_list = list(
  make_option(c("-c", "--cisQTLs"), type="character",
              help="eQTLProbesFDR0.05-ProbeLevel.txt.gz of cis-eQTLs", metavar="character"),
  make_option(c("-t", "--transQTLs"), type="character",
              help="eQTLProbesFDR0.05-ProbeLevel.txt.gz of trans-eQTLs", metavar="character"),
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
  make_option(c("-b", "--biosExpression"), type="character",
              help="Path to file expression data from BIOS", metavar="character"),
  make_option(c("-z", "--eqtlGenBins"), type="character",
              help="Path to file with binned data from eqtlGen", metavar="character"),
  make_option(c("-g", "--gtf"), type="character",
              help="GTF file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- list()
opt$cisQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/eQTLProbesFDR0.05-ProbeLevel.CIS.txt.gz"
opt$transQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/eQTLProbesFDR0.05-ProbeLevel.TRANS.txt.gz"
opt$expression <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.mean_and_sd.txt"
opt$proteinCodingGenes <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/protein_coding_genes_ensembl84.txt"
opt$pLI <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"
opt$out_path <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/figures/freeze2/"
opt$ensemble_IDs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/ensembl_transcript_gene_IDs.txt"
opt$biosExpression <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/gene_read_counts_BIOS_and_LLD_passQC.tsv.SampleSelection.ProbesWithZeroVarianceRemoved.TMM.txt.gz"
opt$eqtlGenBins <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/BIOS_expression_summary_bins_20180515.txt"
opt$gtf <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/gencode.v24.chr_patch_hapl_scaff.annotation.genesOnly.gtf"
##### read in data ####
cis_qtl <- read.table(gzfile(opt$cisQTLs), header=T, sep='\t')
trans_qtl <- read.table(gzfile(opt$transQTLs), header=T, sep='\t')
protein_coding_genes <- read.table(opt$proteinCodingGenes, sep='\t', header=F)
pLI <- read.table(opt$pLI, sep='\t', header=T)
ensemble_IDs <- fread(opt$ensemble_IDs)
BIOS_expression <- data.frame(fread(opt$biosExpression))
expression <- read.table(opt$expression, header=T, sep='\t', row.names=1)
eqtlGen_bins <- fread(opt$eqtlGenBins )
gtf <- fread(opt$gtf, header = F)
gtf$gene_id <- str_match(gtf$V9, 'gene_id "(ENSG.+?)"')[,2]
# make sure all are FDR < 0.05
cis_qtl <- cis_qtl[cis_qtl$FDR < 0.05,]
trans_qtl <- trans_qtl[trans_qtl$FDR < 0.05,]

# genes still have version numbers, remove
rownames(expression) <- sapply(strsplit(rownames(expression),"\\."), `[`, 1)
cis_qtl$ProbeName <- sapply(strsplit(as.character(cis_qtl$ProbeName),"\\."), `[`, 1)
trans_qtl$ProbeName <- sapply(strsplit(as.character(trans_qtl$ProbeName),"\\."), `[`, 1)

# select only protein coding genes
gtf_only_autosomal <- gtf[gtf$V1 %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
                                        'chr9','chr10','chr11','chr12','chr13','chr14','chr15',
                                        'chr16','chr17','chr18','chr19','chr20','chr21','chr22')]
gtf_only_autosomal$gene_id <- gsub('\\.[0-9]+','', gtf_only_autosomal$gene_id)
expression_protCoding <- expression[rownames(expression) %in% protein_coding_genes$V1,]
cis_qtl_protCoding <- cis_qtl[cis_qtl$ProbeName %in% protein_coding_genes$V1,]
trans_qtl_protCoding <- trans_qtl[trans_qtl$ProbeName %in% protein_coding_genes$V1,]
expression_protCoding_autosomal <- expression_protCoding[rownames(expression_protCoding) %in% gtf_only_autosomal$gene_id,]

# select eQTLgen genes
expression_eqtlGen_genes <- expression[rownames(expression) %in% eqtlGen_bins$ENSG,]
#####


##### make expression bins #####

expression$expression_bin <- cut_number(expression$mean_expression, 10)
expression_protCoding$expression_bin <- cut_number(expression_protCoding$mean_expression, 10)
expression_protCoding$expression_bin_order <- mapvalues(expression_protCoding$expression_bin, 
                                                        c("[4.75e-05,0.666]","(0.666,2.37]","(2.37,4.61]","(4.61,6.45]","(6.45,7.61]","(7.61,8.44]",
                                                          "(8.44,9.15]","(9.15,9.84]","(9.84,10.7]","(10.7,18.2]"),
                                                        c(1, 2, 3,4,5,6,7,8,9,10))
expression_eqtlGen_genes$expression_bin <- cut_number(expression_eqtlGen_genes$mean_expression, 10)
expression_eqtlGen_genes$expression_bin_order <- mapvalues(expression_eqtlGen_genes$expression_bin, 
                                                           c("[8.35e-05,1.53]","(1.53,2.92]","(2.92,4.47]","(4.47,6.11]","(6.11,7.39]","(7.39,8.29]",
                                                             "(8.29,9.04]","(9.04,9.75]","(9.75,10.6]","(10.6,18.2]"),
                                                           c(1, 2, 3,4,5,6,7,8,9,10))
expression_eqtlGen_genes$expression_bin_order_eqtlgen <- eqtlGen_bins[match(rownames(expression_eqtlGen_genes), eqtlGen_bins$ENSG),]$cis_exp_mean_bin
expression_eqtlGen_genes$expression_bin_order <- gsub('bin','',expression_eqtlGen_genes$expression_bin_order_eqtlgen)


expression_protCoding_autosomal$expression_bin <- cut_number(expression_protCoding_autosomal$mean_expression, 10)
expression_protCoding_autosomal$expression_bin_order <- mapvalues(expression_protCoding_autosomal$expression_bin, 
                                                           c("[4.75e-05,0.77]","(0.77,2.5]","(2.5,4.7]","(4.7,6.51]","(6.51,7.65]","(7.65,8.47]","(8.47,9.17]",
                                                             "(9.17,9.85]","(9.85,10.7]","(10.7,16]"),
                                                           c(1, 2, 3,4,5,6,7,8,9,10))
#### For BIOS, calcualte mean expression and SD per gene ####
rownames(BIOS_expression) <- BIOS_expression$V1
BIOS_expression$V1 <- NULL
BIOS_expression_matrix <- BIOS_expression
BIOS_expression_matrix[BIOS_expression_matrix < 0] <- 0
BIOS_expression_matrix <- data.frame(BIOS_expression_matrix)
rownames(BIOS_expression_matrix) <- rownames(BIOS_expression)
BIOS_expression_matrix_log2 <- log2(BIOS_expression_matrix+1)


BIOS_expression_mean <- rowMeans(BIOS_expression_matrix_log2)
BIOS_expression_sd <- rowSds(as.matrix(BIOS_expression_matrix_log2)) # same result as before.
BIOS_expression_mean_and_sd <- data.frame(BIOS_expression_mean, BIOS_expression_sd)
colnames(BIOS_expression_mean_and_sd) <- c('mean_expression','SD')
BIOS_expression_mean_and_sd$expression_bin <- cut_number(BIOS_expression_mean_and_sd$mean_expression, 10)
BIOS_expression_mean_and_sd$expression_bin_order = mapvalues(BIOS_expression_mean_and_sd$expression_bin, 
                                                             c("[0.000114,0.00183]","(0.00183,0.00503]","(0.00503,0.0157]","(0.0157,0.0559]","(0.0559,0.197]",
                                                               "(0.197,0.717]","(0.717,2.4]","(2.4,5.74]","(5.74,8.59]","(8.59,17.2]"),
                                                             c(1, 2, 3,4,5,6,7,8,9,10))
####
get_gene_per_bin <- function(expression){
  # per expression bin calculate the proportion of cis/trans eQTLs and get info such as pLI and length per bin
  gene_per_bin <- data.frame()
  gene_lengths_per_bin <- data.frame()
  for(bin in sort(unique(expression$expression_bin_order))){
    expression_current_bin <- expression[expression$expression_bin_order==bin,]
    genes_current_bin <- rownames(expression_current_bin)
    n_genes <- length(genes_current_bin)
    n_cis_QTL <- sum(genes_current_bin %in% cis_qtl$ProbeName)
    n_trans_QTL <- sum(genes_current_bin %in% trans_qtl$ProbeName)
  
    mean_sd_cis <- round(mean(expression_current_bin[rownames(expression_current_bin) %in% cis_qtl$ProbeName,]$SD),2)
    mean_sd_not_cis <- round(mean(expression_current_bin[!rownames(expression_current_bin) %in% cis_qtl$ProbeName,]$SD),2)
    
    
    gene_lengths <- protein_coding_genes[match(rownames(expression_current_bin), protein_coding_genes$V1),]
    gene_lengths <- data.frame(V1=gene_lengths, gene_length=0)
    is_cis_qtl <- gene_lengths$V1 %in% cis_qtl$ProbeName
    is_trans_qtl <- gene_lengths$V1 %in% trans_qtl$ProbeName
    
    pLI$transcript <- gsub('\\.[0-9]+$','',pLI$transcript)
    pLI$gene <- ensemble_IDs[match(pLI$transcript, ensemble_IDs$`Transcript stable ID`),]$`Gene stable ID`
    df <- data.frame('bin'=expression_current_bin$expression_bin, 
                     'bin_number'=expression_current_bin$expression_bin_order,
                     'length'=gene_lengths$gene_length,
                     'is_cis_qtl'=is_cis_qtl,
                     'is_trans_qtl'=is_trans_qtl,
                     'gene'=gene_lengths$V1,
                     'expression'=expression_current_bin[gene_lengths$V1,]$mean_expression,
                     'sd'=expression_current_bin[gene_lengths$V1,]$mean_expression,
                     'pLI'=pLI[match(gene_lengths$V1, pLI$gene),]$pLI)
    gene_lengths_per_bin <- rbind(gene_lengths_per_bin, df)
    
    df <- data.frame('n_genes'=c((n_cis_QTL/n_genes)*100, 
                               (n_trans_QTL/n_genes)*100, 
                               ((n_genes-n_cis_QTL)/n_genes)*100, 
                               ((n_genes-n_trans_QTL)/n_genes)*100),
                   'bin'=c(bin, bin, bin, bin),
                   'qtl_type'=c('cis','trans','cis','trans'),
                   'genes_are_QTL'=c('yes','yes','no','no'))
    gene_per_bin <- rbind(gene_per_bin, df)
  }
  returnList <- list("gene_per_bin" = gene_per_bin, "gene_lengths_per_bin" = gene_lengths_per_bin)
  return(returnList)
}

bin_results_protCoding<- get_gene_per_bin(expression_protCoding)
bin_results<- get_gene_per_bin(expression)
bin_results_eqtlGenGenes <- get_gene_per_bin(expression_eqtlGen_genes)
bin_results_protCoding_autosomal <- get_gene_per_bin(expression_protCoding_autosomal)
#####
 
make_plots <- function(gene_per_bin, gene_lengths_per_bin, out_dir){
  #### make the histogram ####
  ggplot(gene_per_bin, aes(bin, n_genes, fill=genes_are_QTL))+
    geom_bar(stat='identity')+
    scale_fill_manual(values=c('grey95','lightblue'))+
    theme_pubr(base_size=18)+ 
    guides(fill=FALSE)+
    facet_wrap(~qtl_type)+ 
    xlab('Expression levels')+
    ylab('Proprotion QTLs')+
    scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                       labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))+
    theme(axis.text= element_text(colour='grey70'))+
    scale_x_discrete(breaks = c('[0.00778,0.774]','(10.9,19.3]'),
                   labels=c('low','high'))
  ggsave(paste0(out_dir,'proportion_of_QTL_per_bin.png'),width=8, height=5)  
  
  
  gene_per_bin$bin <- factor(gene_per_bin$bin, levels=c('not_tested_in_cis','1','2','3','4','5','6','7','8','9','10'))
  ggplot(gene_per_bin[gene_per_bin$qtl_type=='cis',], aes(bin, n_genes, fill=genes_are_QTL))+
        geom_bar(stat='identity')+
        scale_fill_manual(values=c('grey95','lightblue'))+
        theme_pubr(base_size=18)+ 
        guides(fill=FALSE)+
        xlab('Average brain gene expression')+
        ylab('Proportion of genes showing cis-eQTL effect')+
        scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                         labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))+
        theme(axis.text= element_text(colour='grey70'),
              axis.text.x = element_text(angle = 90))#+))+

    #scale_x_discrete(breaks = c('[0.00778,0.774]','(10.9,19.3]'),
    #                 labels=c('low','high'))
  ggsave(paste0(out_dir,'proportion_of_QTL_per_bin_only_cis.png'),width=8, height=5)  
  
  ####

  #### per bin plot the length of the genes for the QTL genes and non-QTL genes ####
  p1 <- ggplot(gene_lengths_per_bin, aes(log(length), fill=is_cis_qtl))+
    geom_histogram(position="identity", alpha=0.5)+
    facet_wrap(~bin, ncol=2)+
    theme_pubr(base_size=18)+
    scale_fill_manual(values=c('grey70','blue'))+
    labs(fill="is cis QTL")
  ggsave(paste0(out_dir,'length_per_expression_bin_cis.pdf'),width=8, height=12)  
  
  
  p2 <- ggplot(gene_lengths_per_bin, aes(log(length), fill=is_trans_qtl))+
    geom_histogram(position="identity", alpha=0.5)+
    facet_wrap(~bin, ncol=2)+
    theme_pubr(base_size=18)+
    scale_fill_manual(values=c('grey70','blue'))+
    labs(fill="is trans QTL")
  ggsave(paste0(out_dir,'length_per_expression_bin_trans.pdf'),width=8, height=12)  
  
  pdf(paste0(out_dir,'length_per_expression_bin.pdf'),width=12, height=12)  
  grid.arrange(p1, p2, nrow = 1)
  dev.off()
  ####
  
  
  #### writing genes from highest bin for GO analysis ####
  highest_bin <- gene_lengths_per_bin[gene_lengths_per_bin$bin == "(10.9,19.3]",]
  write.table(as.character(highest_bin[highest_bin$is_cis_qtl,]$gene), 'cis_qtl_genes_highest_bin.txt',
              quote=F, row.names=F)
  write.table(as.character(highest_bin[!highest_bin$is_cis_qtl,]$gene), 'not_cis_qtl_genes_highest_bin.txt',
              quote=F, row.names=F)
  ####
  
  #### plot sd per bin for cis vs non cis ####
  p1 <- ggplot(gene_lengths_per_bin, aes(sd, fill=is_cis_qtl))+
    geom_histogram(position='identity')+
    theme_pubr(base_size=18)+ 
    scale_fill_manual(values=c('grey95','lightblue'))+
    facet_wrap(~bin, ncol=2,scale='free_x')
  
  p2 <- ggplot(gene_lengths_per_bin, aes(sd, fill=is_trans_qtl))+
    geom_histogram(position='identity')+
    theme_pubr(base_size=18)+ 
    scale_fill_manual(values=c('grey95','lightblue'))+
    facet_wrap(~bin, ncol=2,scale='free_x')
  pdf(paste0(out_dir,'sd_per_expression_bin.pdf'),width=12, height=12)  
  grid.arrange(p1, p2, nrow = 1)
  dev.off()
  ####
  
  #### Make pLI plot ####
  gene_lengths_per_bin$is_cis_qtl <- factor(gene_lengths_per_bin$is_cis_qtl, levels=c(TRUE,FALSE))
  ggplot(gene_lengths_per_bin, aes(is_cis_qtl,pLI, fill=is_cis_qtl))+
    geom_violin()+
    stat_summary(fun.y = median, geom = "point", shape=23, size = 3)+
    theme_bw(base_size=18)+
    facet_grid(~bin)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top')+
    scale_fill_manual(values=c('lightblue','grey95'))+
    xlab('')
  ggsave(paste0(out_dir,'pLI_per_bin.png'),height=5, width=12)
  
  ####
}

make_plots(bin_results_protCoding$gene_per_bin, bin_results_protCoding$gene_lengths_per_bin, paste0(opt$out_path,'proteinCoding/'))
make_plots(bin_results$gene_per_bin, bin_results$gene_lengths_per_bin, paste0(opt$out_path,'all/'))
make_plots(bin_results_eqtlGenGenes$gene_per_bin, bin_results_eqtlGenGenes$gene_lengths_per_bin, paste0(opt$out_path,'eqtlGen_genes/'))
make_plots(bin_results_protCoding_autosomal$gene_per_bin, bin_results_protCoding_autosomal$gene_lengths_per_bin, paste0(opt$out_path,'proteinCoding_autosomal/'))

expression_protCoding$name <- 'MetaBrain'
BIOS_expression_mean_and_sd$name <- 'BIOS'
expression_metaBrain_bios <- rbind(expression_protCoding,BIOS_expression_mean_and_sd)
ggplot(expression_metaBrain_bios, aes(expression_bin_order, SD,fill=name))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  xlab('Expression bin')
ggsave(paste0(opt$out_path,'proteinCoding/sd_comparison.png'),width=8, height=8)
ggplot(expression_metaBrain_bios, aes(expression_bin_order, mean_expression,fill=name))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  xlab('Expression bin')
ggsave(paste0(opt$out_path, 'proteinCoding/mean_comparison.png'),width=8, height=8)


write.table(bin_results_eqtlGenGenes$gene_lengths_per_bin, file='2019-06-05-metabrain-eQTLgenGenes-bins.txt',sep='\t',quote=F,row.names=F)

