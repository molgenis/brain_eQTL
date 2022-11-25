library("optparse")
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(topGO)
library(ALL)
library(dplyr)
library(data.table)

# Get command line arguments 
option_list = list(
  make_option(c("-f", "--fdrQTLs"), type="character",
              help="DER-08a_hg38_eQTL.significant.txt of psychEncode", metavar="character"),
  make_option(c("-b", "--bonferonniQTLs"), type="character",
              help="DER-08b_hg38_eQTL.bonferroni.txt.gz of psychEncode", metavar="character"),
  make_option(c("-e", "--expression"), type="character",
              help="File with mean expression per gene", metavar="character"),
  make_option(c("-p", "--proteinCodingGenes"), type="character",
              help="File with protein coding genes", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#opt$fdrQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/DER-08a_hg38_eQTL.significant.txt"
#opt$bonferonniQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/DER-08b_hg38_eQTL.bonferroni.txt.gz"
#opt$expression <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/DER-02_PEC_Gene_expression_matrix_TPM.txt"
#opt$proteinCodingGenes <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/protein_coding_genes_ensembl84.txt"


##### read in data ####
fdr_qtl <- read.table(gzfile(opt$fdrQTLs), header=T, sep='\t')
bonferonni_qtl <- read.table(gzfile(opt$bonferonniQTLs), header=T, sep='\t')
protein_coding_genes <- read.table(opt$proteinCodingGenes, sep='\t', header=T)

# make sure all are FDR < 0.05, and select top snps
fdr_qtl <- fdr_qtl[fdr_qtl$FDR < 0.05,]
fdr_qtl <- fdr_qtl  %>% 
  group_by(gene_id) %>% 
  dplyr::slice(which.min(FDR))

bonferonni_qtl <- bonferonni_qtl[bonferonni_qtl$FDR < 0.05,]
bonferonni_qtl <- bonferonni_qtl  %>% 
  group_by(gene_id) %>% 
  dplyr::slice(which.min(FDR))


expression <- data.frame(fread(opt$expression, header=T, sep='\t'))
rownames(expression) <- expression$GeneName
expression$GeneName <- NULL

# genes still have version numbers, remove
fdr_qtl$gene_id <- sub('\\..*', '', fdr_qtl$gene_id)
bonferonni_qtl$gene_id <- sub('\\..*', '', bonferonni_qtl$gene_id)

# select only protein coding genes
expression_prot_coding <- expression[rownames(expression) %in% protein_coding_genes$Ensembl.Gene.ID,]
fdr_qtl <- fdr_qtl[fdr_qtl$gene_id %in% protein_coding_genes$Ensembl.Gene.ID,]
bonferonni_qtl <- bonferonni_qtl[bonferonni_qtl$gene_id %in% protein_coding_genes$Ensembl.Gene.ID,]
#####

##### per expression bin calculate the proportion of cis/trans eQTLs #####
mean_expression <- data.frame('mean_expression'=rowSums(expression_prot_coding))

mean_expression$mean_expression <- log(mean_expression$mean_expression+0.0000000000000000000001)
mean_expression$expression_bin <- cut_number(mean_expression$mean_expression, 10)

gene_per_bin <- data.frame()
gene_lengths_per_bin <- data.frame()
for(bin in sort(unique(mean_expression$expression_bin))){
  expression_current_bin <- mean_expression[mean_expression$expression_bin==bin,]
  genes_current_bin <- rownames(expression_current_bin)
  n_genes <- length(genes_current_bin)
  n_fdr_qtl <- sum(genes_current_bin %in% fdr_qtl$gene_id)
  n_bonferonni_qtl <- sum(genes_current_bin %in% bonferonni_qtl$gene_id)
  
  
  gene_lengths <- protein_coding_genes[match(rownames(expression_current_bin), protein_coding_genes$Ensembl.Gene.ID),]
  is_fdr_qtl <- gene_lengths$Ensembl.Gene.ID %in% fdr_qtl$gene_id
  is_bonferonni_qtl <- gene_lengths$Ensembl.Gene.ID %in% bonferonni_qtl$gene_id
  
  df <- data.frame('bin'=bin, 'length'=gene_lengths$Transcript.length..including.UTRs.and.CDS.,
                   'is_fdr_qtl'=is_fdr_qtl,'is_bonferonni_qtl'=is_bonferonni_qtl,'gene'=gene_lengths$Ensembl.Gene.ID,
                   'expression'=expression_current_bin[gene_lengths$Ensembl.Gene.ID,]$mean_expression)

  
  df <- data.frame('n_genes'=c((n_fdr_qtl/n_genes)*100, 
                               (n_bonferonni_qtl/n_genes)*100, 
                               ((n_genes-n_fdr_qtl)/n_genes)*100, 
                               ((n_genes-n_bonferonni_qtl)/n_genes)*100),
                   'bin'=c(bin, bin, bin, bin),
                   'qtl_type'=c('fdr','bonferonni','fdr','bonferonni'),
                   'genes_are_QTL'=c('yes','yes','no','no'))
  gene_per_bin <- rbind(gene_per_bin, df)
}
#####
 

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
ggsave('figures/psychEncode_proortion_of_QTL_per_bin_proteinCoding_only.png',width=8, height=5)  


ggplot(gene_per_bin[gene_per_bin$qtl_type=='bonferonni',], aes(bin, n_genes, fill=genes_are_QTL))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=c('grey95','lightblue'))+
  theme_pubr(base_size=18)+ 
  guides(fill=FALSE)+
  xlab('Average brain gene expression')+
  ylab('Proportion of genes showing\n cis-eQTL effect')+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))+
  theme(axis.text= element_text(colour='grey70'))+
  scale_x_discrete(breaks = c('[0.00778,0.774]','(10.9,19.3]'),
                   labels=c('low','high'))
ggsave('figures/psychEncode_proortion_of_QTL_per_bin_proteinCoding_only_bonferonni.png',width=8, height=5)  

####



