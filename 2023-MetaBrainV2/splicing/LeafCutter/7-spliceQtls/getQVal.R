library('qvalue')

args = commandArgs(trailingOnly=TRUE)

input<-args[1]
output<-args[2]

fdr=0.05

data<-read.table(input,header=T,fill=T)
# using lambda of 0.85. Why?
# drop qval if it's already there
cols<-names(data) 
cols=cols[cols != "Qval"] 
cols=cols[cols != "qval"] 
data <- data[, cols]
cols<-names(data) 
# print(cols)

data <- data[!is.na(data$BetaAdjustedMetaP),]

data$qval<-qvalue(data$BetaAdjustedMetaP,lambda=0.85)$qvalues 

# got this from: https://github.com/francois-a/fastqtl/blob/master/R/calculateSignificanceFastQTL.R
# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(data[data$qval > fdr, 'BetaAdjustedMetaP'])[1]  # smallest p-value above FDR
print(ub)
lb <- -sort(-data[data$qval <= fdr, 'BetaAdjustedMetaP'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("  * min p-value threshold @ FDR ", fdr, ": ", pthreshold, "\n", sep="")
data$PvalueNominalThreshold <- signif(qbeta(pthreshold, data[, 'BetaDistAlpha'], data[, 'BetaDistBeta'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

# reorder some cols
cols<-names(data)
#print(cols)
cols=cols[cols != "qval"]
cols=cols[cols != "PvalueNominalThreshold"]
cols<-c(cols,"PvalueNominalThreshold")
cols<-c(cols,"qval")

data <- data[, cols]

data <- data[order(data$qval),]
# print(head(data))

write.table(data, file=output,sep='\t', quote=F, row.names=F)
