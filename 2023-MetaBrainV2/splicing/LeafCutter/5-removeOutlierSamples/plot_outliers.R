library(ggplot2)
setwd("/groups/umcg-biogen/tmp01/output/2021-FreezeThree/2021-02-18-splicing/2022-09-22-all-samples-leafcutter/")

basalganglia <- read.table("7-removeOutlierSamples/output/Basalganglia-EUR/pca_PCs.txt", header = T)
cerebellum <- read.table("7-removeOutlierSamples/output/Cerebellum-EUR/pca_PCs.txt", header = T)
cortexEur <- read.table("7-removeOutlierSamples/output/Cortex-EUR/pca_PCs.txt", header = T)
cortexAfr <- read.table("7-removeOutlierSamples/output/Cortex-AFR/pca_PCs.txt", header = T)
cortexEas <- read.table("7-removeOutlierSamples/output/Cortex-EAS/pca_PCs.txt", header = T)
hippocampus <- read.table("7-removeOutlierSamples/output/Hippocampus-EUR/pca_PCs.txt", header = T)
spinalcord <- read.table("7-removeOutlierSamples/output/Spinalcord-EUR/pca_PCs.txt", header = T)

basalganglia_scaled <- data.frame(apply(basalganglia, 2, scale))
ggplot(basalganglia_scaled, aes(PC1,PC2)) + geom_point() + 
  geom_vline(xintercept=3, linetype="dashed", color = "red") + 
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22)) +
  ggtitle("Basalganglia-EUR")

cerebellum_scaled <- data.frame(apply(cerebellum, 2, scale))
ggplot(cerebellum_scaled, aes(PC1,PC2)) + geom_point() + 
  geom_vline(xintercept=3, linetype="dashed", color = "red") + 
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22)) +
  ggtitle("Cerebellum-EUR")

cortexEur_scaled <- data.frame(apply(cortexEur, 2, scale))
ggplot(cortexEur_scaled, aes(PC1,PC2)) + geom_point() + 
  geom_vline(xintercept=3, linetype="dashed", color = "red") + 
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22)) +
  ggtitle("Cortex-EUR")

cortexAfr_scaled <- data.frame(apply(cortexAfr, 2, scale))
ggplot(cortexAfr_scaled, aes(PC1,PC2)) + geom_point() + 
  geom_vline(xintercept=3, linetype="dashed", color = "red") + 
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22)) +
  ggtitle("Cortex-AFR")

hippocampus_scaled <- data.frame(apply(hippocampus, 2, scale))
ggplot(hippocampus_scaled, aes(PC1,PC2)) + geom_point() + 
  geom_vline(xintercept=3, linetype="dashed", color = "red") + 
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22)) +
  ggtitle("Hippocampus-EUR")

spinalcord_scaled <- data.frame(apply(spinalcord, 2, scale))
ggplot(hippocampus_scaled, aes(PC1,PC2)) + geom_point() + 
  geom_vline(xintercept=3, linetype="dashed", color = "red") + 
  geom_vline(xintercept=-3, linetype="dashed", color = "red") +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_hline(yintercept=-3, linetype="dashed", color = "red") +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size=22)) +
  ggtitle("Spinalcord-EUR")
