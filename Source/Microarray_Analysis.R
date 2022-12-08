# Install Deps: install.packages(c("GEOquery", "limma", "pheatmap", "ggplot2", "reshape2", "plyr"))
#               install.packages("Rtsne")

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(plyr)
library(MASS)
library(Rtsne)

setwd("~/Documents/Uni_Work/Introduction_to_Bioinformatics/Project/Phase_1")
series <- "GSE48558"
platform <- "GPL6244"

gset <- getGEO("GSE48558", GSEMatrix=TRUE, AnnotGPL=TRUE, destdir="Data/")
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gr_raw_string <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
                        "XXXXXXXXXXXXXXXXXX1X1XXX1X1111X1XX11XX11X1X1X1X1X1",
                        "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
                        "11111111111111111111")
gr_raw_string <- strsplit(gr_raw_string, "")[[1]]
filter <- sort(union(which(gr_raw_string %in% "0"),
                     which(gr_raw_string %in% "1")))
gr_raw_string <- gr_raw_string[filter]
class_mapper = c('0'="AML", '1'="Healthy")
gr <- as.character(class_mapper[as.factor(gr_raw_string)])
gset <- gset[,filter]

### Quality Control Phase
ex <- exprs(gset)          # Result is already in a logarithmic scale, as checked with min and max.

max(ex)
min(ex)

pdf("Results/Boxplot.pdf", width=64)
boxplot(ex)
dev.off()                  # Data is already normalized.

pdf("Results/Gene_Means.pdf", width=64)
barplot(rowMeans(ex))
dev.off()                  # Should zero out all of the means for PCA

ex.scale <- t(scale(t(ex), scale=FALSE))
pdf("Results/Gene_Means_Zeroed.pdf", width=64)
barplot(rowMeans(ex.scale))
dev.off()

pdf("Results/Corr_Heatmap.pdf", width=15, height=15)
pheatmap(cor(ex), labels_row=gr)
dev.off()
ind_to_remove <- match("GSM1180835", gset$geo_accession)
ex <- ex[,-ind_to_remove]
gr <- gr[-ind_to_remove]
gset <- gset[,-ind_to_remove]
pdf("Results/Corr_Heatmap_No_Outliers.pdf", width=15, height=15)
pheatmap(cor(ex), labels_row=gr, labels_col=gset$source_name_ch1)
dev.off()
ex.scale <- t(scale(t(ex), scale=FALSE))

### Dimension Reduction
###### PCA
pc <- prcomp(ex)
pdf("Results/PCA/PCA_Genes_Not_Zeroed.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()
pc <- prcomp(ex.scale)
pdf("Results/PCA/PCA_Genes.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()
pc_samples <- data.frame(pc$r[,1:3], Group=gr)
pdf("Results/PCA/PCA_Samples.pdf")
ggplot(pc_samples, aes(PC1, PC2, color=Group)) + geom_point(size=3)
dev.off()
###### MDS
distance_matrix <- dist(t(ex.scale))
fit <- cmdscale(distance_matrix, eig=TRUE, k=2)
mds_samples <- data.frame(fit$points, Group=gr)
pdf("Results/MDS/Metric_MDS_Samples.pdf")
ggplot(mds_samples, aes(X1, X2, color=Group)) + geom_point(size=3)
dev.off()
fit <- isoMDS(distance_matrix, k=2)
mds_samples <- data.frame(fit$points, Group=gr)
pdf("Results/MDS/Non-Metric_MDS_Samples.pdf")
ggplot(mds_samples, aes(X1, X2, color=Group)) + geom_point(size=3)
dev.off()
###### t-SNE
tsne_samples <- Rtsne(t(ex.scale), perplexity=15, check_duplicates=FALSE)
tsne_samples <- data.frame(tsne_samples$Y, Group=gr)
pdf("Results/t-SNE/t-SNE_Samples.pdf")
ggplot(tsne_samples, aes(X1, X2, color=Group)) + geom_point(size=3)
dev.off()

### Group Correlation: This is basically already done in the above pheatmap, after the removal of the outlier
