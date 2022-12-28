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

setwd("~/Documents/Uni_Work/Introduction_to_Bioinformatics/Project/Phase_2")
series <- "GSE48558"
platform <- "GPL6244"
pval_threshold <- 0.05

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

### Only take the relevant group of samples for analysis
indices_to_remove <- sort(match(c("GSM1180912", "GSM1180893", "GSM1180889",
                                  "GSM1180887", "GSM1180892", "GSM1180888",
                                  "GSM1180891", "GSM1180890", "GSM1180826",
                                  "GSM1180794", "GSM1180824", "GSM1180790",
                                  "GSM1180838", "GSM1180834", "GSM1180839",
                                  "GSM1180831", "GSM1180847", "GSM1180843",
                                  "GSM1180820", "GSM1180845", "GSM1180841",
                                  "GSM1180818", "GSM1180829"),
                                gset$geo_accession))
ex <- ex[,-indices_to_remove]
gr <- gr[-indices_to_remove]
gset <- gset[,-indices_to_remove]
pdf("Results/Corr_Heatmap_Relevant_Samples.pdf", width=15, height=15)
pheatmap(cor(ex), labels_row=gr)
dev.off()

### Differential Expression Analysis
gset$description <- gr
fac_gr <- factor(gr)
gset$description <- fac_gr
design <- model.matrix(~description + 0, gset)
colnames(design) <- levels(fac_gr)

fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(AML - Healthy, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
results <- topTable(fit2, adjust="fdr", sort.by='B', number=Inf)

results <- subset(results, select=c("Gene.symbol", "Gene.ID", "adj.P.Val", "logFC"))
write.table(results, file="Results/Differential_Expression_AML_vs_Healthy.txt", row.names=FALSE, sep='\t', quote=FALSE)

### Save the names of the significantly differentiated genes
aml.up <- subset(results, logFC > 1 & adj.P.Val < pval_threshold)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file="Results/AML_vs_Healthy_Up_Genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
aml.down <- subset(results, logFC < -1 & adj.P.Val < pval_threshold)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file="Results/AML_vs_Healthy_Down_Genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
