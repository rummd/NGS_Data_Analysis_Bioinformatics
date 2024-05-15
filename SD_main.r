if(!require("BiocManager",quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("pheatmap")

###setting directory
#setwd("D:/roject")

###DESeq2 Preliminary Analysis
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("tidyverse")

Counts=read.csv("GSE215456_featureCounts_globin.csv")
ColData=read.csv("New_file1.csv")



row.names(Counts) <- Counts[,1]
Counts <- Counts[,-1]

row.names(ColData) <- ColData[,1]

ColData$Design_1 <- factor(ColData$Design_1)

dds <- DESeqDataSetFromMatrix(countData = Counts, colData = ColData, design = ~Design_1)
dds

prdds <- DESeq(dds)
prdds

res <- results(prdds)
res

write.csv(res, "SD_output.csv")

res_new <- results(prdds, alpha = 0.1) 
res_new

write.csv(res_new, "Final_DF_expressed_SD.csv")

###Dispersion Plot
plotDispEsts(prdds, main="Dispersion Estimates")

### MA plot 

plotMA(res_new, ylim=c(-5,5), main = "Treatment vs Diagnosis")

###Volcano Plot

smoc2_res_all <- data.frame(res_new) %>% mutate(threshold = padj < 0.05)

# Create the volcano plot 
ggplot(smoc2_res_all,aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  geom_point() + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))

###Heatmap 

# Determine the size factors to use for normalization
dds_smoc2 <- estimateSizeFactors(dds)

# Extract the normalized counts
smoc2_normalized_counts <- counts(dds_smoc2, normalized=TRUE)

# Transform the normalized counts 
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Extract the matrix of transformed counts
vsd_mat_smoc2 <- assay(vsd_smoc2)

# Compute the correlation values between samples
vsd_cor_smoc2 <- cor(vsd_mat_smoc2) 

# Plot the heatmap
pheatmap(vsd_cor_smoc2, annotation = select(ColData, Design_1))


