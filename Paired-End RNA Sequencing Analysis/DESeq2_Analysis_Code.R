#Check working directory is
#/path/to/data/
getwd()

#########################################
#Load libraries
#########################################

#Data manipulation
library(magrittr)
library(tidyverse)
library(dplyr)

#Perform analysis
library(DESeq2)

#Visualisation
library(apeglm)
library(ggplot2)
library(tinytex)
library(limma)
library(ggrepel)

#Gene name annotation
library(biomaRt)

#Gene set enrichment analysis
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(org.Hs.eg.db)

#########################################
#Data manipulation
#########################################

#Import count data
count.data <- read.table(file="counts.txt", sep="\t", header=TRUE, check.names=F, row.names=1)

#Remove unwanted columns
cleaned_counts = subset(count.data, select = -c(Chr, Start, End, Strand, Length))

#Check dataframe dimensions
dim(cleaned_counts)
head(cleaned_counts)

#Change column/sample names to a simpler term (easier to read)
colnames(cleaned_counts) <- gsub("/mnt/scratch/umnpo/star_alignments/", "", colnames(cleaned_counts))
colnames(cleaned_counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(cleaned_counts))

#########################################
#Create metadata
#########################################

#This data has 4 conditions (PBS, 6h, 24h, and 48h), 5 repeats for each sample (e.g. A1, A2, A3, A4, and A5)
#And lot and passage information for the cells. Lot and passage could be confounding factors so these will be controlled for
colData <- data.frame(
  sample = c("A1","B1","C1","D1",
             "A2","B2","C2","D2",
             "A3","B3","C3","D3",
             "A4","B4","C4","D4",
             "A5","B5","C5","D5"),
  condition = rep(c("PBS", "6h-rBP-1", "24h-rBP-1", "48h-rBP-1"), times = 5),
  lot = c(rep("025", 4), rep("028", 4), rep("019", 4), rep("025", 4), rep("028", 4)),
  passage = c(rep("P3", 12), rep("P4", 8))
)
rownames(colData) <- colData$sample

#Check
rownames(colData)

#Create a name map to make colnames rownmames 
new_names <- rownames(colData)
colnames(cleaned_counts) <- new_names

#Save file
write.csv(as.data.frame(keep_counts),file="keep_counts_test2.csv")

#########################################
#Perform comparisons and DESeq2 analysis
#########################################

#PBS as reference
colData$condition <- factor(colData$condition, levels = c("PBS", "6h-rBP-1", "24h-rBP-1", "48h-rBP-1"))
#Passage and lot stated before condition to control for confounders first
dds <- DESeqDataSetFromMatrix(countData = cleaned_counts, colData = colData, design = ~ passage + lot + condition)
#Prioritise results with at least 10 or more reads in 5 samples (this ensures DEG is present in all repeats)
#Thus higher confidence DEGs
keep <-rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]
dds <- DESeq(dds)
summary(dds)

#DESeq2 automatically performs these, however they can be applied just in case
ddsSF <- estimateSizeFactors(dds)
ddsED <- estimateDispersions(ddsSF)
ddsnb <- nbinomWaldTest(ddsED, maxit = 1000) 

#State comparisons with PBS as reference
#("Condition", "Treatment", "Control")
res_pbs_vs_6h  <- results(ddsnb, contrast = c("condition", "6h-rBP-1", "PBS"))
summary_res1 <-summary(res_pbs_vs_6h)
res_pbs_vs_24h <- results(ddsnb, contrast = c("condition", "24h-rBP-1", "PBS"))
summary_res2 <- summary(res_pbs_vs_24h)
res_pbs_vs_48h <- results(ddsnb, contrast = c("condition", "48h-rBP-1", "PBS"))
summary_res3 <-summary(res_pbs_vs_48h)

#Identify coefficients for MA plots
resultsNames(ddsnb)

#########################################
#Shrink data
#########################################

#Plot MA plots for each comparison to check if data needs shrinking
plotMA(res_pbs_vs_6h, ylim = c(-1.5, 1.5), main="PBS vs 6h (before shrinkage)")
res_pbs_vs_6h_shrunk <- lfcShrink(ddsnb, coef="condition_6h.rBP.1_vs_PBS", type="apeglm")
plotMA(res_pbs_vs_6h_shrunk, ylim=c(-1.5,1.5), main="PBS vs 6h (after shrinkage)")

plotMA(res_pbs_vs_24h, ylim = c(-1.5, 1.5), main="PBS vs 24h (before shrinkage)")
res_pbs_vs_24h_shrunk <- lfcShrink(ddsnb, coef="condition_24h.rBP.1_vs_PBS", type="apeglm")
plotMA(res_pbs_vs_24h_shrunk, ylim=c(-1.5,1.5), main="PBS vs 24h (after shrinkage)")

plotMA(res_pbs_vs_48h, ylim = c(-1.5, 1.5), main="PBS vs 48h (before shrinkage)")
res_pbs_vs_48h_shrunk <- lfcShrink(ddsnb, coef="condition_48h.rBP.1_vs_PBS", type="apeglm")
plotMA(res_pbs_vs_48h_shrunk, ylim=c(-1.5,1.5), main="PBS vs 48h (after shrinkage)")

#########################################
#Identify significant DEGs
#########################################

#Pull out significant (padj < 0.05) differentially expressed genes (DEGs) from shrunk dataset 
resSig_pbs_6h <- res_pbs_vs_6h_shrunk[which(res_pbs_vs_6h_shrunk$padj < 0.05), ]
head(resSig_pbs_6h[order(resSig_pbs_6h$padj), ])
summary(resSig_pbs_6h)
resSig_pbs_6h.df <- as.data.frame(resSig_pbs_6h)
resSig_pbs_6h <- (resSig_pbs_6h.df[order(resSig_pbs_6h.df$padj), ])

#Save DEGs
write.csv(resSig_pbs_6h, "Sig_DEGs_pbs_6h_shrunk.csv", row.names = TRUE)

resSig_pbs_24h <- res_pbs_vs_24h_shrunk[which(res_pbs_vs_24h_shrunk$padj < 0.05), ]
head(resSig_pbs_24h[order(resSig_pbs_24h$padj), ])
summary(resSig_pbs_24h)
resSig_pbs_24h.df <- as.data.frame(resSig_pbs_24h)
resSig_pbs_24h <- (resSig_pbs_24h.df[order(resSig_pbs_24h.df$padj), ])

#Save DEGs
write.csv(resSig_pbs_24h, "Sig_DEGs_pbs_24h_shrunk.csv", row.names = TRUE)

resSig_pbs_48h <- res_pbs_vs_48h_shrunk[which(res_pbs_vs_48h_shrunk$padj < 0.05), ]
head(resSig_pbs_48h[order(resSig_pbs_48h$padj), ])
summary(resSig_pbs_48h)
resSig_pbs_48h.df <- as.data.frame(resSig_pbs_48h)
resSig_pbs_48h <- (resSig_pbs_48h.df[order(resSig_pbs_48h.df$padj), ])

#Save DEGs
write.csv(resSig_pbs_48h, "Sig_DEGs_pbs_48h_shrunk.csv", row.names = TRUE)

#########################################
#Perform comparisons for other factors
#########################################

#Change reference to 6h for other comparisons
#Repeat above steps
dds$condition <- relevel(dds$condition, ref = "6h-rBP-1")
ddsnb <- nbinomWaldTest(dds)
resultsNames(ddsnb)

res_6h_vs_24h <- results(ddsnb, name = "condition_24h.rBP.1_vs_6h.rBP.1")
summary_res4 <- summary(res_6h_vs_24h)

res_6h_vs_48h <- results(ddsnb, name = "condition_48h.rBP.1_vs_6h.rBP.1")
summary_res5 <- summary(res_6h_vs_48h)

plotMA(res_6h_vs_24h, ylim = c(-1.5, 1.5), main="6h vs 24h (before shrinkage)")
res_6h_vs_24h_shrunk <- lfcShrink(ddsnb, coef = "condition_24h.rBP.1_vs_6h.rBP.1",type = "apeglm")
plotMA(res_6h_vs_24h_shrunk, ylim=c(-1.5,1.5), main="6h vs 24h (after shrinkage)")

plotMA(res_6h_vs_48h, ylim = c(-1.5, 1.5), main="6h vs 48h (before shrinkage)")
res_6h_vs_48h_shrunk <- lfcShrink(ddsnb, coef = "condition_48h.rBP.1_vs_6h.rBP.1",type = "apeglm")
plotMA(res_6h_vs_48h_shrunk, ylim=c(-1.5,1.5), main="6h vs 48h (after shrinkage)")

resSig_6h_24h <- res_6h_vs_24h_shrunk[which(res_6h_vs_24h_shrunk$padj < 0.05), ]
head(resSig_6h_24h[order(resSig_6h_24h$padj), ])
summary(resSig_6h_24h)
resSig_6h_24h.df <- as.data.frame(resSig_6h_24h)
resSig_6h_24h <- (resSig_6h_24h.df[order(resSig_6h_24h.df$padj), ])
write.csv(resSig_6h_24h, "Sig_DEGs_6h_24h_shrunk.csv", row.names = TRUE)

resSig_6h_48h <- res_6h_vs_48h_shrunk[which(res_6h_vs_48h_shrunk$padj < 0.05), ]
head(resSig_6h_48h[order(resSig_6h_48h$padj), ])
summary(resSig_6h_48h)
resSig_6h_48h.df <- as.data.frame(resSig_6h_48h)
resSig_6h_48h <- (resSig_6h_48h.df[order(resSig_6h_48h.df$padj), ])
write.csv(resSig_6h_48h, "Sig_DEGs_6h_48h_shrunk.csv", row.names = TRUE)

#Change reference to 24h for other comparisons 
#Repeat above steps
dds$condition <- relevel(dds$condition, ref = "24h-rBP-1")
ddsnb <- nbinomWaldTest(dds)
resultsNames(ddsnb)

res_24h_vs_48h <- results(ddsnb, name = "condition_48h.rBP.1_vs_24h.rBP.1")
summary_res6 <- summary(res_24h_vs_48h)

plotMA(res_24h_vs_48h, ylim = c(-1.5, 1.5), main="24h vs 48h (before shrinkage)")
res_24h_vs_48h_shrunk <- lfcShrink(ddsnb, coef = "condition_48h.rBP.1_vs_24h.rBP.1",type = "apeglm")
plotMA(res_24h_vs_48h_shrunk, ylim=c(-1.5,1.5), main="24h vs 48h (after shrinkage)")

resSig_24h_48h <- res_24h_vs_48h_shrunk[which(res_24h_vs_48h_shrunk$padj < 0.05), ]
head(resSig_24h_48h[order(resSig_24h_48h$padj), ])
summary(resSig_24h_48h)
resSig_24h_48h.df <- as.data.frame(resSig_24h_48h)
resSig_24h_48h <- (resSig_24h_48h.df[order(resSig_24h_48h.df$padj), ])
write.csv(resSig_24h_48h, "Sig_DEGs_24h_48h_shrunk.csv", row.names = TRUE)

#########################################
#Data annotation
#########################################

#Annotate dataframes for gene names instead of ensembl gene IDs
#Create ensembl object
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Annotate all shrunk dataframes first
res_24h_vs_48h_shrunk_df <- as.data.frame(res_24h_vs_48h_shrunk)
res_24h_vs_48h_shrunk_df$ensembl_gene_id <- rownames(res_24h_vs_48h_shrunk_df)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res_24h_vs_48h_shrunk_df$ensembl_gene_id,
  mart = ensembl
)

res_24h_vs_48h_shrunk_df_annot <- merge(res_24h_vs_48h_shrunk_df, gene_map, by = "ensembl_gene_id", all.x = TRUE)

#Done for all shrunk dataframes but not shown

#Annotate all significant DEGs dataframes
resSig_24h_48h.df$ensembl_gene_id <- rownames(resSig_24h_48h.df)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = resSig_24h_48h.df$ensembl_gene_id,
  mart = ensembl
)

resSig_24h_48h_df_annot <- merge(resSig_24h_48h.df, gene_map, by = "ensembl_gene_id", all.x = TRUE)

#Done for all significant DEGs dataframes but not shown

#########################################
#Data visualisation
#########################################

#Transform data using variance stabilising transformation method for PCA plot
vst_transformation <- varianceStabilizingTransformation(dds,blind = FALSE)

##############
#PCA plots
##############

#If PCA plot for publication, can remove batch effects before PCA using below method
#batch_corrected <- assay(vst_transformation)
#batch_corrected <- removeBatchEffect(batch_corrected,
#batch = vst_transformation$lot,
#batch2 = vst_transformation$passage)
#assay(vst_transformation) <- batch_corrected
#plotPCA(vst_transformation, intgroup = "condition") + ylim(-5, 5)

#Otherwise standard PCA will suffice

#Standard PCA
pcaData <- plotPCA(vst_transformation, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$SampleName <- rownames(pcaData)

#Plot PCA using ggplot2
ggplot(pcaData, aes(PC1, PC2, color = condition, label = SampleName)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
dev.off()

##############
#Volcano plot
##############

#pbs vs 6h
volcano_pbs_6h <- res_pbs_vs_6h_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )
dev.off()

label_data <- volcano_pbs_6h %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   # lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_pbs_6h, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = label_data,
    aes(label = hgnc_symbol),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "PBS vs 6h", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

volcano_pbs_24h <- res_pbs_vs_24h_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )

label_data <- volcano_pbs_24h %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   # lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_pbs_24h, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = label_data,
    aes(label = hgnc_symbol),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "PBS vs 24h", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

volcano_pbs_48h <- res_pbs_vs_48h_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )

label_data <- volcano_pbs_48h %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   # lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_pbs_48h, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = label_data,
    aes(label = hgnc_symbol),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "PBS vs 48h", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

volcano_6h_24h <- res_6h_vs_24h_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )

label_data <- volcano_6h_24h %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   # lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_6h_24h, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = label_data,
    aes(label = hgnc_symbol),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "6h vs 24h", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

volcano_6h_48h <- res_6h_vs_48h_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )

label_data <- volcano_6h_48h %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   # lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_6h_48h, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = label_data,
    aes(label = hgnc_symbol),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "6h vs 48h", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

volcano_24h_48h <- res_24h_vs_48h_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )

label_data <- volcano_24h_48h %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   # lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_24h_48h, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = label_data,
    aes(label = hgnc_symbol),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "24h vs 48h", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

#########################################
#Gene set enrichment analysis (GSEA)
#########################################
#ClusterProfiler GSEA general preparation

#pbs vs 6h
#Create vector
geneList <- res_pbs_vs_6h_shrunk_df_annot$log2FoldChange
#names(geneList) <- rownames(res_pbs_vs_6h_shrunk_df_annot$ensembl_gene_id)
names(geneList) <- res_pbs_vs_6h_shrunk_df_annot$ensembl_gene_id  # no rownames()
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

#Run GSEA for BP (Biological Process)
gsea_pbs_6h <- gseGO(geneList=geneList,
                     ont="ALL",
                     keyType = "ENSEMBL",
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = "org.Hs.eg.db",
                     eps = 1e-10,
                     pAdjustMethod = "hochberg")

#Convert to dataframe and save
gsea_res1_summary <- data.frame(gsea_pbs_6h) 
write.csv(gsea_res1_summary, "PBS_6h_GSEA_GO.csv") 

#pbs vs 24h
geneList <- res_pbs_vs_24h_shrunk_df_annot$log2FoldChange
names(geneList) <- res_pbs_vs_24h_shrunk_df_annot$ensembl_gene_id
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

#Run GSEA for BP (Biological Process)
gsea_pbs_24h <- gseGO(geneList=geneList,
                      ont="ALL",
                      keyType = "ENSEMBL",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      verbose = TRUE,
                      OrgDb = "org.Hs.eg.db",
                      eps = 1e-10,
                      pAdjustMethod = "hochberg")

#Convert to dataframe and save
gsea_res2_summary <- data.frame(gsea_pbs_24h) 
write.csv(gsea_res2_summary, "PBS_24h_GSEA_GO.csv") 

#pbs vs 48h
geneList <- res_pbs_vs_48h_shrunk_df_annot$log2FoldChange
names(geneList) <- res_pbs_vs_48h_shrunk_df_annot$ensembl_gene_id
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

#Run GSEA for BP (Biological Process)
gsea_pbs_48h <- gseGO(geneList=geneList,
                      ont="ALL",
                      keyType = "ENSEMBL",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      verbose = TRUE,
                      OrgDb = "org.Hs.eg.db",
                      eps = 1e-10,
                      pAdjustMethod = "hochberg")

#Convert to dataframe and save
gsea_res3_summary <- data.frame(gsea_pbs_48h) 
write.csv(gsea_res3_summary, "PBS_48h_GSEA_GO.csv") 

#6h vs 24h
geneList <- res_6h_vs_24h_shrunk_df_annot$log2FoldChange
names(geneList) <- res_6h_vs_24h_shrunk_df_annot$ensembl_gene_id
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

#Run GSEA for BP (Biological Process)
gsea_6h_24h <- gseGO(geneList=geneList,
                     ont="ALL",
                     keyType = "ENSEMBL",
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = "org.Hs.eg.db",
                     eps = 1e-10,
                     pAdjustMethod = "hochberg")

#Convert to dataframe and save
gsea_res4_summary <- data.frame(gsea_6h_24h) 
write.csv(gsea_res4_summary, "6h_24h_GSEA_GO.csv") 

#6h vs 48h
geneList <- res_6h_vs_48h_shrunk_df_annot$log2FoldChange
names(geneList) <- res_6h_vs_48h_shrunk_df_annot$ensembl_gene_id
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

#Run GSEA for BP (Biological Process)
gsea_6h_48h <- gseGO(geneList=geneList,
                     ont="ALL",
                     keyType = "ENSEMBL",
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = "org.Hs.eg.db",
                     eps = 1e-10,
                     pAdjustMethod = "hochberg")

#Convert to dataframe and save
gsea_res5_summary <- data.frame(gsea_6h_48h) 
write.csv(gsea_res5_summary, "6h_48h_GSEA_GO.csv") 

#24h vs 48h
geneList <- res_24h_vs_48h_shrunk_df_annot$log2FoldChange
names(geneList) <- res_24h_vs_48h_shrunk_df_annot$ensembl_gene_id
geneList <- na.omit(geneList)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

#Run GSEA for BP (Biological Process)
gsea_24h_48h <- gseGO(geneList=geneList,
                      ont="ALL",
                      keyType = "ENSEMBL",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      verbose = TRUE,
                      OrgDb = "org.Hs.eg.db",
                      eps = 1e-10,
                      pAdjustMethod = "hochberg")

#Convert to dataframe and save
gsea_res6_summary <- data.frame(gsea_24h_48h) 
write.csv(gsea_res6_summary, "24h_48h_GSEA_GO.csv")


