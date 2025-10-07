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

#Comparison 1: co-culture vs co-culture w/HBVP
comparison_1 <- c(
  "CC_1_S13_R1", "CC_2_S14_R1", "CC_3_S15_R1", "CC_4_S16_R1", "CC_5_S17_R1", "CC_6_S18_R1",
  "CC_HBVP_1_S19", "CC_HBVP_2_S20", "CC_HBVP_3_S21", "CC_HBVP_4_S22", "CC_HBVP_5_S23", "CC_HBVP_6_S24"
)

comparison_2 <- c(
  "HBVP_C_1_S1", "HBVP_C_2_S3", "HBVP_C_3_S5", "HBVP_C_4_S7", "HBVP_C_5_S9", "HBVP_C_6_S11",
  "HBVP_IR_1_S2", "HBVP_IR_2_S4", "HBVP_IR_3_S6", "HBVP_IR_4_S8", "HBVP_IR_5_S10", "HBVP_IR_6_S12"
)

counts1 <- cleaned_counts[ , comparison_1]
coldata1 <- data.frame(row.names=comparison_1, condition=factor(c(rep("coculture",6), rep("coculture_HBVP",6))))

counts2 <- cleaned_counts[ , comparison_2]
coldata2 <- data.frame(row.names=comparison_2, condition=factor(c(rep("HBVP",6), rep("HBVP_IR",6))))

#Create experiment labels
experiment_1 <- factor(c(rep("exp1", 3), rep("exp2", 3), rep("exp1", 3), rep("exp2", 3)))
experiment_2 <- factor(rep(c("exp1", "exp2", "exp3", "exp4", "exp5", "exp6"), times = 2))

#Create condition
condition_1 <- factor(c(rep("CC", 6), rep("CC_HBVP", 6)))
condition_2 <- factor(c(rep("HBVP", 6), rep("HBVP_IR", 6)))

#Create metadata
coldata1 <- data.frame(
  row.names = comparison_1,
  condition = condition_1,
  experiment = experiment_1
)

coldata2 <- data.frame(
  row.names = comparison_2,
  condition = condition_2,
  experiment = experiment_2
)

rownames(coldata1)
rownames(coldata2)

#########################################
#Perform comparisons and DESeq2 analysis
#########################################

#Create dds with design as condition and experiment
dds1_exp <- DESeqDataSetFromMatrix(countData = counts1, colData = coldata1, design = ~  experiment + condition)
dds2_exp <- DESeqDataSetFromMatrix(countData = counts2, colData = coldata2, design = ~  experiment + condition)

#Prioritise results with at least 10 or more reads in 6 samples (this ensures DEG is present in all repeats)
#Thus higher confidence DEGs
keep1 <-rowSums(counts(dds1_exp) >= 10) >= 6
dds1_exp <- dds1_exp[keep1, ]
dds1_exp <- DESeq(dds1_exp)
summary(dds1_exp)

#DESeq2 automatically performs these, however they can be applied just in case
ddsSF1 <- estimateSizeFactors(dds1_exp)
ddsED1 <- estimateDispersions(ddsSF1)
ddsnb1 <- nbinomWaldTest(ddsED1, maxit = 1000) 

#Check comparisons
resultsNames(ddsnb1)

#State comparison
res_cc_vs_cchbvp  <- results(ddsnb1, contrast = c("condition", "CC_HBVP", "CC"))
summary_res1 <-summary(res_cc_vs_cchbvp)

#Same process again for second experiment
keep2 <-rowSums(counts(dds2_exp) >= 10) >= 6
dds2_exp <- dds2_exp[keep2, ]
dds2_exp <- DESeq(dds2_exp)
summary(dds2_exp)

ddsSF2 <- estimateSizeFactors(dds2_exp)
ddsED2 <- estimateDispersions(ddsSF2)
ddsnb2 <- nbinomWaldTest(ddsED2, maxit = 1000) 

resultsNames(ddsnb2)

res_hbvpir_vs_hbvp  <- results(ddsnb2, contrast = c("condition", "HBVP_IR", "HBVP"))
summary_res2 <-summary(res_hbvpir_vs_hbvp)

#Summarise results 
res1_exp <- results(ddsnb1)
res2_exp <- results(ddsnb2)
summary_res1 <-summary(res1_exp)
summary_res2 <-summary(res2_exp)

#########################################
#Shrink data
#########################################

#Plot MA plots for each comparison to check if data needs shrinking
plotMA(res1_exp, ylim = c(-1.5, 1.5), main="CC_vs_CC_HBVP (before shrinkage)")
res1_shrunk <- lfcShrink(ddsnb1, coef="condition_CC_HBVP_vs_CC", type="apeglm")
plotMA(res1_shrunk, ylim=c(-1.5,1.5), main="CC_vs_CC_HBVP (after shrinkage)")

plotMA(res2_exp, ylim = c(-1.5, 1.5), main="HBVP_vs_HBVP_IR (before shrinkage)")
res2_shrunk <- lfcShrink(ddsnb2, coef="condition_HBVP_IR_vs_HBVP", type="apeglm")
plotMA(res2_shrunk, ylim=c(-1.5,1.5), main="HBVP_vs_HBVP_IR (after shrinkage)")

#########################################
#Identify significant DEGs
#########################################

#Pull out significant (padj < 0.05) differentially expressed genes (DEGs) from shrunk dataset 
resSig_1 <- res1_shrunk[which(res1_shrunk$padj < 0.05), ]
head(resSig_1[order(resSig_1$padj), ])
summary(resSig_1)
resSig_1.df <- as.data.frame(resSig_1)
resSig_1 <- (resSig_1.df[order(resSig_1.df$padj), ])

#Save DEGs
write.csv(resSig_1, "Sig_DEGs_res1_shrunk.csv", row.names = TRUE)

resSig_2 <- res2_shrunk[which(res2_shrunk$padj < 0.05), ]
head(resSig_2[order(resSig_2$padj), ])
summary(resSig_2)
resSig_2.df <- as.data.frame(resSig_2)
resSig_2 <- (resSig_2.df[order(resSig_2.df$padj), ])

#Save DEGs
write.csv(resSig_2, "Sig_DEGs_res2_shrunk.csv", row.names = TRUE)

#########################################
#Data annotation
#########################################

#Annotate dataframes for gene names instead of ensembl gene IDs
#Create ensembl object
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Annotate all shrunk dataframes first
res1_shrunk_df <- as.data.frame(res1_shrunk)
res1_shrunk_df$ensembl_gene_id <- rownames(res1_shrunk_df)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res1_shrunk_df$ensembl_gene_id,
  mart = ensembl
)

res1_shrunk_df_annot <- merge(res1_shrunk_df, gene_map, by = "ensembl_gene_id", all.x = TRUE)

res2_shrunk_df <- as.data.frame(res2_shrunk)
res2_shrunk_df$ensembl_gene_id <- rownames(res2_shrunk_df)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = res2_shrunk_df$ensembl_gene_id,
  mart = ensembl
)

res2_shrunk_df_annot <- merge(res2_shrunk_df, gene_map, by = "ensembl_gene_id", all.x = TRUE)

#Annotate all significant DEGs dataframes
resSig_1.df$ensembl_gene_id <- rownames(resSig_1.df)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = resSig_1.df$ensembl_gene_id,
  mart = ensembl
)

resSig_1_df_annot <- merge(resSig_1.df, gene_map, by = "ensembl_gene_id", all.x = TRUE)

resSig_2.df$ensembl_gene_id <- rownames(resSig_2.df)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = resSig_2.df$ensembl_gene_id,
  mart = ensembl
)

resSig_2_df_annot <- merge(resSig_2.df, gene_map, by = "ensembl_gene_id", all.x = TRUE)

#########################################
#Data visualisation
#########################################

#Transform data using variance stabilising transformation method for PCA plot
vst_transformation1 <- varianceStabilizingTransformation(ddsnb1,blind = FALSE)
vst_transformation2 <- varianceStabilizingTransformation(ddsnb2,blind = FALSE)

##############
#PCA plots
##############

pcaData1 <- plotPCA(vst_transformation1, intgroup = "condition", returnData = TRUE)
percentVar1 <- round(100 * attr(pcaData1, "percentVar"))
pcaData1$SampleName <- rownames(pcaData1)

#Plot PCA using ggplot2
ggplot(pcaData1, aes(PC1, PC2, color = condition, label = SampleName)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  xlab(paste0("PC1: ", percentVar1[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar1[2], "% variance")) +
  theme_minimal()

pcaData2 <- plotPCA(vst_transformation2, intgroup = "condition", returnData = TRUE)
percentVar2 <- round(100 * attr(pcaData2, "percentVar"))
pcaData2$SampleName <- rownames(pcaData2)

ggplot(pcaData2, aes(PC1, PC2, color = condition, label = SampleName)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  xlab(paste0("PC1: ", percentVar2[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar2[2], "% variance")) +
  theme_minimal()

##############
#Volcano plot
##############

#Comparison 1
volcano_res1 <- res1_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )

label_data <- volcano_res1 %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   #lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
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
  labs(title = "CC vs CC_HBVP", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

#Comparison 2
volcano_res2 <- res2_shrunk_df_annot %>%
  as.data.frame() %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    sig = factor(sig, levels = c("Not Significant", "Upregulated", "Downregulated"))
  )

label_data <- volcano_res2 %>%
  filter(sig %in% c("Upregulated", "Downregulated")) %>%
  group_by(sig) %>%
  arrange(padj) %>%   #lowest adjusted p-value first
  slice_head(n = 5) %>%
  ungroup()

ggplot(volcano_res2, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
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
  labs(title = "HBVP vs HBVP_IR", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

#########################################
#Gene set enrichment analysis (GSEA)
#########################################
#ClusterProfiler GSEA general preparation

#Create vectors
#Comparison 1
geneList1 <- res1_shrunk_df_annot$log2FoldChange
names(geneList1) <- res1_shrunk_df_annot$ensembl_gene_id 
geneList1 <- na.omit(geneList1)
geneList1 <- sort(geneList1, decreasing = TRUE)
geneList1 <- geneList1[!duplicated(names(geneList1))]

#Run GSEA for BP (Biological Process)
gsea_res1 <- gseGO(geneList=geneList1,
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
gsea_res1_summary <- data.frame(gsea_res1) 
write.csv(gsea_res1_summary, "res1_GSEA_GO.csv") 

#Comparison 2
geneList2 <- res2_shrunk_df_annot$log2FoldChange
names(geneList2) <- res2_shrunk_df_annot$ensembl_gene_id 
geneList2 <- na.omit(geneList2)
geneList2 <- sort(geneList2, decreasing = TRUE)
geneList2 <- geneList2[!duplicated(names(geneList2))]

#Run GSEA for BP (Biological Process)
gsea_res2 <- gseGO(geneList=geneList2,
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
gsea_res2_summary <- data.frame(gsea_res2) 
write.csv(gsea_res2_summary, "res2_GSEA_GO.csv") 
