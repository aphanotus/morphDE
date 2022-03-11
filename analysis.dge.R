rm(list=ls())
gc(verbose = TRUE)

{
  library(tidyverse)
  library(viridis)
  library(DESeq2)
  library(IHW)
  library(pheatmap)
  library(ggpubr)
}

###################
# Import the Kallisto gene count matrix and metadata
###################
# For details of count pre-processing, see count.preprocessing.R

load("combined.counts.rda", verbose = TRUE)

# Set filtering criteria
# Note that DESeq2 recommends only minimal filtering of low-count genes
gene.read.number.cut.off <- 100

###################
# Annotations
###################
{
  ann <- read.delim("final_annotations_lvl0.tsv")
  keep.ann.cols <- c("Query.Sequence","Description","Species","Origin.Database",
                     "Contaminant","Informative","UniProt.KEGG.Terms","UniProt.GO.Biological",
                     "UniProt.GO.Cellular","UniProt.GO.Molecular","EggNOG.Predicted.Gene",
                     "EggNOG.Description","EggNOG.KEGG.Terms","EggNOG.GO.Biological",
                     "EggNOG.GO.Cellular","EggNOG.GO.Molecular","EggNOG.Protein.Domains")
  ann <- ann[,keep.ann.cols]
  ann$Origin.Database <- sub("/research/drangeli/phase2_BCG_RNAseq/Jhae_GU_trinity_assembly/entap_outfiles/similarity_search/DIAMOND/blastp_GU_final_","",ann$Origin.Database)
  ann$Origin.Database <- sub("\\.out","",ann$Origin.Database)
  xenic.ids <- as.vector(unlist(read.delim("final_annotations_contam_id_list.txt", header = FALSE)))
  xenic.ids <- sub(">","",xenic.ids)
  ann$xenic <- ann$Query.Sequence %in% xenic.ids
}

apply.annotation <- function(df, ann) {
  i <- match(rownames(df), ann$Query.Sequence)
  df$Description <- ann$Description[i]
  df$Species <- ann$Species[i]
  df$Contaminant <- ann$Contaminant[i]
  df$xenic <- ann$xenic[i]
  df$UniProt.KEGG.Terms <- ann$UniProt.KEGG.Terms[i]
  df$UniProt.GO.Biological <- ann$UniProt.GO.Biological[i]
  df$UniProt.GO.Cellular <- ann$UniProt.GO.Cellular[i]
  df$UniProt.GO.Molecular <- ann$UniProt.GO.Molecular[i]
  df$EggNOG.Predicted.Gene <- ann$EggNOG.Predicted.Gene[i]
  df$EggNOG.Description <- ann$EggNOG.Description[i]
  df$EggNOG.KEGG.Terms <- ann$EggNOG.KEGG.Terms[i]
  df$EggNOG.GO.Biological <- ann$EggNOG.GO.Biological[i]
  df$EggNOG.GO.Cellular <- ann$EggNOG.GO.Cellular[i]
  df$EggNOG.GO.Molecular <- ann$EggNOG.GO.Molecular[i]
  df$EggNOG.Protein.Domains <- ann$EggNOG.Protein.Domains[i]
  return(df)
}

###################
# Differential expression analysis
###################

# See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

### Null model

#### Sub-setting the data and metadata
cts.i <- as.matrix(cts[which(rowSums(cts) > gene.read.number.cut.off),])
dim(cts.i)

# Why a null model accounting for batch is critical!
pca.uncorrected <- prcomp(t(cts.i))
(pc.variance <- head(borealis::pcvar(pca.uncorrected)))

pca.uncorrected.by.biol.grp.plot <- pca.uncorrected$x %>%
  as.data.frame() %>%
  ggplot(aes(x=PC1, y=PC2)) +
  theme_bw() + theme(legend.position="bottom") +
  geom_point(aes(color=grp), size = 3, alpha = 0.85) +
  scale_color_viridis(name = NULL, discrete = TRUE, begin = 0, end = 1) +
  labs(x=paste0("PC1 (",pc.variance[1]," variance)"),
       y=paste0("PC2 (",pc.variance[2]," variance)")) +
  coord_fixed()
pca.uncorrected.by.biol.grp.plot

pca.uncorrected.by.batch.plot <- pca.uncorrected$x %>%
  as.data.frame() %>%
  ggplot(aes(x=PC1, y=PC2)) +
  theme_bw() + theme(legend.position="bottom") +
  geom_point(aes(color=as.factor(as.numeric(sample.metadata$plate))), size = 3, alpha = 0.85) +
  scale_color_viridis(name = "batch", discrete = TRUE, begin = 0.2, end = 0.8, option = "magma") +
  labs(x=paste0("PC1 (",pc.variance[1]," variance)"),
       y=paste0("PC2 (",pc.variance[2]," variance)")) +
  coord_fixed()
pca.uncorrected.by.batch.plot

pca.uncorrected.plots <- ggpubr::ggarrange(
  pca.uncorrected.by.biol.grp.plot, pca.uncorrected.by.batch.plot,
  ncol = 2)

ggsave("plots/pca.uncorrected.plots.pdf", pca.uncorrected.plots, width = 14, height = 5, scale = 1)
ggsave("plots/pca.uncorrected.plots.jpg", pca.uncorrected.plots, width = 14, height = 5, scale = 1)

# DESeq2
dds.null <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=sample.metadata,
  design = ~ plate )
dds.null <- DESeq(dds.null)
dds.null <- estimateSizeFactors(dds.null)

# # Extract and save normalized counts
# ctn.null <- counts(dds.null, normalized=TRUE)
# write.csv(ctn.null, 'normalized.counts.null.csv')
# # Use normalized counts for downstream visualization

vsd.null <- vst(dds.null, blind=FALSE)

### hierarchical clustering
vsd.cor.null <- cor(assay(vsd.null))
colnames(sample.metadata)
sample.metadata %>% select(7,3,4,8,5,2) %>%
  pheatmap(vsd.cor.null, annotation = ., labels_row = " ", labels_col = " ")

# PCA
pca.null <- prcomp(t(assay(vsd.null)))
(pc.variance <- head(borealis::pcvar(pca.null)))

pca.null.plot <- pca.null$x %>%
  as.data.frame() %>%
  ggplot(aes(x=PC1, y=PC2)) +
  theme_bw() + theme(legend.position="bottom") +
  geom_point(aes(color=grp), size = 3, alpha = 0.85) +
  scale_color_viridis(name = NULL, discrete = TRUE, begin = 0, end = 1) +
  labs(x=paste0("PC1 (",pc.variance[1]," variance)"),
       y=paste0("PC2 (",pc.variance[2]," variance)")) +
  coord_fixed()
pca.null.plot

ggsave("plots/pca.null.plot.pdf", pca.null.plot, width = 7, height = 5, scale = 1)
ggsave("plots/pca.null.plot.jpg", pca.null.plot, width = 7, height = 5, scale = 0.95)

rm(cts.i, ctn.null, dds.null, pca.uncorrected, pca.uncorrected.by.biol.grp.plot, pca.uncorrected.by.batch.plot, vsd.cor.null, vsd.null, pca.null)

###################
### M v. F for adult gonad
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.gonad.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.adult.gonad.by.sex <- DESeq(dds.adult.gonad.by.sex)
dds.adult.gonad.by.sex <- estimateSizeFactors(dds.adult.gonad.by.sex)

# # Extract and save normalized counts
# ctn.adult.gonad.by.sex <- counts(dds.adult.gonad.by.sex, normalized=TRUE)
# write.csv(ctn.adult.gonad.by.sex, 'normalized.counts.adult.gonad.by.sex.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.sex)
res.adult.gonad.by.sex <- results(dds.adult.gonad.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.sex <- lfcShrink(dds.adult.gonad.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.adult.gonad.by.sex)
summary(res.adult.gonad.by.sex)
# out of 54709 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 35628, 65%
# LFC < 0 (down)     : 10205, 19%
# outliers [1]       : 0, 0%

# plotMA(res.adult.gonad.by.sex, main="apeGLM")
# # interactively detect the row number
# i <- identify(res.adult.gonad.by.sex$baseMean, res.adult.gonad.by.sex$log2FoldChange)
# plotCounts(dds.adult.gonad.by.sex, gene=i[1], intgroup="sex", main="adult.gonad.by.sex")
# plotCounts(dds.adult.gonad.by.sex, gene=which.min(res.adult.gonad.by.sex$padj), intgroup="sex", main="adult.gonad.by.sex")
# # Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.sex <- vst(dds.adult.gonad.by.sex, blind=FALSE)
plotPCA(vsd.adult.gonad.by.sex, intgroup=c("sex"))
plotPCA(vsd.adult.gonad.by.sex, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.sex <- res.adult.gonad.by.sex[order(res.adult.gonad.by.sex$padj),]
deg100.adult.gonad.by.sex <- deg100.adult.gonad.by.sex[1:100,]
deg100.adult.gonad.by.sex <- apply.annotation(deg100.adult.gonad.by.sex, ann)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.sex), "deg100.adult.gonad.by.sex.tsv", quote = FALSE)
save(meta.i, cts.i, dds.adult.gonad.by.sex, res.adult.gonad.by.sex, vsd.adult.gonad.by.sex, deg100.adult.gonad.by.sex,
     file = "adult.gonad.by.sex.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.sex, res.adult.gonad.by.sex, vsd.adult.gonad.by.sex)
# load("adult.gonad.by.sex.rda", verbose = TRUE)

###################
### M v. F for L5 gonad
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "L5", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.L5.gonad.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.L5.gonad.by.sex <- DESeq(dds.L5.gonad.by.sex)
dds.L5.gonad.by.sex <- estimateSizeFactors(dds.L5.gonad.by.sex)

# # Extract and save normalized counts
# ctn.L5.gonad.by.sex <- counts(dds.L5.gonad.by.sex, normalized=TRUE)
# write.csv(ctn.L5.gonad.by.sex, 'normalized.counts.L5.gonad.by.sex.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.sex)
res.L5.gonad.by.sex <- results(dds.L5.gonad.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.sex <- lfcShrink(dds.L5.gonad.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.L5.gonad.by.sex)
summary(res.L5.gonad.by.sex)
# out of 32596 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 7179, 22%
# LFC < 0 (down)     : 2367, 7.3%
# outliers [1]       : 114, 0.35%

# plotMA(res.L5.gonad.by.sex, main="apeGLM")
# # Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.sex <- vst(dds.L5.gonad.by.sex, blind=FALSE)
plotPCA(vsd.L5.gonad.by.sex, intgroup=c("sex"))

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.sex <- res.L5.gonad.by.sex[order(res.L5.gonad.by.sex$padj),]
deg100.L5.gonad.by.sex <- deg100.L5.gonad.by.sex[1:100,]
deg100.L5.gonad.by.sex <- apply.annotation(deg100.L5.gonad.by.sex, ann)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.sex), "deg100.L5.gonad.by.sex.tsv", quote = FALSE)
save(meta.i, cts.i, dds.L5.gonad.by.sex, res.L5.gonad.by.sex, vsd.L5.gonad.by.sex, deg100.L5.gonad.by.sex,
     file = "L5.gonad.by.sex.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.sex, res.L5.gonad.by.sex, vsd.L5.gonad.by.sex)
# load("L5.gonad.by.sex.rda", verbose = TRUE)

###################
### M v. F for adult thorax
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.thorax.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.adult.thorax.by.sex <- DESeq(dds.adult.thorax.by.sex)
dds.adult.thorax.by.sex <- estimateSizeFactors(dds.adult.thorax.by.sex)

# # Extract and save normalized counts
# ctn.adult.thorax.by.sex <- counts(dds.adult.thorax.by.sex, normalized=TRUE)
# write.csv(ctn.adult.thorax.by.sex, 'normalized.counts.adult.thorax.by.sex.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.sex)
res.adult.thorax.by.sex <- results(dds.adult.thorax.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.sex <- lfcShrink(dds.adult.thorax.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.adult.thorax.by.sex)
summary(res.adult.thorax.by.sex)
# out of 32352 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1244, 3.8%
# LFC < 0 (down)     : 745, 2.3%
# outliers [1]       : 0, 0%

# plotMA(res.adult.thorax.by.sex, main="apeGLM")
# # Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.sex <- vst(dds.adult.thorax.by.sex, blind=FALSE)
plotPCA(vsd.adult.thorax.by.sex, intgroup=c("sex"))

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.sex <- res.adult.thorax.by.sex[order(res.adult.thorax.by.sex$padj),]
deg100.adult.thorax.by.sex <- deg100.adult.thorax.by.sex[1:100,]
deg100.adult.thorax.by.sex <- apply.annotation(deg100.adult.thorax.by.sex, ann)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.sex), "deg100.adult.thorax.by.sex.tsv", quote = FALSE)
save(meta.i, cts.i, dds.adult.thorax.by.sex, res.adult.thorax.by.sex, vsd.adult.thorax.by.sex, deg100.adult.thorax.by.sex,
     file = "adult.thorax.by.sex.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.sex, res.adult.thorax.by.sex, vsd.adult.thorax.by.sex)
# load("adult.thorax.by.sex.rda", verbose = TRUE)

###################
### M v. F for L5 thorax
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "L5", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.L5.thorax.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.L5.thorax.by.sex <- DESeq(dds.L5.thorax.by.sex)
dds.L5.thorax.by.sex <- estimateSizeFactors(dds.L5.thorax.by.sex)

# # Extract and save normalized counts
# ctn.L5.thorax.by.sex <- counts(dds.L5.thorax.by.sex, normalized=TRUE)
# write.csv(ctn.L5.thorax.by.sex, 'normalized.counts.L5.thorax.by.sex.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.sex)
res.L5.thorax.by.sex <- results(dds.L5.thorax.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.sex <- lfcShrink(dds.L5.thorax.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.L5.thorax.by.sex)
summary(res.L5.thorax.by.sex)
# out of 24940 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 25, 0.1%
# LFC < 0 (down)     : 21, 0.084%
# outliers [1]       : 55, 0.22%

plotMA(res.L5.thorax.by.sex, main="apeGLM")
# Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.sex <- vst(dds.L5.thorax.by.sex, blind=FALSE)
plotPCA(vsd.L5.thorax.by.sex, intgroup=c("sex"))

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.sex <- res.L5.thorax.by.sex[order(res.L5.thorax.by.sex$padj),]
deg100.L5.thorax.by.sex <- deg100.L5.thorax.by.sex[1:100,]
deg100.L5.thorax.by.sex <- apply.annotation(deg100.L5.thorax.by.sex, ann)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.sex), "deg100.L5.thorax.by.sex.tsv", quote = FALSE)
save(meta.i, cts.i, dds.L5.thorax.by.sex, res.L5.thorax.by.sex, vsd.L5.thorax.by.sex, deg100.L5.thorax.by.sex,
     file = "L5.thorax.by.sex.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.sex, res.L5.thorax.by.sex, vsd.L5.thorax.by.sex)
# load("L5.thorax.by.sex.rda", verbose = TRUE)

###################
### Adult thorax by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.thorax.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph )
dds.adult.thorax.by.morph <- DESeq(dds.adult.thorax.by.morph)
dds.adult.thorax.by.morph <- estimateSizeFactors(dds.adult.thorax.by.morph)

# # Extract and save normalized counts
# ctn.adult.thorax.by.morph <- counts(dds.adult.thorax.by.morph, normalized=TRUE)
# write.csv(ctn.adult.thorax.by.morph, 'normalized.counts.adult.thorax.by.morph.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.morph)
res.adult.thorax.by.morph <- results(dds.adult.thorax.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.morph <- lfcShrink(dds.adult.thorax.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.adult.thorax.by.morph)
summary(res.adult.thorax.by.morph)
# out of 32356 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1635, 5.1%
# LFC < 0 (down)     : 1364, 4.2%
# outliers [1]       : 67, 0.21%

plotMA(res.adult.thorax.by.morph, main="apeGLM")
# # interactively detect the row number
# i <- identify(res.adult.thorax.by.morph$baseMean, res.adult.thorax.by.morph$log2FoldChange)
# res.adult.thorax.by.morph[i[1],]
# rownames(ctn.adult.thorax.by.morph)[i[1]]
# data.frame(
#   counts = ctn.adult.thorax.by.morph[i[1],],
#   morph = meta.i$morph) %>%
#   ggplot(aes(morph,counts)) +
#   geom_violin(fill="gray85", color = NA, alpha = 0.65) +
#   geom_point(position=position_jitter(w=0.1,h=0), alpha = 0.65) +
#   scale_y_log10()
# # Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.morph <- vst(dds.adult.thorax.by.morph, blind=FALSE)
plotPCA(vsd.adult.thorax.by.morph, intgroup=c("morph"))
# The separation of groups in the NW/SE axis is by plate!

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.morph <- res.adult.thorax.by.morph[order(res.adult.thorax.by.morph$padj),]
deg100.adult.thorax.by.morph <- deg100.adult.thorax.by.morph[1:100,]
deg100.adult.thorax.by.morph <- apply.annotation(deg100.adult.thorax.by.morph, ann)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.morph), "deg100.adult.thorax.by.morph.tsv", quote = FALSE)
save(meta.i, cts.i, dds.adult.thorax.by.morph, res.adult.thorax.by.morph, vsd.adult.thorax.by.morph, deg100.adult.thorax.by.morph,
     file = "adult.thorax.by.morph.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.morph, res.adult.thorax.by.morph, vsd.adult.thorax.by.morph)
# load("adult.thorax.by.morph.rda", verbose = TRUE)

###################
### Adult gonad by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.gonad.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph )
dds.adult.gonad.by.morph <- DESeq(dds.adult.gonad.by.morph)
dds.adult.gonad.by.morph <- estimateSizeFactors(dds.adult.gonad.by.morph)

# # Extract and save normalized counts
# ctn.adult.gonad.by.morph <- counts(dds.adult.gonad.by.morph, normalized=TRUE)
# write.csv(ctn.adult.gonad.by.morph, 'normalized.counts.adult.gonad.by.morph.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.morph)
res.adult.gonad.by.morph <- results(dds.adult.gonad.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.morph <- lfcShrink(dds.adult.gonad.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.adult.gonad.by.morph)
summary(res.adult.gonad.by.morph)
# out of 54717 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3, 0.0055%
# LFC < 0 (down)     : 2, 0.0037%
# outliers [1]       : 119, 0.22%

plotMA(res.adult.gonad.by.morph, main="apeGLM")
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.morph <- vst(dds.adult.gonad.by.morph, blind=FALSE)
plotPCA(vsd.adult.gonad.by.morph, intgroup=c("morph"))
plotPCA(vsd.adult.gonad.by.morph, intgroup=c("sex"))

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.morph <- res.adult.gonad.by.morph[order(res.adult.gonad.by.morph$padj),]
deg100.adult.gonad.by.morph <- deg100.adult.gonad.by.morph[1:100,]
deg100.adult.gonad.by.morph <- apply.annotation(deg100.adult.gonad.by.morph, ann)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.morph), "deg100.adult.gonad.by.morph.tsv", quote = FALSE)
save(meta.i, cts.i, dds.adult.gonad.by.morph, res.adult.gonad.by.morph, vsd.adult.gonad.by.morph, deg100.adult.gonad.by.morph,
     file = "adult.gonad.by.morph.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.morph, res.adult.gonad.by.morph, vsd.adult.gonad.by.morph)
# load("adult.gonad.by.morph.rda", verbose = TRUE)

###################
### Adult gonad by food regime
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.gonad.by.food <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph + food_regime)
dds.adult.gonad.by.food <- DESeq(dds.adult.gonad.by.food)
dds.adult.gonad.by.food <- estimateSizeFactors(dds.adult.gonad.by.food)

# # Extract and save normalized counts
# ctn.adult.gonad.by.food <- counts(dds.adult.gonad.by.food, normalized=TRUE)
# write.csv(ctn.adult.gonad.by.food, 'normalized.counts.adult.gonad.by.food.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.food)
res.adult.gonad.by.food <- results(dds.adult.gonad.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.food <- lfcShrink(dds.adult.gonad.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.adult.gonad.by.food)
summary(res.adult.gonad.by.food)
# out of 54724 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 40, 0.073%
# LFC < 0 (down)     : 7, 0.013%
# outliers [1]       : 430, 0.79%

plotMA(res.adult.gonad.by.food, main="apeGLM", ylim = c(-3,5))
# Positive LFC values are associated with low food regime expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.food <- vst(dds.adult.gonad.by.food, blind=FALSE)
plotPCA(vsd.adult.gonad.by.food, intgroup=c("food_regime"))

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.food <- res.adult.gonad.by.food[order(res.adult.gonad.by.food$padj),]
deg100.adult.gonad.by.food <- deg100.adult.gonad.by.food[1:100,]
deg100.adult.gonad.by.food <- apply.annotation(deg100.adult.gonad.by.food, ann)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.food), "deg100.adult.gonad.by.food.tsv", quote = FALSE)
save(meta.i, cts.i, dds.adult.gonad.by.food, res.adult.gonad.by.food, vsd.adult.gonad.by.food, deg100.adult.gonad.by.food,
     file = "adult.gonad.by.food.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.food, res.adult.gonad.by.food, vsd.adult.gonad.by.food)
# load("adult.gonad.by.food.rda", verbose = TRUE)

###################
### Adult thorax by food regime
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.thorax.by.food <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph + food_regime)
dds.adult.thorax.by.food <- DESeq(dds.adult.thorax.by.food)
dds.adult.thorax.by.food <- estimateSizeFactors(dds.adult.thorax.by.food)

# # Extract and save normalized counts
# ctn.adult.thorax.by.food <- counts(dds.adult.thorax.by.food, normalized=TRUE)
# write.csv(ctn.adult.thorax.by.food, 'normalized.counts.adult.thorax.by.food.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.food)
res.adult.thorax.by.food <- results(dds.adult.thorax.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.food <- lfcShrink(dds.adult.thorax.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.adult.thorax.by.food)
summary(res.adult.thorax.by.food)
# out of 32367 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 0.028%
# LFC < 0 (down)     : 6, 0.019%
# outliers [1]       : 161, 0.5%

plotMA(res.adult.thorax.by.food, main="apeGLM", ylim = c(-3,3) )
# Positive LFC values are associated with low food regime expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.food <- vst(dds.adult.thorax.by.food, blind=FALSE)
plotPCA(vsd.adult.thorax.by.food, intgroup=c("food_regime"))
plotPCA(vsd.adult.thorax.by.food, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.food <- res.adult.thorax.by.food[order(res.adult.thorax.by.food$padj),]
deg100.adult.thorax.by.food <- deg100.adult.thorax.by.food[1:100,]
deg100.adult.thorax.by.food <- apply.annotation(deg100.adult.thorax.by.food, ann)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.food), "deg100.adult.thorax.by.food.tsv", quote = FALSE)
save(meta.i, cts.i, dds.adult.thorax.by.food, res.adult.thorax.by.food, vsd.adult.thorax.by.food, deg100.adult.thorax.by.food,
     file = "adult.thorax.by.food.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.food, res.adult.thorax.by.food, vsd.adult.thorax.by.food)
# load("adult.thorax.by.food.rda", verbose = TRUE)

###################
### L5 gonad by food regime
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "L5", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(food_regime,food_regime,length)))

dds.L5.gonad.by.food <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + food_regime )
dds.L5.gonad.by.food <- DESeq(dds.L5.gonad.by.food)
dds.L5.gonad.by.food <- estimateSizeFactors(dds.L5.gonad.by.food)

# # Extract and save normalized counts
# ctn.L5.gonad.by.food <- counts(dds.L5.gonad.by.food, normalized=TRUE)
# write.csv(ctn.L5.gonad.by.food, 'normalized.counts.L5.gonad.by.food.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.food)
res.L5.gonad.by.food <- results(dds.L5.gonad.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.food <- lfcShrink(dds.L5.gonad.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.L5.gonad.by.food)
summary(res.L5.gonad.by.food)
# out of 32603 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3, 0.0092%
# LFC < 0 (down)     : 24, 0.074%
# outliers [1]       : 151, 0.46%

plotMA(res.L5.gonad.by.food, main="apeGLM", ylim=c(-3,5))
plotCounts(dds.L5.gonad.by.food, gene=order(res.L5.gonad.by.food$padj)[2], intgroup="food_regime")
res.L5.gonad.by.food$log2FoldChange[order(res.L5.gonad.by.food$padj)[2]]
# Positive LFC values are associated with low food expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.food <- vst(dds.L5.gonad.by.food, blind=FALSE)
plotPCA(vsd.L5.gonad.by.food, intgroup=c("food_regime"))
plotPCA(vsd.L5.gonad.by.food, intgroup=c("sex"))

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.food <- res.L5.gonad.by.food[order(res.L5.gonad.by.food$padj),]
deg100.L5.gonad.by.food <- deg100.L5.gonad.by.food[1:100,]
deg100.L5.gonad.by.food <- apply.annotation(deg100.L5.gonad.by.food, ann)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.food), "deg100.L5.gonad.by.food.tsv", quote = FALSE)
save(meta.i, cts.i, dds.L5.gonad.by.food, res.L5.gonad.by.food, vsd.L5.gonad.by.food, deg100.L5.gonad.by.food,
     file = "L5.gonad.by.food.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.food, res.L5.gonad.by.food, vsd.L5.gonad.by.food)
# load("L5.gonad.by.food.rda", verbose = TRUE)

###################
### L5 thorax by food regime
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "L5", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(food_regime,food_regime,length)))

dds.L5.thorax.by.food <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + food_regime )
dds.L5.thorax.by.food <- DESeq(dds.L5.thorax.by.food)
dds.L5.thorax.by.food <- estimateSizeFactors(dds.L5.thorax.by.food)

# # Extract and save normalized counts
# ctn.L5.thorax.by.food <- counts(dds.L5.thorax.by.food, normalized=TRUE)
# write.csv(ctn.L5.thorax.by.food, 'normalized.counts.L5.thorax.by.food.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.food)
res.L5.thorax.by.food <- results(dds.L5.thorax.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.food <- lfcShrink(dds.L5.thorax.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.L5.thorax.by.food)
# summary(res.L5.thorax.by.food)
# out of 24956 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 2, 0.008%
# outliers [1]       : 244, 0.98%

plotMA(res.L5.thorax.by.food, main="apeGLM", ylim = c(-2.5,3.5))
plotCounts(dds.L5.thorax.by.food, gene=order(res.L5.thorax.by.food$padj)[2], intgroup="food_regime")
res.L5.thorax.by.food$log2FoldChange[order(res.L5.thorax.by.food$padj)[2]]
# Positive LFC values are associated with low food expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.food <- vst(dds.L5.thorax.by.food, blind=FALSE)
plotPCA(vsd.L5.thorax.by.food, intgroup=c("food_regime"))
plotPCA(vsd.L5.thorax.by.food, intgroup=c("sex"))

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.food <- res.L5.thorax.by.food[order(res.L5.thorax.by.food$padj),]
deg100.L5.thorax.by.food <- deg100.L5.thorax.by.food[1:100,]
deg100.L5.thorax.by.food <- apply.annotation(deg100.L5.thorax.by.food, ann)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.food), "deg100.L5.thorax.by.food.tsv", quote = FALSE)
save(meta.i, cts.i, dds.L5.thorax.by.food, res.L5.thorax.by.food, vsd.L5.thorax.by.food, deg100.L5.thorax.by.food,
     file = "L5.thorax.by.food.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.food, res.L5.thorax.by.food, vsd.L5.thorax.by.food)
# load("L5.thorax.by.food.rda", verbose = TRUE)

###################
### Ovaries by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad", sex =="f")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))

dds.ovaries.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph )
dds.ovaries.by.morph <- DESeq(dds.ovaries.by.morph)
dds.ovaries.by.morph <- estimateSizeFactors(dds.ovaries.by.morph)

# # Extract and save normalized counts
# ctn.ovaries.by.morph <- counts(dds.ovaries.by.morph, normalized=TRUE)
# write.csv(ctn.ovaries.by.morph, 'normalized.counts.ovaries.by.morph.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.ovaries.by.morph)
res.ovaries.by.morph <- results(dds.ovaries.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.ovaries.by.morph <- lfcShrink(dds.ovaries.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.ovaries.by.morph)
summary(res.ovaries.by.morph)
# out of 26403 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5, 0.019%
# LFC < 0 (down)     : 4, 0.015%
# outliers [1]       : 25, 0.095%

plotMA(res.ovaries.by.morph, main="apeGLM", ylim = c(-4,5))
plotCounts(dds.ovaries.by.morph, gene=order(res.ovaries.by.morph$padj)[8], intgroup="morph")
res.ovaries.by.morph$log2FoldChange[order(res.ovaries.by.morph$padj)[8]]
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.ovaries.by.morph <- vst(dds.ovaries.by.morph, blind=FALSE)
plotPCA(vsd.ovaries.by.morph, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.ovaries.by.morph <- res.ovaries.by.morph[order(res.ovaries.by.morph$padj),]
deg100.ovaries.by.morph <- deg100.ovaries.by.morph[1:100,]
deg100.ovaries.by.morph <- apply.annotation(deg100.ovaries.by.morph, ann)

# Save results
write.table(as.data.frame(deg100.ovaries.by.morph), "deg100.ovaries.by.morph.tsv", quote = FALSE)
save(meta.i, cts.i, dds.ovaries.by.morph, res.ovaries.by.morph, vsd.ovaries.by.morph, deg100.ovaries.by.morph,
     file = "ovaries.by.morph.rda")
rm(meta.i, cts.i, dds.ovaries.by.morph, res.ovaries.by.morph, vsd.ovaries.by.morph)
# load("ovaries.by.morph.rda", verbose = TRUE)

###################
### Testes by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad", sex =="m")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))

dds.testes.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph )
dds.testes.by.morph <- DESeq(dds.testes.by.morph)
dds.testes.by.morph <- estimateSizeFactors(dds.testes.by.morph)

# # Extract and save normalized counts
# ctn.testes.by.morph <- counts(dds.testes.by.morph, normalized=TRUE)
# write.csv(ctn.testes.by.morph, 'normalized.counts.testes.by.morph.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.testes.by.morph)
res.testes.by.morph <- results(dds.testes.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.testes.by.morph <- lfcShrink(dds.testes.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.testes.by.morph)
summary(res.testes.by.morph)
# out of 41078 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.0024%
# LFC < 0 (down)     : 7, 0.017%
# outliers [1]       : 217, 0.53%

plotMA(res.testes.by.morph, main="apeGLM", ylim = c(-4,4.5))
plotCounts(dds.testes.by.morph, gene=order(res.testes.by.morph$padj)[6], intgroup="morph")
res.testes.by.morph$log2FoldChange[order(res.testes.by.morph$padj)[6]]
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.testes.by.morph <- vst(dds.testes.by.morph, blind=FALSE)
plotPCA(vsd.testes.by.morph, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.testes.by.morph <- res.testes.by.morph[order(res.testes.by.morph$padj),]
deg100.testes.by.morph <- deg100.testes.by.morph[1:100,]
deg100.testes.by.morph <- apply.annotation(deg100.testes.by.morph, ann)

# Save results
write.table(as.data.frame(deg100.testes.by.morph), "deg100.testes.by.morph.tsv", quote = FALSE)
save(meta.i, cts.i, dds.testes.by.morph, res.testes.by.morph, vsd.testes.by.morph, deg100.testes.by.morph,
     file = "testes.by.morph.rda")
rm(meta.i, cts.i, dds.testes.by.morph, res.testes.by.morph, vsd.testes.by.morph)
# load("testes.by.morph.rda", verbose = TRUE)

###################
### Thorax by stage
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(stage,stage,length)))

dds.thorax.by.stage <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + stage )
dds.thorax.by.stage <- DESeq(dds.thorax.by.stage)
dds.thorax.by.stage <- estimateSizeFactors(dds.thorax.by.stage)

# # Extract and save normalized counts
# ctn.thorax.by.stage <- counts(dds.thorax.by.stage, normalized=TRUE)
# write.csv(ctn.thorax.by.stage, 'normalized.counts.thorax.by.stage.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.thorax.by.stage)
res.thorax.by.stage <- results(dds.thorax.by.stage, contrast=c("stage","L5","adult"), alpha=0.05, filterFun=ihw)
res.thorax.by.stage <- lfcShrink(dds.thorax.by.stage, coef="stage_L5_vs_adult", type="apeglm", res = res.thorax.by.stage)
summary(res.thorax.by.stage)
# out of 43635 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5830, 13%
# LFC < 0 (down)     : 4286, 9.8%
# outliers [1]       : 30, 0.069%

plotMA(res.thorax.by.stage, main="apeGLM", ylim=c(-10,10))
plotCounts(dds.thorax.by.stage, gene=which.min(res.thorax.by.stage$padj), intgroup="stage")
res.thorax.by.stage$log2FoldChange[which.min(res.thorax.by.stage$padj)]
# Positive LFC values are associated with L5 expression bias

# Use variance stabilizing transformations only for ML applications
vsd.thorax.by.stage <- vst(dds.thorax.by.stage, blind=FALSE)
plotPCA(vsd.thorax.by.stage, intgroup=c("stage"))
plotPCA(vsd.thorax.by.stage, intgroup=c("morph"))
plotPCA(vsd.thorax.by.stage, intgroup=c("food_regime"))

# We can order our results table by the smallest p value:
deg100.thorax.by.stage <- res.thorax.by.stage[order(res.thorax.by.stage$padj),]
deg100.thorax.by.stage <- deg100.thorax.by.stage[1:100,]
deg100.thorax.by.stage <- apply.annotation(deg100.thorax.by.stage, ann)

# Save results
write.table(as.data.frame(deg100.thorax.by.stage), "deg100.thorax.by.stage.tsv", quote = FALSE)
save(meta.i, cts.i, dds.thorax.by.stage, res.thorax.by.stage, vsd.thorax.by.stage, deg100.thorax.by.stage,
     file = "thorax.by.stage.rda")
rm(meta.i, cts.i, dds.thorax.by.stage, res.thorax.by.stage, vsd.thorax.by.stage)
# load("thorax.by.stage.rda", verbose = TRUE)

###################
### Ovaries by stage
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(tissue != "thorax" & sex == "f")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(stage,stage,length)))

dds.ovaries.by.stage <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + stage )
dds.ovaries.by.stage <- DESeq(dds.ovaries.by.stage)
dds.ovaries.by.stage <- estimateSizeFactors(dds.ovaries.by.stage)

# # Extract and save normalized counts
# ctn.ovaries.by.stage <- counts(dds.ovaries.by.stage, normalized=TRUE)
# write.csv(ctn.ovaries.by.stage, 'normalized.counts.ovaries.by.stage.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.ovaries.by.stage)
res.ovaries.by.stage <- results(dds.ovaries.by.stage, contrast=c("stage","L5","adult"), alpha=0.05, filterFun=ihw)
res.ovaries.by.stage <- lfcShrink(dds.ovaries.by.stage, coef="stage_L5_vs_adult", type="apeglm", res = res.ovaries.by.stage)
summary(res.ovaries.by.stage)
# out of 33765 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 16564, 49%
# LFC < 0 (down)     : 8016, 24%
# outliers [1]       : 38, 0.11%

plotMA(res.ovaries.by.stage, main="apeGLM", ylim=c(-10,10))
plotCounts(dds.ovaries.by.stage, gene=which.min(res.ovaries.by.stage$padj), intgroup="stage")
res.ovaries.by.stage$log2FoldChange[which.min(res.ovaries.by.stage$padj)]
# Positive LFC values are associated with L5 expression bias

# Use variance stabilizing transformations only for ML applications
vsd.ovaries.by.stage <- vst(dds.ovaries.by.stage, blind=FALSE)
plotPCA(vsd.ovaries.by.stage, intgroup=c("stage"))
plotPCA(vsd.ovaries.by.stage, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.ovaries.by.stage <- res.ovaries.by.stage[order(res.ovaries.by.stage$padj),]
deg100.ovaries.by.stage <- deg100.ovaries.by.stage[1:100,]
deg100.ovaries.by.stage <- apply.annotation(deg100.ovaries.by.stage, ann)

# Save results
write.table(as.data.frame(deg100.ovaries.by.stage), "deg100.ovaries.by.stage.tsv", quote = FALSE)
save(meta.i, cts.i, dds.ovaries.by.stage, res.ovaries.by.stage, vsd.ovaries.by.stage, deg100.ovaries.by.stage,
     file = "ovaries.by.stage.rda")
rm(meta.i, cts.i, dds.ovaries.by.stage, res.ovaries.by.stage, vsd.ovaries.by.stage)
# load("ovaries.by.stage.rda", verbose = TRUE)

###################
### Testes by stage
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(tissue != "thorax" & sex == "m")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(stage,stage,length)))

dds.testes.by.stage <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + stage )
dds.testes.by.stage <- DESeq(dds.testes.by.stage)
dds.testes.by.stage <- estimateSizeFactors(dds.testes.by.stage)

# # Extract and save normalized counts
# ctn.testes.by.stage <- counts(dds.testes.by.stage, normalized=TRUE)
# write.csv(ctn.testes.by.stage, 'normalized.counts.testes.by.stage.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.testes.by.stage)
res.testes.by.stage <- results(dds.testes.by.stage, contrast=c("stage","L5","adult"), alpha=0.05, filterFun=ihw)
res.testes.by.stage <- lfcShrink(dds.testes.by.stage, coef="stage_L5_vs_adult", type="apeglm", res = res.testes.by.stage)
summary(res.testes.by.stage)
# out of 52226 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 16981, 33%
# LFC < 0 (down)     : 23090, 44%
# outliers [1]       : 42, 0.08%

plotMA(res.testes.by.stage, main="apeGLM", ylim=c(-10,10))
plotCounts(dds.testes.by.stage, gene=which.min(res.testes.by.stage$padj), intgroup="stage")
res.testes.by.stage$log2FoldChange[which.min(res.testes.by.stage$padj)]
# Positive LFC values are associated with L5 expression bias

# Use variance stabilizing transformations only for ML applications
vsd.testes.by.stage <- vst(dds.testes.by.stage, blind=FALSE)
plotPCA(vsd.testes.by.stage, intgroup=c("stage"))
plotPCA(vsd.testes.by.stage, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.testes.by.stage <- res.testes.by.stage[order(res.testes.by.stage$padj),]
deg100.testes.by.stage <- deg100.testes.by.stage[1:100,]
deg100.testes.by.stage <- apply.annotation(deg100.testes.by.stage, ann)

# Save results
write.table(as.data.frame(deg100.testes.by.stage), "deg100.testes.by.stage.tsv", quote = FALSE)
save(meta.i, cts.i, dds.testes.by.stage, res.testes.by.stage, vsd.testes.by.stage, deg100.testes.by.stage,
     file = "testes.by.stage.rda")
rm(meta.i, cts.i, dds.testes.by.stage, res.testes.by.stage, vsd.testes.by.stage)
# load("testes.by.stage.rda", verbose = TRUE) 

###################
# Save the results
###################
save(ann, apply.annotation, gene.read.number.cut.off, 
     dds.adult.gonad.by.food, dds.adult.gonad.by.morph, dds.adult.gonad.by.sex, dds.adult.thorax.by.food, dds.adult.thorax.by.morph, dds.adult.thorax.by.sex, dds.L5.gonad.by.food, dds.L5.gonad.by.sex, dds.L5.thorax.by.food, dds.L5.thorax.by.sex, dds.ovaries.by.morph, dds.ovaries.by.stage, dds.testes.by.morph, dds.testes.by.stage, dds.thorax.by.stage, 
     file = "dge.analysis.dds.rda")
save(ann, apply.annotation, gene.read.number.cut.off, 
     res.adult.gonad.by.food, res.adult.gonad.by.morph, res.adult.gonad.by.sex, res.adult.thorax.by.food, res.adult.thorax.by.morph, res.adult.thorax.by.sex, res.L5.gonad.by.food, res.L5.gonad.by.sex, res.L5.thorax.by.food, res.L5.thorax.by.sex, res.ovaries.by.morph, res.ovaries.by.stage, res.testes.by.morph, res.testes.by.stage, res.thorax.by.stage,
     file = "dge.analysis.res.rda")
save(deg100.adult.gonad.by.food, deg100.adult.gonad.by.morph, deg100.adult.gonad.by.sex, deg100.adult.thorax.by.food, deg100.adult.thorax.by.morph, deg100.adult.thorax.by.sex, deg100.L5.gonad.by.sex, deg100.L5.thorax.by.sex, deg100.L5.gonad.by.food, deg100.L5.thorax.by.food,
     deg100.ovaries.by.morph, deg100.testes.by.morph, deg100.thorax.by.stage, deg100.ovaries.by.stage, deg100.testes.by.stage,
     file = "dge.analysis.deg100s.rda")

###################
# Enrichment Analysis
###################
# Enrichment Analysis for GO, KEGG and EggNOG Protein Domain terms 
# was performed using the hypergeometric test as implemented in stats::phyper
# 
# source('enrichment.R')

load("enrichment.pvalues.rda",verbose = TRUE)

###################
# Volcano Plots
###################
add.manual.annotations <- function(res) {
  manual.ann <- read.delim("manual.annotations.tsv", header = FALSE)
  colnames(manual.ann) <- c("id","name")
  res$manual.ann <- NA
  for (i in 1:(dim(manual.ann)[1])) {
    id.i <- sub("\\|","\\\\\\|",manual.ann$id[i])
    res$manual.ann[grep(id.i,rownames(res))] <- manual.ann$name[i]
  }
  return(res)
}

volcano.plot <- function(
  res, pCutoff = 1e-3, FCcutoff = 1, title = NULL, caption = NULL,
  text.col = "grey15", max.overlaps = 100) 
{
  require(ggplot2)
  require(magrittr)
  require(ggrepel)
  require(paletteer)
  res <- apply.annotation(res, ann)
  
  res$highlight <- (abs(res$log2FoldChange)>FCcutoff & res$padj<pCutoff ) 
  if (!("manual.ann" %in% colnames(res))) {
    cat("Adding annotations.\n")
    res <- add.manual.annotations(res)
  }
  # Over-write those not highlighted in this plot
  plot.labels <- res$manual.ann
  plot.labels[which(!res$highlight)] <- NA
  
  res %>% 
    ggplot(aes(x=log2FoldChange, y=-log10(padj), color=highlight)) +
    theme_bw() +
    theme(
      legend.position="none",
      panel.grid.minor = element_blank() ) +
    geom_point(alpha=0.65) +
    geom_hline(yintercept = -log10(pCutoff), color = "gray75") +
    geom_vline(xintercept = c(-FCcutoff,FCcutoff), color = "gray75") +
    geom_text_repel(aes(label=plot.labels), color=text.col, max.overlaps = max.overlaps) +
    scale_color_manual(values=c("grey40","#A60021")) +
    labs(x="log2 fold change", y = "-log10 p",
         title = title, caption = caption)
} # End function

{
  vplots <- list()
  
  vplots[[1]] <- 
    volcano.plot(res.L5.thorax.by.sex[,1:5], pCutoff = 0.05, 
                 title = "Expression in juvenile thorax ~ batch + sex",
                 caption = "Positive LFC values indicate male expression bias")
  
  vplots[[2]] <- 
    volcano.plot(res.adult.thorax.by.sex[,1:5], 
                 title = "Expression in adult thorax ~ batch + sex",
                 caption = "Positive LFC values indicate male expression bias")
  
  vplots[[3]] <- 
    volcano.plot(res.adult.thorax.by.morph[,1:5], 
                 title = "Expression in adult thorax ~ batch + sex + morph",
                 caption = "Positive LFC values indicate long-wing expression bias")
  
  
  vplots[[4]] <- 
    volcano.plot(res.L5.thorax.by.food[,1:5], pCutoff = 0.05, 
                 title = "Expression in juvenile thorax ~ batch + sex + food regime",
                 caption = "Positive LFC values indicate low food expression bias")
  
  vplots[[5]] <- 
    volcano.plot(res.adult.thorax.by.food[,1:5], pCutoff = 0.05, 
                 title = "Expression in adult thorax ~ batch + sex + morph + food regime",
                 caption = "Positive LFC values indicate low food expression bias")
  
  vplots[[6]] <- 
    volcano.plot(res.L5.gonad.by.sex[,1:5], max.overlaps = 50,
                 title = "Expression in juvenile gonad ~ batch + sex",
                 caption = "Positive LFC values indicate male expression bias")
  
  vplots[[7]] <- 
    volcano.plot(res.adult.gonad.by.sex[,1:5], max.overlaps = 20,
                 title = "Expression in adult gonad ~ batch + sex",
                 caption = "Positive LFC values indicate male expression bias")
  
  vplots[[8]] <- 
    volcano.plot(res.adult.gonad.by.morph[,1:5], pCutoff = 0.05, 
                 title = "Expression in adult gonad ~ batch + sex + morph",
                 caption = "Positive LFC values indicate long-wing expression bias")
  
  vplots[[9]] <- 
    volcano.plot(res.L5.gonad.by.food[,1:5], pCutoff = 0.05, 
                 title = "Expression in juvenile gonad ~ batch + sex + food regime",
                 caption = "Positive LFC values indicate low food expression bias")
  
  vplots[[10]] <- 
    volcano.plot(res.adult.gonad.by.food[,1:5], pCutoff = 0.05, 
                 title = "Expression in adult gonad ~ batch + sex + morph + food regime",
                 caption = "Positive LFC values indicate low food expression bias")
}

bigplot <- ggarrange(plotlist = vplots, ncol = 5, nrow = 2, labels="AUTO")

ggsave("plots/volcano.plots.jpg", plot = bigplot, width = 25, height = 10, scale = 1.2)




