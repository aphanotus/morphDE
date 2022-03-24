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

add.manual.annotations <- function(res) {
  manual.ann <- read.delim("manual.annotations.tsv", header = FALSE)
  colnames(manual.ann) <- c("id","name")
  res.ids <- sub("\\|","\\\\\\|",rownames(res))
  res$manual.ann <-
    unlist(lapply(res.ids, function (x) {
      s <- manual.ann$name[grep(x,manual.ann$id)]
      if (!isTRUE(nchar(s)>0)) { s <- NA }
      return(s)
    }))
  return(res)
}

###################
# Differential expression analysis
###################

# See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

### Null model

#### Sub-setting the data and metadata
cts.i <- as.matrix(cts[which(rowSums(cts) > gene.read.number.cut.off),])
dim(cts.i)

pca.uncorrected <- prcomp(t(cts.i))
# https://github.com/aphanotus/borealis
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

pca.uncorrected.and.null.plots <- ggpubr::ggarrange(
  pca.uncorrected.by.biol.grp.plot, pca.uncorrected.by.batch.plot, pca.null.plot,
  ncol = 3)

ggsave("plots/pca.uncorrected.and.null.plots.jpg", pca.uncorrected.and.null.plots, width = 21, height = 5, scale = 1)
ggsave("plots/pca.uncorrected.and.null.plots.pdf", pca.uncorrected.and.null.plots, width = 21, height = 5, scale = 1)

rm(cts.i, ctn.null, dds.null, pca.uncorrected, pca.uncorrected.by.biol.grp.plot, pca.uncorrected.by.batch.plot, vsd.cor.null, vsd.null, pca.null)

###################
### Adult gonad by sex
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
deg100.adult.gonad.by.sex <- add.manual.annotations(deg100.adult.gonad.by.sex)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.sex), "model.output/deg100.adult.gonad.by.sex.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.gonad.by.sex, res.adult.gonad.by.sex, deg100.adult.gonad.by.sex,
     file = "model.output/adult.gonad.by.sex.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.sex, res.adult.gonad.by.sex, vsd.adult.gonad.by.sex)
# load("adult.gonad.by.sex.rda", verbose = TRUE)

###################
### L5 gonad by sex
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
deg100.L5.gonad.by.sex <- add.manual.annotations(deg100.L5.gonad.by.sex)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.sex), "model.output/deg100.L5.gonad.by.sex.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.gonad.by.sex, res.L5.gonad.by.sex, deg100.L5.gonad.by.sex,
     file = "model.output/L5.gonad.by.sex.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.sex, res.L5.gonad.by.sex, vsd.L5.gonad.by.sex)
# load("L5.gonad.by.sex.rda", verbose = TRUE)

###################
### Adult thorax by sex
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
deg100.adult.thorax.by.sex <- add.manual.annotations(deg100.adult.thorax.by.sex)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.sex), "model.output/deg100.adult.thorax.by.sex.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.thorax.by.sex, res.adult.thorax.by.sex, deg100.adult.thorax.by.sex,
     file = "model.output/adult.thorax.by.sex.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.sex, res.adult.thorax.by.sex, vsd.adult.thorax.by.sex)
# load("adult.thorax.by.sex.rda", verbose = TRUE)

###################
### L5 thorax by sex
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
deg100.L5.thorax.by.sex <- add.manual.annotations(deg100.L5.thorax.by.sex)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.sex), "model.output/deg100.L5.thorax.by.sex.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.thorax.by.sex, res.L5.thorax.by.sex, deg100.L5.thorax.by.sex,
     file = "model.output/L5.thorax.by.sex.rda")
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
deg100.adult.thorax.by.morph <- add.manual.annotations(deg100.adult.thorax.by.morph)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.morph), "model.output/deg100.adult.thorax.by.morph.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.thorax.by.morph, res.adult.thorax.by.morph, deg100.adult.thorax.by.morph,
     file = "model.output/adult.thorax.by.morph.rda")
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
deg100.adult.gonad.by.morph <- add.manual.annotations(deg100.adult.gonad.by.morph)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.morph), "model.output/deg100.adult.gonad.by.morph.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.gonad.by.morph, res.adult.gonad.by.morph, deg100.adult.gonad.by.morph,
     file = "model.output/adult.gonad.by.morph.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.morph, res.adult.gonad.by.morph, vsd.adult.gonad.by.morph)
# load("adult.gonad.by.morph.rda", verbose = TRUE)

###################
### Adult ovaries by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad", sex =="f")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))

dds.adult.ovaries.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph )
dds.adult.ovaries.by.morph <- DESeq(dds.adult.ovaries.by.morph)
dds.adult.ovaries.by.morph <- estimateSizeFactors(dds.adult.ovaries.by.morph)

# # Extract and save normalized counts
# ctn.ovaries.by.morph <- counts(dds.adult.ovaries.by.morph, normalized=TRUE)
# write.csv(ctn.ovaries.by.morph, 'normalized.counts.ovaries.by.morph.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.ovaries.by.morph)
res.adult.ovaries.by.morph <- results(dds.adult.ovaries.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.adult.ovaries.by.morph <- lfcShrink(dds.adult.ovaries.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.adult.ovaries.by.morph)
summary(res.adult.ovaries.by.morph)
# out of 26403 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5, 0.019%
# LFC < 0 (down)     : 4, 0.015%
# outliers [1]       : 25, 0.095%

plotMA(res.adult.ovaries.by.morph, main="apeGLM", ylim = c(-4,5))
plotCounts(dds.adult.ovaries.by.morph, gene=order(res.adult.ovaries.by.morph$padj)[8], intgroup="morph")
res.adult.ovaries.by.morph$log2FoldChange[order(res.adult.ovaries.by.morph$padj)[8]]
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.ovaries.by.morph <- vst(dds.adult.ovaries.by.morph, blind=FALSE)
plotPCA(vsd.adult.ovaries.by.morph, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.ovaries.by.morph <- res.adult.ovaries.by.morph[order(res.adult.ovaries.by.morph$padj),]
deg100.adult.ovaries.by.morph <- deg100.adult.ovaries.by.morph[1:100,]
deg100.adult.ovaries.by.morph <- apply.annotation(deg100.adult.ovaries.by.morph, ann)
deg100.adult.ovaries.by.morph <- add.manual.annotations(deg100.adult.ovaries.by.morph)

# Save results
write.table(as.data.frame(deg100.adult.ovaries.by.morph), "model.output/deg100.adult.ovaries.by.morph.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.ovaries.by.morph, res.adult.ovaries.by.morph, deg100.adult.ovaries.by.morph,
     file = "model.output/ovaries.by.morph.rda")
rm(meta.i, cts.i, dds.adult.ovaries.by.morph, res.adult.ovaries.by.morph, vsd.adult.ovaries.by.morph)
# load("ovaries.by.morph.rda", verbose = TRUE)

###################
### Adult testes by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad", sex =="m")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))

dds.adult.testes.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph )
dds.adult.testes.by.morph <- DESeq(dds.adult.testes.by.morph)
dds.adult.testes.by.morph <- estimateSizeFactors(dds.adult.testes.by.morph)

# # Extract and save normalized counts
# ctn.testes.by.morph <- counts(dds.adult.testes.by.morph, normalized=TRUE)
# write.csv(ctn.testes.by.morph, 'normalized.counts.testes.by.morph.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.testes.by.morph)
res.adult.testes.by.morph <- results(dds.adult.testes.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.adult.testes.by.morph <- lfcShrink(dds.adult.testes.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.adult.testes.by.morph)
summary(res.adult.testes.by.morph)
# out of 41078 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.0024%
# LFC < 0 (down)     : 7, 0.017%
# outliers [1]       : 217, 0.53%

plotMA(res.adult.testes.by.morph, main="apeGLM", ylim = c(-4,4.5))
plotCounts(dds.adult.testes.by.morph, gene=order(res.adult.testes.by.morph$padj)[6], intgroup="morph")
res.adult.testes.by.morph$log2FoldChange[order(res.adult.testes.by.morph$padj)[6]]
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.testes.by.morph <- vst(dds.adult.testes.by.morph, blind=FALSE)
plotPCA(vsd.adult.testes.by.morph, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.testes.by.morph <- res.adult.testes.by.morph[order(res.adult.testes.by.morph$padj),]
deg100.adult.testes.by.morph <- deg100.adult.testes.by.morph[1:100,]
deg100.adult.testes.by.morph <- apply.annotation(deg100.adult.testes.by.morph, ann)
deg100.adult.testes.by.morph <- add.manual.annotations(deg100.adult.testes.by.morph)

# Save results
write.table(as.data.frame(deg100.adult.testes.by.morph), "model.output/deg100.adult.testes.by.morph.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.testes.by.morph, res.adult.testes.by.morph, deg100.adult.testes.by.morph,
     file = "model.output/testes.by.morph.rda")
rm(meta.i, cts.i, dds.adult.testes.by.morph, res.adult.testes.by.morph, vsd.adult.testes.by.morph)
# load("testes.by.morph.rda", verbose = TRUE)

###################
### Adult thorax by wing PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue=="thorax", !is.na(wingPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))
hist(meta.i$wingPC1) # higher PC1 values are associated with more short-wing shapes

dds.adult.thorax.by.wingPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + wingPC1 )
dds.adult.thorax.by.wingPC1 <- DESeq(dds.adult.thorax.by.wingPC1)
dds.adult.thorax.by.wingPC1 <- estimateSizeFactors(dds.adult.thorax.by.wingPC1)

# # Extract and save normalized counts
# ctn.adult.thorax.by.wingPC1 <- counts(dds.adult.thorax.by.wingPC1, normalized=TRUE)
# write.csv(ctn.adult.thorax.by.wingPC1, 'normalized.counts.adult.thorax.by.wingPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.wingPC1)
res.adult.thorax.by.wingPC1 <- results(dds.adult.thorax.by.wingPC1, name=c("wingPC1"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.wingPC1 <- lfcShrink(dds.adult.thorax.by.wingPC1, coef="wingPC1", type="apeglm", res = res.adult.thorax.by.wingPC1)
summary(res.adult.thorax.by.wingPC1)
# out of 32236 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1457, 4.5%
# LFC < 0 (down)     : 1627, 5%
# outliers [1]       : 0, 0%

plotMA(res.adult.thorax.by.wingPC1, main="apeGLM", ylim=c(-22,13))
plot(plotCounts(dds.adult.thorax.by.wingPC1, gene=order(res.adult.thorax.by.wingPC1$padj)[4], intgroup="wingPC1", returnData=TRUE))
res.adult.thorax.by.wingPC1$log2FoldChange[order(res.adult.thorax.by.wingPC1$padj)[4]]
# Positive LFC values are associated with positive wingPC1 values (short-wing bias)

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.wingPC1 <- vst(dds.adult.thorax.by.wingPC1, blind=FALSE)
plotPCA(vsd.adult.thorax.by.wingPC1, intgroup=c("wingPC1"))
plotPCA(vsd.adult.thorax.by.wingPC1, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.wingPC1 <- res.adult.thorax.by.wingPC1[order(res.adult.thorax.by.wingPC1$padj),]
deg100.adult.thorax.by.wingPC1 <- deg100.adult.thorax.by.wingPC1[1:100,]
deg100.adult.thorax.by.wingPC1 <- apply.annotation(deg100.adult.thorax.by.wingPC1, ann)
deg100.adult.thorax.by.wingPC1 <- add.manual.annotations(deg100.adult.thorax.by.wingPC1)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.wingPC1), "model.output/deg100.adult.thorax.by.wingPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.thorax.by.wingPC1, res.adult.thorax.by.wingPC1, deg100.adult.thorax.by.wingPC1,
     file = "model.output/adult.thorax.by.wingPC1.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.wingPC1, res.adult.thorax.by.wingPC1, vsd.adult.thorax.by.wingPC1)
# load("adult.thorax.by.wingPC1.rda", verbose = TRUE)

###################
### Adult gonad by wing PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue=="gonad", !is.na(wingPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))
hist(meta.i$wingPC1) # higher PC1 values are associated with more short-wing shapes

dds.adult.gonad.by.wingPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + wingPC1 )
dds.adult.gonad.by.wingPC1 <- DESeq(dds.adult.gonad.by.wingPC1)
dds.adult.gonad.by.wingPC1 <- estimateSizeFactors(dds.adult.gonad.by.wingPC1)

# # Extract and save normalized counts
# ctn.adult.gonad.by.wingPC1 <- counts(dds.adult.gonad.by.wingPC1, normalized=TRUE)
# write.csv(ctn.adult.gonad.by.wingPC1, 'normalized.counts.adult.gonad.by.wingPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.wingPC1)
res.adult.gonad.by.wingPC1 <- results(dds.adult.gonad.by.wingPC1, name=c("wingPC1"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.wingPC1 <- lfcShrink(dds.adult.gonad.by.wingPC1, coef="wingPC1", type="apeglm", res = res.adult.gonad.by.wingPC1)
summary(res.adult.gonad.by.wingPC1)
# out of 54574 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 11, 0.02%
# LFC < 0 (down)     : 15, 0.027%
# outliers [1]       : 0, 0%

plotMA(res.adult.gonad.by.wingPC1, main="apeGLM", ylim=c(-5.5,9))
plot(plotCounts(dds.adult.gonad.by.wingPC1, gene=order(res.adult.gonad.by.wingPC1$padj)[2], intgroup="wingPC1", returnData=TRUE))
res.adult.gonad.by.wingPC1$log2FoldChange[order(res.adult.gonad.by.wingPC1$padj)[2]]
# Positive LFC values are associated with positive wingPC1 values (short-wing bias)

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.wingPC1 <- vst(dds.adult.gonad.by.wingPC1, blind=FALSE)
plotPCA(vsd.adult.gonad.by.wingPC1, intgroup=c("wingPC1"))
plotPCA(vsd.adult.gonad.by.wingPC1, intgroup=c("morph"))
plotPCA(vsd.adult.gonad.by.wingPC1, intgroup=c("population"))
plotPCA(vsd.adult.gonad.by.wingPC1, intgroup=c("plate"))

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.wingPC1 <- res.adult.gonad.by.wingPC1[order(res.adult.gonad.by.wingPC1$padj),]
deg100.adult.gonad.by.wingPC1 <- deg100.adult.gonad.by.wingPC1[1:100,]
deg100.adult.gonad.by.wingPC1 <- apply.annotation(deg100.adult.gonad.by.wingPC1, ann)
deg100.adult.gonad.by.wingPC1 <- add.manual.annotations(deg100.adult.gonad.by.wingPC1)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.wingPC1), "model.output/deg100.adult.gonad.by.wingPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.gonad.by.wingPC1, res.adult.gonad.by.wingPC1, deg100.adult.gonad.by.wingPC1,
     file = "model.output/adult.gonad.by.wingPC1.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.wingPC1, res.adult.gonad.by.wingPC1, vsd.adult.gonad.by.wingPC1)
# load("adult.gonad.by.wingPC1.rda", verbose = TRUE)

###################
### Adult ovaries by wing PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue=="gonad", sex=="f", !is.na(wingPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))
hist(meta.i$wingPC1) # higher PC1 values are associated with more short-wing shapes

dds.adult.ovaries.by.wingPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + wingPC1 )
dds.adult.ovaries.by.wingPC1 <- DESeq(dds.adult.ovaries.by.wingPC1)
dds.adult.ovaries.by.wingPC1 <- estimateSizeFactors(dds.adult.ovaries.by.wingPC1)

# # Extract and save normalized counts
# ctn.ovaries.by.wingPC1 <- counts(dds.adult.ovaries.by.wingPC1, normalized=TRUE)
# write.csv(ctn.ovaries.by.wingPC1, 'normalized.counts.ovaries.by.wingPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.ovaries.by.wingPC1)
res.adult.ovaries.by.wingPC1 <- results(dds.adult.ovaries.by.wingPC1, name=c("wingPC1"), alpha=0.05, filterFun=ihw)
res.adult.ovaries.by.wingPC1 <- lfcShrink(dds.adult.ovaries.by.wingPC1, coef="wingPC1", type="apeglm", res = res.adult.ovaries.by.wingPC1)
summary(res.adult.ovaries.by.wingPC1)
# out of 26241 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 10, 0.038%
# LFC < 0 (down)     : 10, 0.038%
# outliers [1]       : 0, 0%

plotMA(res.adult.ovaries.by.wingPC1, main="apeGLM")
plot(plotCounts(dds.adult.ovaries.by.wingPC1, gene=order(res.adult.ovaries.by.wingPC1$padj)[2], intgroup="wingPC1", returnData=TRUE))
res.adult.ovaries.by.wingPC1$log2FoldChange[order(res.adult.ovaries.by.wingPC1$padj)[2]]
# Positive LFC values are associated with positive wingPC1 values (short-wing bias)

# Use variance stabilizing transformations only for ML applications
vsd.adult.ovaries.by.wingPC1 <- vst(dds.adult.ovaries.by.wingPC1, blind=FALSE)
plotPCA(vsd.adult.ovaries.by.wingPC1, intgroup=c("wingPC1"))

# We can order our results table by the smallest p value:
deg100.adult.ovaries.by.wingPC1 <- res.adult.ovaries.by.wingPC1[order(res.adult.ovaries.by.wingPC1$padj),]
deg100.adult.ovaries.by.wingPC1 <- deg100.adult.ovaries.by.wingPC1[1:100,]
deg100.adult.ovaries.by.wingPC1 <- apply.annotation(deg100.adult.ovaries.by.wingPC1, ann)
deg100.adult.ovaries.by.wingPC1 <- add.manual.annotations(deg100.adult.ovaries.by.wingPC1)

# Save results
write.table(as.data.frame(deg100.adult.ovaries.by.wingPC1), "model.output/deg100.adult.ovaries.by.wingPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.ovaries.by.wingPC1, res.adult.ovaries.by.wingPC1, deg100.adult.ovaries.by.wingPC1,
     file = "model.output/adult.ovaries.by.wingPC1.rda")
rm(meta.i, cts.i, dds.adult.ovaries.by.wingPC1, res.adult.ovaries.by.wingPC1, vsd.adult.ovaries.by.wingPC1)
# load("adult.ovaries.by.wingPC1.rda", verbose = TRUE)

###################
### Adult testes by wing PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue=="gonad", sex=="m", !is.na(wingPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))
hist(meta.i$wingPC1) # higher PC1 values are associated with more short-wing shapes

dds.adult.testes.by.wingPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + wingPC1 )
dds.adult.testes.by.wingPC1 <- DESeq(dds.adult.testes.by.wingPC1)
dds.adult.testes.by.wingPC1 <- estimateSizeFactors(dds.adult.testes.by.wingPC1)

# # Extract and save normalized counts
# ctn.testes.by.wingPC1 <- counts(dds.adult.testes.by.wingPC1, normalized=TRUE)
# write.csv(ctn.testes.by.wingPC1, 'normalized.counts.testes.by.wingPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.testes.by.wingPC1)
res.adult.testes.by.wingPC1 <- results(dds.adult.testes.by.wingPC1, name=c("wingPC1"), alpha=0.05, filterFun=ihw)
res.adult.testes.by.wingPC1 <- lfcShrink(dds.adult.testes.by.wingPC1, coef="wingPC1", type="apeglm", res = res.adult.testes.by.wingPC1)
summary(res.adult.testes.by.wingPC1)
# out of 41085 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 7, 0.017%
# LFC < 0 (down)     : 45, 0.11%
# outliers [1]       : 0, 0%

plotMA(res.adult.testes.by.wingPC1, main="apeGLM")
plot(plotCounts(dds.adult.testes.by.wingPC1, gene=order(res.adult.testes.by.wingPC1$padj)[1], intgroup="wingPC1", returnData=TRUE))
res.adult.testes.by.wingPC1$log2FoldChange[order(res.adult.testes.by.wingPC1$padj)[1]]
# Positive LFC values are associated with positive wingPC1 values (short-wing bias)

# Use variance stabilizing transformations only for ML applications
vsd.adult.testes.by.wingPC1 <- vst(dds.adult.testes.by.wingPC1, blind=FALSE)
plotPCA(vsd.adult.testes.by.wingPC1, intgroup=c("wingPC1"))
plotPCA(vsd.adult.testes.by.wingPC1, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.testes.by.wingPC1 <- res.adult.testes.by.wingPC1[order(res.adult.testes.by.wingPC1$padj),]
deg100.adult.testes.by.wingPC1 <- deg100.adult.testes.by.wingPC1[1:100,]
deg100.adult.testes.by.wingPC1 <- apply.annotation(deg100.adult.testes.by.wingPC1, ann)
deg100.adult.testes.by.wingPC1 <- add.manual.annotations(deg100.adult.testes.by.wingPC1)

# Save results
write.table(as.data.frame(deg100.adult.testes.by.wingPC1), "model.output/deg100.adult.testes.by.wingPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.testes.by.wingPC1, res.adult.testes.by.wingPC1, deg100.adult.testes.by.wingPC1,
     file = "model.output/adult.testes.by.wingPC1.rda")
rm(meta.i, cts.i, dds.adult.testes.by.wingPC1, res.adult.testes.by.wingPC1, vsd.adult.testes.by.wingPC1)
# load("adult.testes.by.wingPC1.rda", verbose = TRUE)

###################
### L5 thorax by wingpad PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue=="thorax", !is.na(wingpadPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(food_regime,food_regime,length)))
hist(meta.i$wingpadPC1) # higher PC1 values are associated with more short-wing shapes

dds.L5.thorax.by.wingpadPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + wingpadPC1 )
dds.L5.thorax.by.wingpadPC1 <- DESeq(dds.L5.thorax.by.wingpadPC1)
dds.L5.thorax.by.wingpadPC1 <- estimateSizeFactors(dds.L5.thorax.by.wingpadPC1)

# # Extract and save normalized counts
# ctn.L5.thorax.by.wingpadPC1 <- counts(dds.L5.thorax.by.wingpadPC1, normalized=TRUE)
# write.csv(ctn.L5.thorax.by.wingpadPC1, 'normalized.counts.L5.thorax.by.wingpadPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.wingpadPC1)
res.L5.thorax.by.wingpadPC1 <- results(dds.L5.thorax.by.wingpadPC1, name=c("wingpadPC1"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.wingpadPC1 <- lfcShrink(dds.L5.thorax.by.wingpadPC1, coef="wingpadPC1", type="apeglm", res = res.L5.thorax.by.wingpadPC1)
summary(res.L5.thorax.by.wingpadPC1)
# out of 24958 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.004%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%

plotMA(res.L5.thorax.by.wingpadPC1, main="apeGLM")
plot(plotCounts(dds.L5.thorax.by.wingpadPC1, gene=order(res.L5.thorax.by.wingpadPC1$padj)[4], intgroup="wingpadPC1", returnData=TRUE))
res.L5.thorax.by.wingpadPC1$log2FoldChange[order(res.L5.thorax.by.wingpadPC1$padj)[4]]
# Positive LFC values are associated with positive wingpadPC1 values

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.wingpadPC1 <- vst(dds.L5.thorax.by.wingpadPC1, blind=FALSE)
plotPCA(vsd.L5.thorax.by.wingpadPC1, intgroup=c("wingpadPC1"))

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.wingpadPC1 <- res.L5.thorax.by.wingpadPC1[order(res.L5.thorax.by.wingpadPC1$padj),]
deg100.L5.thorax.by.wingpadPC1 <- deg100.L5.thorax.by.wingpadPC1[1:100,]
deg100.L5.thorax.by.wingpadPC1 <- apply.annotation(deg100.L5.thorax.by.wingpadPC1, ann)
deg100.L5.thorax.by.wingpadPC1 <- add.manual.annotations(deg100.L5.thorax.by.wingpadPC1)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.wingpadPC1), "model.output/deg100.L5.thorax.by.wingpadPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.thorax.by.wingpadPC1, res.L5.thorax.by.wingpadPC1, deg100.L5.thorax.by.wingpadPC1,
     file = "model.output/L5.thorax.by.wingpadPC1.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.wingpadPC1, res.L5.thorax.by.wingpadPC1, vsd.L5.thorax.by.wingpadPC1)
# load("L5.thorax.by.wingpadPC1.rda", verbose = TRUE)

###################
### L5 gonad by wingpad PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue!="thorax", !is.na(wingpadPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(food_regime,food_regime,length)))
hist(meta.i$wingpadPC1) # higher PC1 values are associated with more short-wing shapes

dds.L5.gonad.by.wingpadPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + wingpadPC1 )
dds.L5.gonad.by.wingpadPC1 <- DESeq(dds.L5.gonad.by.wingpadPC1)
dds.L5.gonad.by.wingpadPC1 <- estimateSizeFactors(dds.L5.gonad.by.wingpadPC1)

# # Extract and save normalized counts
# ctn.L5.gonad.by.wingpadPC1 <- counts(dds.L5.gonad.by.wingpadPC1, normalized=TRUE)
# write.csv(ctn.L5.gonad.by.wingpadPC1, 'normalized.counts.L5.gonad.by.wingpadPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.wingpadPC1)
res.L5.gonad.by.wingpadPC1 <- results(dds.L5.gonad.by.wingpadPC1, name=c("wingpadPC1"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.wingpadPC1 <- lfcShrink(dds.L5.gonad.by.wingpadPC1, coef="wingpadPC1", type="apeglm", res = res.L5.gonad.by.wingpadPC1)
summary(res.L5.gonad.by.wingpadPC1)
# out of 32603 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.0031%
# LFC < 0 (down)     : 5, 0.015%
# outliers [1]       : 0, 0%

plotMA(res.L5.gonad.by.wingpadPC1, main="apeGLM")
plot(plotCounts(dds.L5.gonad.by.wingpadPC1, gene=order(res.L5.gonad.by.wingpadPC1$padj)[4], intgroup="wingpadPC1", returnData=TRUE))
res.L5.gonad.by.wingpadPC1$log2FoldChange[order(res.L5.gonad.by.wingpadPC1$padj)[4]]
# Positive LFC values are associated with positive wingpadPC1 values

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.wingpadPC1 <- vst(dds.L5.gonad.by.wingpadPC1, blind=FALSE)
plotPCA(vsd.L5.gonad.by.wingpadPC1, intgroup=c("wingpadPC1"))

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.wingpadPC1 <- res.L5.gonad.by.wingpadPC1[order(res.L5.gonad.by.wingpadPC1$padj),]
deg100.L5.gonad.by.wingpadPC1 <- deg100.L5.gonad.by.wingpadPC1[1:100,]
deg100.L5.gonad.by.wingpadPC1 <- apply.annotation(deg100.L5.gonad.by.wingpadPC1, ann)
deg100.L5.gonad.by.wingpadPC1 <- add.manual.annotations(deg100.L5.gonad.by.wingpadPC1)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.wingpadPC1), "model.output/deg100.L5.gonad.by.wingpadPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.gonad.by.wingpadPC1, res.L5.gonad.by.wingpadPC1, deg100.L5.gonad.by.wingpadPC1,
     file = "model.output/L5.gonad.by.wingpadPC1.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.wingpadPC1, res.L5.gonad.by.wingpadPC1, vsd.L5.gonad.by.wingpadPC1)
# load("L5.gonad.by.wingpadPC1.rda", verbose = TRUE)

###################
### Adult thorax by thorax shape PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue=="thorax", !is.na(txPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))
hist(meta.i$txPC1) # higher PC1 values are associated with more short-wing shapes

dds.adult.thorax.by.txPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + txPC1 )
dds.adult.thorax.by.txPC1 <- DESeq(dds.adult.thorax.by.txPC1)
dds.adult.thorax.by.txPC1 <- estimateSizeFactors(dds.adult.thorax.by.txPC1)

# # Extract and save normalized counts
# ctn.adult.thorax.by.txPC1 <- counts(dds.adult.thorax.by.txPC1, normalized=TRUE)
# write.csv(ctn.adult.thorax.by.txPC1, 'normalized.counts.adult.thorax.by.txPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.txPC1)
res.adult.thorax.by.txPC1 <- results(dds.adult.thorax.by.txPC1, name=c("txPC1"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.txPC1 <- lfcShrink(dds.adult.thorax.by.txPC1, coef="txPC1", type="apeglm", res = res.adult.thorax.by.txPC1)
summary(res.adult.thorax.by.txPC1)
# out of 32236 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 169, 0.52%
# LFC < 0 (down)     : 37, 0.11%
# outliers [1]       : 0, 0%

plotMA(res.adult.thorax.by.txPC1, main="apeGLM", ylim=c(-30,30))
plot(plotCounts(dds.adult.thorax.by.txPC1, gene=order(res.adult.thorax.by.txPC1$padj)[1], intgroup="txPC1", returnData=TRUE))
res.adult.thorax.by.txPC1$log2FoldChange[order(res.adult.thorax.by.txPC1$padj)[1]]
# Positive LFC values are associated with positive txPC1 values (wider-thorax bias)

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.txPC1 <- vst(dds.adult.thorax.by.txPC1, blind=FALSE)
plotPCA(vsd.adult.thorax.by.txPC1, intgroup=c("txPC1"))
plotPCA(vsd.adult.thorax.by.txPC1, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.txPC1 <- res.adult.thorax.by.txPC1[order(res.adult.thorax.by.txPC1$padj),]
deg100.adult.thorax.by.txPC1 <- deg100.adult.thorax.by.txPC1[1:100,]
deg100.adult.thorax.by.txPC1 <- apply.annotation(deg100.adult.thorax.by.txPC1, ann)
deg100.adult.thorax.by.txPC1 <- add.manual.annotations(deg100.adult.thorax.by.txPC1)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.txPC1), "model.output/deg100.adult.thorax.by.txPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.thorax.by.txPC1, res.adult.thorax.by.txPC1, deg100.adult.thorax.by.txPC1,
     file = "model.output/adult.thorax.by.txPC1.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.txPC1, res.adult.thorax.by.txPC1, vsd.adult.thorax.by.txPC1)
# load("adult.thorax.by.txPC1.rda", verbose = TRUE)

###################
### Adult gonad by thorax shape PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue=="gonad", !is.na(txPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph,morph,length)))
hist(meta.i$txPC1) # higher PC1 values are associated with more short-wing shapes

dds.adult.gonad.by.txPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + txPC1 )
dds.adult.gonad.by.txPC1 <- DESeq(dds.adult.gonad.by.txPC1)
dds.adult.gonad.by.txPC1 <- estimateSizeFactors(dds.adult.gonad.by.txPC1)

# # Extract and save normalized counts
# ctn.adult.gonad.by.txPC1 <- counts(dds.adult.gonad.by.txPC1, normalized=TRUE)
# write.csv(ctn.adult.gonad.by.txPC1, 'normalized.counts.adult.gonad.by.txPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.txPC1)
res.adult.gonad.by.txPC1 <- results(dds.adult.gonad.by.txPC1, name=c("txPC1"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.txPC1 <- lfcShrink(dds.adult.gonad.by.txPC1, coef="txPC1", type="apeglm", res = res.adult.gonad.by.txPC1)
summary(res.adult.gonad.by.txPC1)
# out of 54574 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%

plotMA(res.adult.gonad.by.txPC1, main="apeGLM", ylim=c(-30,30))
plot(plotCounts(dds.adult.gonad.by.txPC1, gene=order(res.adult.gonad.by.txPC1$padj)[1], intgroup="txPC1", returnData=TRUE))
res.adult.gonad.by.txPC1$log2FoldChange[order(res.adult.gonad.by.txPC1$padj)[1]]
# Positive LFC values are associated with positive txPC1 values (wider-gonad bias)

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.txPC1 <- vst(dds.adult.gonad.by.txPC1, blind=FALSE)
plotPCA(vsd.adult.gonad.by.txPC1, intgroup=c("txPC1"))
plotPCA(vsd.adult.gonad.by.txPC1, intgroup=c("morph"))

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.txPC1 <- res.adult.gonad.by.txPC1[order(res.adult.gonad.by.txPC1$padj),]
deg100.adult.gonad.by.txPC1 <- deg100.adult.gonad.by.txPC1[1:100,]
deg100.adult.gonad.by.txPC1 <- apply.annotation(deg100.adult.gonad.by.txPC1, ann)
deg100.adult.gonad.by.txPC1 <- add.manual.annotations(deg100.adult.gonad.by.txPC1)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.txPC1), "model.output/deg100.adult.gonad.by.txPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.gonad.by.txPC1, res.adult.gonad.by.txPC1, deg100.adult.gonad.by.txPC1,
     file = "model.output/adult.gonad.by.txPC1.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.txPC1, res.adult.gonad.by.txPC1, vsd.adult.gonad.by.txPC1)
# load("adult.gonad.by.txPC1.rda", verbose = TRUE)

###################
### L5 thorax by thorax shape PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue=="thorax", !is.na(txPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
hist(meta.i$txPC1) # higher PC1 values are associated with more short-wing shapes

dds.L5.thorax.by.txPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + txPC1 )
dds.L5.thorax.by.txPC1 <- DESeq(dds.L5.thorax.by.txPC1)
dds.L5.thorax.by.txPC1 <- estimateSizeFactors(dds.L5.thorax.by.txPC1)

# # Extract and save normalized counts
# ctn.L5.thorax.by.txPC1 <- counts(dds.L5.thorax.by.txPC1, normalized=TRUE)
# write.csv(ctn.L5.thorax.by.txPC1, 'normalized.counts.L5.thorax.by.txPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.txPC1)
res.L5.thorax.by.txPC1 <- results(dds.L5.thorax.by.txPC1, name=c("txPC1"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.txPC1 <- lfcShrink(dds.L5.thorax.by.txPC1, coef="txPC1", type="apeglm", res = res.L5.thorax.by.txPC1)
summary(res.L5.thorax.by.txPC1)
# out of 24549 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%

plotMA(res.L5.thorax.by.txPC1, main="apeGLM")
plot(plotCounts(dds.L5.thorax.by.txPC1, gene=order(res.L5.thorax.by.txPC1$padj)[1], intgroup="txPC1", returnData=TRUE))
res.L5.thorax.by.txPC1$log2FoldChange[order(res.L5.thorax.by.txPC1$padj)[1]]
# Positive LFC values are associated with positive txPC1 values (narrower-mesonotum bias)

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.txPC1 <- vst(dds.L5.thorax.by.txPC1, blind=FALSE)
plotPCA(vsd.L5.thorax.by.txPC1, intgroup=c("txPC1"))
plotPCA(vsd.L5.thorax.by.txPC1, intgroup=c("food_regime"))

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.txPC1 <- res.L5.thorax.by.txPC1[order(res.L5.thorax.by.txPC1$padj),]
deg100.L5.thorax.by.txPC1 <- deg100.L5.thorax.by.txPC1[1:100,]
deg100.L5.thorax.by.txPC1 <- apply.annotation(deg100.L5.thorax.by.txPC1, ann)
deg100.L5.thorax.by.txPC1 <- add.manual.annotations(deg100.L5.thorax.by.txPC1)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.txPC1), "model.output/deg100.L5.thorax.by.txPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.thorax.by.txPC1, res.L5.thorax.by.txPC1, deg100.L5.thorax.by.txPC1,
     file = "model.output/L5.thorax.by.txPC1.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.txPC1, res.L5.thorax.by.txPC1, vsd.L5.thorax.by.txPC1)
# load("L5.thorax.by.txPC1.rda", verbose = TRUE)

###################
### L5 gonad by thorax shape PC1
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue!="thorax", !is.na(txPC1))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
hist(meta.i$txPC1) # higher PC1 values are associated with more short-wing shapes

dds.L5.gonad.by.txPC1 <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + txPC1 )
dds.L5.gonad.by.txPC1 <- DESeq(dds.L5.gonad.by.txPC1)
dds.L5.gonad.by.txPC1 <- estimateSizeFactors(dds.L5.gonad.by.txPC1)

# # Extract and save normalized counts
# ctn.L5.gonad.by.txPC1 <- counts(dds.L5.gonad.by.txPC1, normalized=TRUE)
# write.csv(ctn.L5.gonad.by.txPC1, 'normalized.counts.L5.gonad.by.txPC1.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.txPC1)
res.L5.gonad.by.txPC1 <- results(dds.L5.gonad.by.txPC1, name=c("txPC1"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.txPC1 <- lfcShrink(dds.L5.gonad.by.txPC1, coef="txPC1", type="apeglm", res = res.L5.gonad.by.txPC1)
summary(res.L5.gonad.by.txPC1)
# out of 32155 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%

plotMA(res.L5.gonad.by.txPC1, main="apeGLM")
plot(plotCounts(dds.L5.gonad.by.txPC1, gene=order(res.L5.gonad.by.txPC1$padj)[1], intgroup="txPC1", returnData=TRUE))
res.L5.gonad.by.txPC1$log2FoldChange[order(res.L5.gonad.by.txPC1$padj)[1]]
# Positive LFC values are associated with positive txPC1 values (narrower-mesonotum bias)

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.txPC1 <- vst(dds.L5.gonad.by.txPC1, blind=FALSE)
plotPCA(vsd.L5.gonad.by.txPC1, intgroup=c("txPC1"))
plotPCA(vsd.L5.gonad.by.txPC1, intgroup=c("food_regime"))

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.txPC1 <- res.L5.gonad.by.txPC1[order(res.L5.gonad.by.txPC1$padj),]
deg100.L5.gonad.by.txPC1 <- deg100.L5.gonad.by.txPC1[1:100,]
deg100.L5.gonad.by.txPC1 <- apply.annotation(deg100.L5.gonad.by.txPC1, ann)
deg100.L5.gonad.by.txPC1 <- add.manual.annotations(deg100.L5.gonad.by.txPC1)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.txPC1), "model.output/deg100.L5.gonad.by.txPC1.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.gonad.by.txPC1, res.L5.gonad.by.txPC1, deg100.L5.gonad.by.txPC1,
     file = "model.output/L5.gonad.by.txPC1.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.txPC1, res.L5.gonad.by.txPC1, vsd.L5.gonad.by.txPC1)
# load("L5.gonad.by.txPC1.rda", verbose = TRUE)

###################
### Adult thorax by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)

dds.adult.thorax.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph + food_density )
dds.adult.thorax.by.food_density <- DESeq(dds.adult.thorax.by.food_density)
dds.adult.thorax.by.food_density <- estimateSizeFactors(dds.adult.thorax.by.food_density)

# # Extract and save normalized counts
# ctn.adult.thorax.by.food_density <- counts(dds.adult.thorax.by.food_density, normalized=TRUE)
# write.csv(ctn.adult.thorax.by.food_density, 'normalized.counts.adult.thorax.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.food_density)
res.adult.thorax.by.food_density <- results(dds.adult.thorax.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.food_density <- lfcShrink(dds.adult.thorax.by.food_density, coef="food_density", type="apeglm", res = res.adult.thorax.by.food_density)
summary(res.adult.thorax.by.food_density)
# out of 32368 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 0.028%
# LFC < 0 (down)     : 13, 0.04%
# outliers [1]       : 145, 0.45%

plotMA(res.adult.thorax.by.food_density, main="apeGLM", ylim=c(-4,5))
plot(plotCounts(dds.adult.thorax.by.food_density, gene=order(res.adult.thorax.by.food_density$padj)[6], intgroup="food_density", returnData=TRUE))
res.adult.thorax.by.food_density$log2FoldChange[order(res.adult.thorax.by.food_density$padj)[6]]
# Positive LFC values are associated with higher food density expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.food_density <- vst(dds.adult.thorax.by.food_density, blind=FALSE)
plotPCA(vsd.adult.thorax.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.food_density <- res.adult.thorax.by.food_density[order(res.adult.thorax.by.food_density$padj),]
deg100.adult.thorax.by.food_density <- deg100.adult.thorax.by.food_density[1:100,]
deg100.adult.thorax.by.food_density <- apply.annotation(deg100.adult.thorax.by.food_density, ann)
deg100.adult.thorax.by.food_density <- add.manual.annotations(deg100.adult.thorax.by.food_density)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.food_density), "model.output/deg100.adult.thorax.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.thorax.by.food_density, res.adult.thorax.by.food_density, deg100.adult.thorax.by.food_density,
     file = "model.output/adult.thorax.by.food_density.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.food_density, res.adult.thorax.by.food_density, vsd.adult.thorax.by.food_density)
# load("ovaries.by.food_density.rda", verbose = TRUE)

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
deg100.adult.thorax.by.food <- add.manual.annotations(deg100.adult.thorax.by.food)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.food), "model.output/deg100.adult.thorax.by.food.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.thorax.by.food, res.adult.thorax.by.food, deg100.adult.thorax.by.food,
     file = "model.output/adult.thorax.by.food_regime.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.food, res.adult.thorax.by.food, vsd.adult.thorax.by.food)
# load("adult.thorax.by.food_regime.rda", verbose = TRUE)

###################
### Adult gonad by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)

dds.adult.gonad.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph + food_density )
dds.adult.gonad.by.food_density <- DESeq(dds.adult.gonad.by.food_density)
dds.adult.gonad.by.food_density <- estimateSizeFactors(dds.adult.gonad.by.food_density)

# # Extract and save normalized counts
# ctn.adult.gonad.by.food_density <- counts(dds.adult.gonad.by.food_density, normalized=TRUE)
# write.csv(ctn.adult.gonad.by.food_density, 'normalized.counts.adult.gonad.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.food_density)
res.adult.gonad.by.food_density <- results(dds.adult.gonad.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.food_density <- lfcShrink(dds.adult.gonad.by.food_density, coef="food_density", type="apeglm", res = res.adult.gonad.by.food_density)
summary(res.adult.gonad.by.food_density)
# out of 54952 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 0.016%
# LFC < 0 (down)     : 49, 0.089%
# outliers [1]       : 165, 0.3%

plotMA(res.adult.gonad.by.food_density, main="apeGLM", ylim=c(-8,3))
plot(plotCounts(dds.adult.gonad.by.food_density, gene=order(res.adult.gonad.by.food_density$padj)[1], intgroup="food_density", returnData=TRUE))
res.adult.gonad.by.food_density$log2FoldChange[order(res.adult.gonad.by.food_density$padj)[6]]
# Positive LFC values are associated with higher food density expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.food_density <- vst(dds.adult.gonad.by.food_density, blind=FALSE)
plotPCA(vsd.adult.gonad.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.food_density <- res.adult.gonad.by.food_density[order(res.adult.gonad.by.food_density$padj),]
deg100.adult.gonad.by.food_density <- deg100.adult.gonad.by.food_density[1:100,]
deg100.adult.gonad.by.food_density <- apply.annotation(deg100.adult.gonad.by.food_density, ann)
deg100.adult.gonad.by.food_density <- add.manual.annotations(deg100.adult.gonad.by.food_density)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.food_density), "model.output/deg100.adult.gonad.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.gonad.by.food_density, res.adult.gonad.by.food_density, deg100.adult.gonad.by.food_density,
     file = "model.output/adult.gonad.by.food_density.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.food_density, res.adult.gonad.by.food_density, vsd.adult.gonad.by.food_density)
# load("adult.gonad.by.food_density.rda", verbose = TRUE)

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
deg100.adult.gonad.by.food <- add.manual.annotations(deg100.adult.gonad.by.food)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.food), "model.output/deg100.adult.gonad.by.food.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.gonad.by.food, res.adult.gonad.by.food, deg100.adult.gonad.by.food,
     file = "model.output/adult.gonad.by.food_regime.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.food, res.adult.gonad.by.food, vsd.adult.gonad.by.food)
# load("adult.gonad.by.food_regime.rda", verbose = TRUE)

###################
### Adult ovaries by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", sex == "f", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)

dds.adult.ovaries.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph + food_density )
dds.adult.ovaries.by.food_density <- DESeq(dds.adult.ovaries.by.food_density)
dds.adult.ovaries.by.food_density <- estimateSizeFactors(dds.adult.ovaries.by.food_density)

# # Extract and save normalized counts
# ctn.adult.ovaries.by.food_density <- counts(dds.adult.ovaries.by.food_density, normalized=TRUE)
# write.csv(ctn.adult.ovaries.by.food_density, 'normalized.counts.adult.ovaries.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.ovaries.by.food_density)
res.adult.ovaries.by.food_density <- results(dds.adult.ovaries.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.adult.ovaries.by.food_density <- lfcShrink(dds.adult.ovaries.by.food_density, coef="food_density", type="apeglm", res = res.adult.ovaries.by.food_density)
summary(res.adult.ovaries.by.food_density)
# out of 26681 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 21, 0.079%
# LFC < 0 (down)     : 16, 0.06%
# outliers [1]       : 0, 0%

plotMA(res.adult.ovaries.by.food_density, main="apeGLM", ylim=c(-8,3))
plot(plotCounts(dds.adult.ovaries.by.food_density, gene=order(res.adult.ovaries.by.food_density$padj)[1], intgroup="food_density", returnData=TRUE))
res.adult.ovaries.by.food_density$log2FoldChange[order(res.adult.ovaries.by.food_density$padj)[6]]
# Positive LFC values are associated with higher food density expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.ovaries.by.food_density <- vst(dds.adult.ovaries.by.food_density, blind=FALSE)
plotPCA(vsd.adult.ovaries.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.adult.ovaries.by.food_density <- res.adult.ovaries.by.food_density[order(res.adult.ovaries.by.food_density$padj),]
deg100.adult.ovaries.by.food_density <- deg100.adult.ovaries.by.food_density[1:100,]
deg100.adult.ovaries.by.food_density <- apply.annotation(deg100.adult.ovaries.by.food_density, ann)
deg100.adult.ovaries.by.food_density <- add.manual.annotations(deg100.adult.ovaries.by.food_density)

# Save results
write.table(as.data.frame(deg100.adult.ovaries.by.food_density), "model.output/deg100.adult.ovaries.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.ovaries.by.food_density, res.adult.ovaries.by.food_density, deg100.adult.ovaries.by.food_density,
     file = "model.output/adult.ovaries.by.food_density.rda")
rm(meta.i, cts.i, dds.adult.ovaries.by.food_density, res.adult.ovaries.by.food_density, vsd.adult.ovaries.by.food_density)
# load("ovaries.by.food_density.rda", verbose = TRUE)

###################
### Adult testes by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", sex == "m", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)
hist(meta.i$food_density)

dds.adult.testes.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph + food_density )
dds.adult.testes.by.food_density <- DESeq(dds.adult.testes.by.food_density)
dds.adult.testes.by.food_density <- estimateSizeFactors(dds.adult.testes.by.food_density)

# # Extract and save normalized counts
# ctn.adult.testes.by.food_density <- counts(dds.adult.testes.by.food_density, normalized=TRUE)
# write.csv(ctn.adult.testes.by.food_density, 'normalized.counts.adult.testes.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.testes.by.food_density)
res.adult.testes.by.food_density <- results(dds.adult.testes.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.adult.testes.by.food_density <- lfcShrink(dds.adult.testes.by.food_density, coef="food_density", type="apeglm", res = res.adult.testes.by.food_density)
summary(res.adult.testes.by.food_density)
# out of 41085 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 0.022%
# LFC < 0 (down)     : 4, 0.0097%
# outliers [1]       : 350, 0.85%

plotMA(res.adult.testes.by.food_density, main="apeGLM", ylim=c(-5,1.5))
plot(plotCounts(dds.adult.testes.by.food_density, gene=order(res.adult.testes.by.food_density$padj)[1], intgroup="food_density", returnData=TRUE))
res.adult.testes.by.food_density$log2FoldChange[order(res.adult.testes.by.food_density$padj)[6]]
# Positive LFC values are associated with higher food density expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.testes.by.food_density <- vst(dds.adult.testes.by.food_density, blind=FALSE)
plotPCA(vsd.adult.testes.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.adult.testes.by.food_density <- res.adult.testes.by.food_density[order(res.adult.testes.by.food_density$padj),]
deg100.adult.testes.by.food_density <- deg100.adult.testes.by.food_density[1:100,]
deg100.adult.testes.by.food_density <- apply.annotation(deg100.adult.testes.by.food_density, ann)
deg100.adult.testes.by.food_density <- add.manual.annotations(deg100.adult.testes.by.food_density)

# Save results
write.table(as.data.frame(deg100.adult.testes.by.food_density), "model.output/deg100.adult.testes.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.testes.by.food_density, res.adult.testes.by.food_density, deg100.adult.testes.by.food_density,
     file = "model.output/adult.testes.by.food_density.rda")
rm(meta.i, cts.i, dds.adult.testes.by.food_density, res.adult.testes.by.food_density, vsd.adult.testes.by.food_density)
# load("adult.testes.by.food_density.rda", verbose = TRUE)

###################
### L5 thorax by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)
hist(meta.i$food_density)

dds.L5.thorax.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + food_density )
dds.L5.thorax.by.food_density <- DESeq(dds.L5.thorax.by.food_density)
dds.L5.thorax.by.food_density <- estimateSizeFactors(dds.L5.thorax.by.food_density)

# # Extract and save normalized counts
# ctn.L5.thorax.by.food_density <- counts(dds.L5.thorax.by.food_density, normalized=TRUE)
# write.csv(ctn.L5.thorax.by.food_density, 'normalized.counts.L5.thorax.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.food_density)
res.L5.thorax.by.food_density <- results(dds.L5.thorax.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.food_density <- lfcShrink(dds.L5.thorax.by.food_density, coef="food_density", type="apeglm", res = res.L5.thorax.by.food_density)
summary(res.L5.thorax.by.food_density)
# out of 24958 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 10, 0.04%
# LFC < 0 (down)     : 6, 0.024%
# outliers [1]       : 229, 0.92%

plotMA(res.L5.thorax.by.food_density, main="apeGLM", ylim=c(-4,2))
plot(plotCounts(dds.L5.thorax.by.food_density, gene=order(res.L5.thorax.by.food_density$padj)[1], intgroup="food_density", returnData=TRUE))
res.L5.thorax.by.food_density$log2FoldChange[order(res.L5.thorax.by.food_density$padj)[1]]
# Positive LFC values are associated with higher food density expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.food_density <- vst(dds.L5.thorax.by.food_density, blind=FALSE)
plotPCA(vsd.L5.thorax.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.food_density <- res.L5.thorax.by.food_density[order(res.L5.thorax.by.food_density$padj),]
deg100.L5.thorax.by.food_density <- deg100.L5.thorax.by.food_density[1:100,]
deg100.L5.thorax.by.food_density <- apply.annotation(deg100.L5.thorax.by.food_density, ann)
deg100.L5.thorax.by.food_density <- add.manual.annotations(deg100.L5.thorax.by.food_density)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.food_density), "model.output/deg100.L5.thorax.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.thorax.by.food_density, res.L5.thorax.by.food_density, deg100.L5.thorax.by.food_density,
     file = "model.output/L5.thorax.by.food_density.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.food_density, res.L5.thorax.by.food_density, vsd.L5.thorax.by.food_density)
# load("L5.thorax.by.food_density.rda", verbose = TRUE)

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
deg100.L5.gonad.by.food <- add.manual.annotations(deg100.L5.gonad.by.food)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.food), "model.output/deg100.L5.gonad.by.food.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.gonad.by.food, res.L5.gonad.by.food, deg100.L5.gonad.by.food,
     file = "model.output/L5.gonad.by.food_regime.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.food, res.L5.gonad.by.food, vsd.L5.gonad.by.food)
# load("L5.gonad.by.food_regime.rda", verbose = TRUE)

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
summary(res.L5.thorax.by.food)
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
deg100.L5.thorax.by.food <- add.manual.annotations(deg100.L5.thorax.by.food)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.food), "model.output/deg100.L5.thorax.by.food.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.thorax.by.food, res.L5.thorax.by.food, deg100.L5.thorax.by.food,
     file = "model.output/L5.thorax.by.food_regime.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.food, res.L5.thorax.by.food, vsd.L5.thorax.by.food)
# load("L5.thorax.by.food.rda", verbose = TRUE)

###################
### L5 gonad by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)
hist(meta.i$food_density)

dds.L5.gonad.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + food_density )
dds.L5.gonad.by.food_density <- DESeq(dds.L5.gonad.by.food_density)
dds.L5.gonad.by.food_density <- estimateSizeFactors(dds.L5.gonad.by.food_density)

# # Extract and save normalized counts
# ctn.L5.gonad.by.food_density <- counts(dds.L5.gonad.by.food_density, normalized=TRUE)
# write.csv(ctn.L5.gonad.by.food_density, 'normalized.counts.L5.gonad.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.food_density)
res.L5.gonad.by.food_density <- results(dds.L5.gonad.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.food_density <- lfcShrink(dds.L5.gonad.by.food_density, coef="food_density", type="apeglm", res = res.L5.gonad.by.food_density)
summary(res.L5.gonad.by.food_density)
# out of 32603 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4, 0.012%
# LFC < 0 (down)     : 3, 0.0092%
# outliers [1]       : 188, 0.58%

plotMA(res.L5.gonad.by.food_density, main="apeGLM", ylim=c(-3.5,3))
plot(plotCounts(dds.L5.gonad.by.food_density, gene=order(res.L5.gonad.by.food_density$padj)[1], intgroup="food_density", returnData=TRUE))
res.L5.gonad.by.food_density$log2FoldChange[order(res.L5.gonad.by.food_density$padj)[1]]
# Positive LFC values are associated with higher food density expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.food_density <- vst(dds.L5.gonad.by.food_density, blind=FALSE)
plotPCA(vsd.L5.gonad.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.food_density <- res.L5.gonad.by.food_density[order(res.L5.gonad.by.food_density$padj),]
deg100.L5.gonad.by.food_density <- deg100.L5.gonad.by.food_density[1:100,]
deg100.L5.gonad.by.food_density <- apply.annotation(deg100.L5.gonad.by.food_density, ann)
deg100.L5.gonad.by.food_density <- add.manual.annotations(deg100.L5.gonad.by.food_density)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.food_density), "model.output/deg100.L5.gonad.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.gonad.by.food_density, res.L5.gonad.by.food_density, deg100.L5.gonad.by.food_density,
     file = "model.output/L5.gonad.by.food_density.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.food_density, res.L5.gonad.by.food_density, vsd.L5.gonad.by.food_density)
# load("gonad.by.food_density.rda", verbose = TRUE)

###################
### L5 ovaries by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue != "thorax" & sex == "f")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(food_regime,food_regime,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)

dds.L5.ovaries.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + food_density )
dds.L5.ovaries.by.food_density <- DESeq(dds.L5.ovaries.by.food_density)
dds.L5.ovaries.by.food_density <- estimateSizeFactors(dds.L5.ovaries.by.food_density)

# # Extract and save normalized counts
# ctn.L5.ovaries.by.food_density <- counts(dds.L5.ovaries.by.food_density, normalized=TRUE)
# write.csv(ctn.L5.ovaries.by.food_density, 'normalized.counts.L5.ovaries.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.ovaries.by.food_density)
res.L5.ovaries.by.food_density <- results(dds.L5.ovaries.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.L5.ovaries.by.food_density <- lfcShrink(dds.L5.ovaries.by.food_density, coef="food_density", type="apeglm", res = res.L5.ovaries.by.food_density)
summary(res.L5.ovaries.by.food_density)
# out of 17635 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 0.051%
# LFC < 0 (down)     : 7, 0.04%
# outliers [1]       : 68, 0.39%

plotMA(res.L5.ovaries.by.food_density, main="apeGLM")
plotMA(res.L5.ovaries.by.food_density, main="apeGLM", ylim=c(-3.5,2.5))
plot(plotCounts(dds.L5.ovaries.by.food_density, gene=order(res.L5.ovaries.by.food_density$padj)[8], intgroup="food_density", returnData=TRUE))
res.L5.ovaries.by.food_density$log2FoldChange[order(res.L5.ovaries.by.food_density$padj)[8]]
# Positive LFC values are associated with low food expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.ovaries.by.food_density <- vst(dds.L5.ovaries.by.food_density, blind=FALSE)
plotPCA(vsd.L5.ovaries.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.L5.ovaries.by.food_density <- res.L5.ovaries.by.food_density[order(res.L5.ovaries.by.food_density$padj),]
deg100.L5.ovaries.by.food_density <- deg100.L5.ovaries.by.food_density[1:100,]
deg100.L5.ovaries.by.food_density <- apply.annotation(deg100.L5.ovaries.by.food_density, ann)
deg100.L5.ovaries.by.food_density <- add.manual.annotations(deg100.L5.ovaries.by.food_density)

# Save results
write.table(as.data.frame(deg100.L5.ovaries.by.food_density), "model.output/deg100.L5.ovaries.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.ovaries.by.food_density, res.L5.ovaries.by.food_density, deg100.L5.ovaries.by.food_density,
     file = "model.output/L5.ovaries.by.food_density.rda")
rm(meta.i, cts.i, dds.L5.ovaries.by.food_density, res.L5.ovaries.by.food_density, vsd.L5.ovaries.by.food_density)
# load("L5.ovaries.by.food_density.rda", verbose = TRUE)


###################
### L5 ovaries by food regime
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue != "thorax" & sex == "f")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(food_regime,food_regime,length)))

dds.L5.ovaries.by.food <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + food_regime )
dds.L5.ovaries.by.food <- DESeq(dds.L5.ovaries.by.food)
dds.L5.ovaries.by.food <- estimateSizeFactors(dds.L5.ovaries.by.food)

# # Extract and save normalized counts
# ctn.L5.ovaries.by.food <- counts(dds.L5.ovaries.by.food, normalized=TRUE)
# write.csv(ctn.L5.ovaries.by.food, 'normalized.counts.L5.ovaries.by.food.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.ovaries.by.food)
res.L5.ovaries.by.food <- results(dds.L5.ovaries.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.L5.ovaries.by.food <- lfcShrink(dds.L5.ovaries.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.L5.ovaries.by.food)
summary(res.L5.ovaries.by.food)
# out of 17633 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4, 0.023%
# LFC < 0 (down)     : 4, 0.023%
# outliers [1]       : 117, 0.66%

plotMA(res.L5.ovaries.by.food, main="apeGLM", ylim=c(-4,8.5))
plotCounts(dds.L5.ovaries.by.food, gene=order(res.L5.ovaries.by.food$padj)[3], intgroup="food_regime")
res.L5.ovaries.by.food$log2FoldChange[order(res.L5.ovaries.by.food$padj)[3]]
# Positive LFC values are associated with low food expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.ovaries.by.food <- vst(dds.L5.ovaries.by.food, blind=FALSE)
plotPCA(vsd.L5.ovaries.by.food, intgroup=c("food_regime"))

# We can order our results table by the smallest p value:
deg100.L5.ovaries.by.food <- res.L5.ovaries.by.food[order(res.L5.ovaries.by.food$padj),]
deg100.L5.ovaries.by.food <- deg100.L5.ovaries.by.food[1:100,]
deg100.L5.ovaries.by.food <- apply.annotation(deg100.L5.ovaries.by.food, ann)
deg100.L5.ovaries.by.food <- add.manual.annotations(deg100.L5.ovaries.by.food)

# Save results
write.table(as.data.frame(deg100.L5.ovaries.by.food), "model.output/deg100.L5.ovaries.by.food.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.ovaries.by.food, res.L5.ovaries.by.food, deg100.L5.ovaries.by.food,
     file = "model.output/L5.ovaries.by.food_regime.rda")
rm(meta.i, cts.i, dds.L5.ovaries.by.food, res.L5.ovaries.by.food, vsd.L5.ovaries.by.food)
# load("L5.ovaries.by.food_regime.rda", verbose = TRUE)

###################
### L5 testes by food density
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", sex == "m", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

meta.i$food_density <- log10(meta.i$seeds/meta.i$cohort)
hist(meta.i$food_density)

dds.L5.testes.by.food_density <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + food_density )
dds.L5.testes.by.food_density <- DESeq(dds.L5.testes.by.food_density)
dds.L5.testes.by.food_density <- estimateSizeFactors(dds.L5.testes.by.food_density)

# # Extract and save normalized counts
# ctn.L5.testes.by.food_density <- counts(dds.L5.testes.by.food_density, normalized=TRUE)
# write.csv(ctn.L5.testes.by.food_density, 'normalized.counts.L5.testes.by.food_density.csv')
# # Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.testes.by.food_density)
res.L5.testes.by.food_density <- results(dds.L5.testes.by.food_density, name=c("food_density"), alpha=0.05, filterFun=ihw)
res.L5.testes.by.food_density <- lfcShrink(dds.L5.testes.by.food_density, coef="food_density", type="apeglm", res = res.L5.testes.by.food_density)
summary(res.L5.testes.by.food_density)
# out of 24508 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 29, 0.12%
# LFC < 0 (down)     : 20, 0.082%
# outliers [1]       : 278, 1.1%

plotMA(res.L5.testes.by.food_density, main="apeGLM", ylim=c(-7,5))
plot(plotCounts(dds.L5.testes.by.food_density, gene=order(res.L5.testes.by.food_density$padj)[1], intgroup="food_density", returnData=TRUE))
res.L5.testes.by.food_density$log2FoldChange[order(res.L5.testes.by.food_density$padj)[1]]
# Positive LFC values are associated with higher food density expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.testes.by.food_density <- vst(dds.L5.testes.by.food_density, blind=FALSE)
plotPCA(vsd.L5.testes.by.food_density, intgroup=c("food_density"))

# We can order our results table by the smallest p value:
deg100.L5.testes.by.food_density <- res.L5.testes.by.food_density[order(res.L5.testes.by.food_density$padj),]
deg100.L5.testes.by.food_density <- deg100.L5.testes.by.food_density[1:100,]
deg100.L5.testes.by.food_density <- apply.annotation(deg100.L5.testes.by.food_density, ann)
deg100.L5.testes.by.food_density <- add.manual.annotations(deg100.L5.testes.by.food_density)

# Save results
write.table(as.data.frame(deg100.L5.testes.by.food_density), "model.output/deg100.L5.testes.by.food_density.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.testes.by.food_density, res.L5.testes.by.food_density, deg100.L5.testes.by.food_density,
     file = "model.output/L5.testes.by.food_density.rda")
rm(meta.i, cts.i, dds.L5.testes.by.food_density, res.L5.testes.by.food_density, vsd.L5.testes.by.food_density)
# load("L5.testes.by.food_density.rda", verbose = TRUE)

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
deg100.thorax.by.stage <- add.manual.annotations(deg100.thorax.by.stage)

# Save results
write.table(as.data.frame(deg100.thorax.by.stage), "model.output/deg100.thorax.by.stage.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.thorax.by.stage, res.thorax.by.stage, deg100.thorax.by.stage,
     file = "model.output/thorax.by.stage.rda")
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
deg100.ovaries.by.stage <- add.manual.annotations(deg100.ovaries.by.stage)

# Save results
write.table(as.data.frame(deg100.ovaries.by.stage), "model.output/deg100.ovaries.by.stage.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.ovaries.by.stage, res.ovaries.by.stage, deg100.ovaries.by.stage,
     file = "model.output/ovaries.by.stage.rda")
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
deg100.testes.by.stage <- add.manual.annotations(deg100.testes.by.stage)

# Save results
write.table(as.data.frame(deg100.testes.by.stage), "model.output/deg100.testes.by.stage.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.testes.by.stage, res.testes.by.stage, deg100.testes.by.stage,
     file = "model.output/testes.by.stage.rda")
rm(meta.i, cts.i, dds.testes.by.stage, res.testes.by.stage, vsd.testes.by.stage)
# load("testes.by.stage.rda", verbose = TRUE) 

###################
# Save the results
###################
save(ann, apply.annotation, gene.read.number.cut.off, 
     dds.adult.gonad.by.food,dds.adult.gonad.by.food_density,dds.adult.gonad.by.girth,dds.adult.gonad.by.morph,dds.adult.gonad.by.sex,dds.adult.gonad.by.txPC1,dds.adult.gonad.by.wingPC1,dds.adult.ovaries.by.food_density,dds.adult.ovaries.by.girth,dds.adult.ovaries.by.morph,dds.adult.ovaries.by.wingPC1,dds.adult.testes.by.food_density,dds.adult.testes.by.girth,dds.adult.testes.by.morph,dds.adult.testes.by.wingPC1,dds.adult.thorax.by.food,dds.adult.thorax.by.food_density,dds.adult.thorax.by.girth,dds.adult.thorax.by.morph,dds.adult.thorax.by.sex,dds.adult.thorax.by.txPC1,dds.adult.thorax.by.wingPC1,dds.L5.gonad.by.food,dds.L5.gonad.by.food_density,dds.L5.gonad.by.girth,dds.L5.gonad.by.sex,dds.L5.gonad.by.txPC1,dds.L5.gonad.by.wingpadPC1,dds.L5.ovaries.by.food,dds.L5.ovaries.by.food_density,dds.L5.ovaries.by.girth,dds.L5.testes.by.food_density,dds.L5.testes.by.girth,dds.L5.thorax.by.food,dds.L5.thorax.by.food_density,dds.L5.thorax.by.girth,dds.L5.thorax.by.sex,dds.L5.thorax.by.txPC1,dds.L5.thorax.by.wingpadPC1,dds.ovaries.by.stage,dds.testes.by.stage,dds.thorax.by.stage,
     file = "model.output/dge.analysis.dds.rda")
rm(dds.adult.gonad.by.food,dds.adult.gonad.by.food_density,dds.adult.gonad.by.girth,dds.adult.gonad.by.morph,dds.adult.gonad.by.sex,dds.adult.gonad.by.txPC1,dds.adult.gonad.by.wingPC1,dds.adult.ovaries.by.food_density,dds.adult.ovaries.by.girth,dds.adult.ovaries.by.morph,dds.adult.ovaries.by.wingPC1,dds.adult.testes.by.food_density,dds.adult.testes.by.girth,dds.adult.testes.by.morph,dds.adult.testes.by.wingPC1,dds.adult.thorax.by.food,dds.adult.thorax.by.food_density,dds.adult.thorax.by.girth,dds.adult.thorax.by.morph,dds.adult.thorax.by.sex,dds.adult.thorax.by.txPC1,dds.adult.thorax.by.wingPC1,dds.L5.gonad.by.food,dds.L5.gonad.by.food_density,dds.L5.gonad.by.girth,dds.L5.gonad.by.sex,dds.L5.gonad.by.txPC1,dds.L5.gonad.by.wingpadPC1,dds.L5.ovaries.by.food,dds.L5.ovaries.by.food_density,dds.L5.ovaries.by.girth,dds.L5.testes.by.food_density,dds.L5.testes.by.girth,dds.L5.thorax.by.food,dds.L5.thorax.by.food_density,dds.L5.thorax.by.girth,dds.L5.thorax.by.sex,dds.L5.thorax.by.txPC1,dds.L5.thorax.by.wingpadPC1,dds.ovaries.by.stage,dds.testes.by.stage,dds.thorax.by.stage)
save(ann, apply.annotation, gene.read.number.cut.off, 
     res.adult.gonad.by.food,res.adult.gonad.by.food_density,res.adult.gonad.by.girth,res.adult.gonad.by.morph,res.adult.gonad.by.sex,res.adult.gonad.by.txPC1,res.adult.gonad.by.wingPC1,res.adult.ovaries.by.food_density,res.adult.ovaries.by.girth,res.adult.ovaries.by.morph,res.adult.ovaries.by.wingPC1,res.adult.testes.by.food_density,res.adult.testes.by.girth,res.adult.testes.by.morph,res.adult.testes.by.wingPC1,res.adult.thorax.by.food,res.adult.thorax.by.food_density,res.adult.thorax.by.girth,res.adult.thorax.by.morph,res.adult.thorax.by.sex,res.adult.thorax.by.txPC1,res.adult.thorax.by.wingPC1,res.L5.gonad.by.food,res.L5.gonad.by.food_density,res.L5.gonad.by.girth,res.L5.gonad.by.sex,res.L5.gonad.by.txPC1,res.L5.gonad.by.wingpadPC1,res.L5.ovaries.by.food,res.L5.ovaries.by.food_density,res.L5.ovaries.by.girth,res.L5.testes.by.food_density,res.L5.testes.by.girth,res.L5.thorax.by.food,res.L5.thorax.by.food_density,res.L5.thorax.by.girth,res.L5.thorax.by.sex,res.L5.thorax.by.txPC1,res.L5.thorax.by.wingpadPC1,res.ovaries.by.stage,res.testes.by.stage,res.thorax.by.stage,
     file = "model.output/dge.analysis.res.rda")
save(deg100.adult.gonad.by.food,deg100.adult.gonad.by.food_density,deg100.adult.gonad.by.girth,deg100.adult.gonad.by.morph,deg100.adult.gonad.by.sex,deg100.adult.gonad.by.txPC1,deg100.adult.gonad.by.wingPC1,deg100.adult.ovaries.by.food_density,deg100.adult.ovaries.by.girth,deg100.adult.ovaries.by.morph,deg100.adult.ovaries.by.wingPC1,deg100.adult.testes.by.food_density,deg100.adult.testes.by.girth,deg100.adult.testes.by.morph,deg100.adult.testes.by.wingPC1,deg100.adult.thorax.by.food,deg100.adult.thorax.by.food_density,deg100.adult.thorax.by.girth,deg100.adult.thorax.by.morph,deg100.adult.thorax.by.sex,deg100.adult.thorax.by.txPC1,deg100.adult.thorax.by.wingPC1,deg100.L5.gonad.by.food,deg100.L5.gonad.by.food_density,deg100.L5.gonad.by.girth,deg100.L5.gonad.by.sex,deg100.L5.gonad.by.txPC1,deg100.L5.gonad.by.wingpadPC1,deg100.L5.ovaries.by.food,deg100.L5.ovaries.by.food_density,deg100.L5.ovaries.by.girth,deg100.L5.testes.by.food_density,deg100.L5.testes.by.girth,deg100.L5.thorax.by.food,deg100.L5.thorax.by.food_density,deg100.L5.thorax.by.girth,deg100.L5.thorax.by.sex,deg100.L5.thorax.by.txPC1,deg100.L5.thorax.by.wingpadPC1,deg100.ovaries.by.stage,deg100.testes.by.stage,deg100.thorax.by.stage,
     file = "model.output/dge.analysis.deg100s.rda")

# Record all the annotated DEGs into one file
x <- unlist(lapply(ls()[grep("^res\\.",ls())], function (x) {
  cat(x,"\n")
  res <- get(x)
  if (sum(res$padj<0.05, na.rm = TRUE) < 1) {
    s <- "No DEGs!"
    names(s) <- paste0("\n### ",x)
    return(s)
  } 
  x <- sub("^res\\.","",x)
  x <- gsub("\\."," ",x)
  x <- gsub("_"," ",x)
  x <- gsub("txPC1","thoraxPC1",x)
  x <- gsub("PC1"," PC1",x)
  x <- sub("food$","food regime",x)
  res <- res[which(res$padj<0.05),]
  res <- res[order(res$padj),]
  res$baseMean <- signif(res$baseMean, 5)
  res$log2FoldChange <- signif(res$log2FoldChange, 5)
  res$padj <- signif(res$padj, 5)
  res <- apply.annotation(res[,1:5], ann)
  res <- add.manual.annotations(res)
  res$rank <- 1:dim(res)[1]
  COLS <- c(22,1,2,5,21,14:15,20,8,9)
  s <- apply(res[,COLS], 1, paste0, collapse = "\t")
  s <- s[which(!grepl("\\\tNA\\\t\\\t\\\t\\\t\\\tFALSE",s))]
  s <- s[which(!grepl("\\\tNA\\\t\\\t\\\t\\\tNo\\\tFALSE",s))]
  s <- paste0(paste0(x,"\t",names(s),"\t",s), collapse = "\n")
  names(s) <- paste0("\n### ",x)
  return(s)
}) )
write(paste0(c("model","transcript","DEG rank", "baseMean", "log2FoldChange", "padj","manual.ann","EggNOG.Predicted.Gene","EggNOG.Description","EggNOG.Protein.Domains","Contaminant","xenic"),collapse="\t"),
      file = "annotated.deg.list.tsv", sep = "\n", append = FALSE)
write.table(tibble::enframe(x),
            file = "annotated.deg.list.tsv", sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

annotated.deg.list <- read.delim("annotated.deg.list.tsv", blank.lines.skip = TRUE, comment.char = "#")

# Separated by model and filtering for abs(LFC) > 1
{
  x <- annotated.deg.list %>% filter(grepl("by sex",model)) %>% filter(abs(log2FoldChange) > 1) 
  write.table(x, file = "annotated.degs.by.sex.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
  
  x <- annotated.deg.list %>% filter(grepl("by morph",model)) %>% filter(abs(log2FoldChange) > 1) 
  write.table(x, file = "annotated.degs.by.morph.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = FALSE)
  x <- annotated.deg.list %>% filter(grepl("by wing",model)) %>% filter(abs(log2FoldChange) > 1) 
  write.table(x, file = "annotated.degs.by.morph.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  x <- annotated.deg.list %>% filter(grepl("by thorax",model)) %>% filter(abs(log2FoldChange) > 1) 
  write.table(x, file = "annotated.degs.by.morph.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
  x <- annotated.deg.list %>% filter(grepl("by food",model)) %>% filter(abs(log2FoldChange) > 1) 
  write.table(x, file = "annotated.degs.by.food.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = FALSE)
  
  x <- annotated.deg.list %>% filter(grepl("by stage",model)) %>% filter(abs(log2FoldChange) > 1) 
  write.table(x, file = "annotated.degs.by.stage.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = FALSE)
}

# An example query of the annotated DEG list
query.string <- "dsx"
annotated.deg.list %>% 
  filter(grepl("by sex",model)) %>% 
  filter(grepl(query.string,manual.ann, ignore.case = TRUE) | 
         grepl(query.string,EggNOG.Predicted.Gene, ignore.case = TRUE) |
         grepl(query.string,EggNOG.Description, ignore.case = TRUE) |
         grepl(query.string,EggNOG.Protein.Domains, ignore.case = TRUE)) %>% 
  select(-EggNOG.Description, -EggNOG.Protein.Domains, -Contaminant, -xenic) %>% 
  arrange(transcript)


###################
# Volcano Plots
###################

volcano.plot <- function(
  res, pCutoff = 1e-3, FCcutoff = 1, reverseLFC = FALSE,
  NegLFClabel = "LFC < 0", PosLFClabel = "LFC > 0",
  n = NULL,
  title = NULL, caption = NULL,
  text.col = "grey15", max.overlaps = 100) 
{
  require(magrittr)
  require(ggplot2)
  require(ggrepel)
  
  if (reverseLFC) { res$log2FoldChange <- -res$log2FoldChange }
  res <- apply.annotation(as.data.frame(res), ann)
  res$highlight <- (abs(res$log2FoldChange)>FCcutoff & res$padj<pCutoff ) 
  if (!("manual.ann" %in% colnames(res))) {
    cat("Adding annotations.\n")
    res <- add.manual.annotations(res)
  }
  # Over-write those not highlighted in this plot
  plot.labels <- res$manual.ann
  plot.labels[which(!res$highlight)] <- NA
  
  # Test reporting the dataset
  if (is.null(n)) { deg.text <- "" } else { deg.text <- paste0("Samples: ",n,"; ") }
  deg.text <- paste0(deg.text,"Transcripts: ",dim(res)[1],"\n")
  # Text reporting the number of up and down regulated DEGs
  degs.up <- sum(res$padj < pCutoff & res$log2FoldChange > 0, na.rm = TRUE)
  degs.down <- sum(res$padj < pCutoff & res$log2FoldChange < 0, na.rm = TRUE)
  degs.up.pct <- paste0(" (",signif((degs.up/dim(res)[1])*100,2),"%)")
  degs.down.pct <- paste0(" (",signif((degs.down/dim(res)[1])*100,2),"%)")
  deg.text <- paste0(deg.text,
    NegLFClabel,"-biased: ",degs.down,degs.down.pct,
    ", ",PosLFClabel,"-biased: ",degs.up,degs.up.pct,
    ", at p-value cutoff ",pCutoff )
  deg.text <- gsub("0 \\(0%\\)","0",deg.text)

  # Text indicating right and left biases in the LFC
  NegLFClabel <- grobTree(textGrob(NegLFClabel, x=0.05,  y=0.925, hjust=0, gp=gpar(col="grey70", fontsize=12, fontface="bold")))
  PosLFClabel  <- grobTree(textGrob(PosLFClabel,  x=0.95,  y=0.925, hjust=1, gp=gpar(col="grey70", fontsize=12, fontface="bold")))
  
  res %>% 
    ggplot(aes(x=log2FoldChange, y=-log10(padj), color=highlight)) +
    theme_bw() +
    theme(
      plot.title = element_text(size=12),
      legend.position="none",
      panel.grid.minor = element_blank() ) +
    geom_vline(xintercept = 0, color = "grey50") +
    annotation_custom(NegLFClabel) +
    annotation_custom(PosLFClabel) +
    geom_point(size = 0.5, alpha=0.65) +
    geom_text_repel(aes(label=plot.labels), size = 2.5,
                    color = text.col, segment.color = text.col,
                    max.overlaps = max.overlaps,
                    max.time = 1, max.iter = 2e5,
                    box.padding	= 0.125, point.padding = 0,
                    force = 1.25, force_pull = 1 ) +
    scale_color_manual(values=c("grey40","#A60021")) +
    labs(x="log2 fold change", y = "-log10 p",
         title = title, caption = deg.text)
} # End function

# Volcano plots by sex
{
  volcano.plots.by.sex <- list()
  volcano.plots.by.sex[[1]] <- 
    volcano.plot(res.L5.thorax.by.sex[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in juvenile thorax ~ batch + sex",
                 NegLFClabel = "male", PosLFClabel = "female", n = 51)
  volcano.plots.by.sex[[2]] <- 
    volcano.plot(res.adult.thorax.by.sex[,1:5], reverseLFC = TRUE,
                 title = "Expression in adult thorax ~ batch + sex",
                 NegLFClabel = "male", PosLFClabel = "female", n = 91)
  volcano.plots.by.sex[[3]] <- 
    volcano.plot(res.L5.gonad.by.sex[,1:5], max.overlaps = 35, reverseLFC = TRUE,
                 title = "Expression in juvenile gonad ~ batch + sex",
                 NegLFClabel = "male", PosLFClabel = "female", n = 49)
  volcano.plots.by.sex[[4]] <- 
    volcano.plot(res.adult.gonad.by.sex[,1:5], max.overlaps = 20, reverseLFC = TRUE,
                 title = "Expression in adult gonad ~ batch + sex",
                 NegLFClabel = "male", PosLFClabel = "female", n = 89)
  volcano.plots.by.sex.arranged <- ggarrange(plotlist = volcano.plots.by.sex, ncol = 2, nrow = 2, labels="AUTO")
  ggsave("plots/volcano.plots.by.sex.jpg", plot = volcano.plots.by.sex.arranged, width = 6.5, height = 6.5, scale = 1.5)
  ggsave("plots/volcano.plots.by.sex.pdf", plot = volcano.plots.by.sex.arranged, width = 6.5, height = 6.5, scale = 1.5)
}

# Volcano plots by wingPC1
{
  volcano.plots.by.wingPC1 <- list()
  volcano.plots.by.wingPC1[[1]] <- 
    volcano.plot(res.adult.thorax.by.wingPC1[,1:5], reverseLFC = TRUE, max.overlaps = 26,
                 title = "Expression in adult thorax ~ batch + sex + wing shape",
                 NegLFClabel = "short wing", PosLFClabel = "long wing", n=90)
  volcano.plots.by.wingPC1[[2]] <-
    volcano.plot(res.adult.gonad.by.wingPC1[,1:5], reverseLFC = TRUE, pCutoff = 0.05,
                 title = "Expression in adult gonad ~ batch + sex + wing shape",
                 NegLFClabel = "short wing", PosLFClabel = "long wing", n=88)
  # volcano.plots.by.wingPC1[[3]] <- 
  #   volcano.plot(res.adult.ovaries.by.wingPC1[,1:5], reverseLFC = TRUE, pCutoff = 0.05, 
  #                title = "Expression in adult ovaries ~ batch + wing shape",
  #                NegLFClabel = "short wing", PosLFClabel = "long wing", n=44)
  # volcano.plots.by.wingPC1[[4]] <- 
  #   volcano.plot(res.adult.testes.by.wingPC1[,1:5], reverseLFC = TRUE, pCutoff = 0.05, 
  #                title = "Expression in adult testes ~ batch + wing shape",
  #                NegLFClabel = "short wing", PosLFClabel = "long wing", n=44)
  # volcano.plots.by.wingPC1.arranged <- ggarrange(plotlist = volcano.plots.by.wingPC1, ncol = 2, nrow = 2, labels="AUTO")
  # ggsave("plots/volcano.plots.by.wingPC1.jpg", plot = volcano.plots.by.wingPC1.arranged, width = 6.5, height = 6.5, scale = 1.5)
  # ggsave("plots/volcano.plots.by.wingPC1.pdf", plot = volcano.plots.by.wingPC1.arranged, width = 6.5, height = 6.5, scale = 1.5)
}

# Volcano plots by txPC1
{
  volcano.plots.by.txPC1 <- list()
  volcano.plots.by.txPC1[[1]] <- 
    volcano.plot(res.adult.thorax.by.txPC1[,1:5], pCutoff = 0.05, 
                 title = "Expression in adult thorax ~ batch + sex + thorax shape",
                 NegLFClabel = "narrower thorax", PosLFClabel = "wider thorax", n=90)
  volcano.plots.by.txPC1[[2]] <-
    volcano.plot(res.adult.gonad.by.txPC1[,1:5], pCutoff = 0.05,
                 title = "Expression in adult gonad ~ batch + sex + thorax shape",
                 NegLFClabel = "narrower thorax", PosLFClabel = "wider thorax", n=88)
  # volcano.plots.by.txPC1.arranged <- ggarrange(plotlist = volcano.plots.by.txPC1, ncol = 2, nrow = 1, labels="AUTO")
  # ggsave("plots/volcano.plots.by.txPC1.jpg", plot = volcano.plots.by.txPC1.arranged, width = 6.5, height = 3.25, scale = 1.5)
  # ggsave("plots/volcano.plots.by.txPC1.pdf", plot = volcano.plots.by.txPC1.arranged, width = 6.5, height = 3.25, scale = 1.5)
}

# Volcano plots by morph
{
  volcano.plots.by.morph <- list()
  volcano.plots.by.morph[[1]] <- 
    volcano.plot(res.adult.thorax.by.morph[,1:5], max.overlaps = 30, 
                 title = "Expression in adult thorax ~ batch + sex + morph",
                 NegLFClabel = "short wing", PosLFClabel = "long wing", n=91)
  volcano.plots.by.morph[[2]] <-
    volcano.plot(res.adult.gonad.by.morph[,1:5], pCutoff = 0.05,
                 title = "Expression in adult gonad ~ batch + sex + morph",
                 NegLFClabel = "short wing", PosLFClabel = "long wing", n=89)
  # volcano.plots.by.morph[[3]] <- 
  #   volcano.plot(res.adult.ovaries.by.morph[,1:5], pCutoff = 0.05, 
  #                title = "Expression in adult ovaries ~ batch + morph",
  #                NegLFClabel = "short wing", PosLFClabel = "long wing", n=45)
  # volcano.plots.by.morph[[4]] <- 
  #   volcano.plot(res.adult.testes.by.morph[,1:5], pCutoff = 0.05, 
  #                title = "Expression in adult testes ~ batch + morph",
  #                NegLFClabel = "short wing", PosLFClabel = "long wing", n=44)
  # volcano.plots.by.morph.arranged <- ggarrange(plotlist = volcano.plots.by.morph, ncol = 2, nrow = 2, labels="AUTO")
  # ggsave("plots/volcano.plots.by.morph.jpg", plot = volcano.plots.by.morph.arranged, width = 6.5, height = 6.5, scale = 1.5)
  # ggsave("plots/volcano.plots.by.morph.pdf", plot = volcano.plots.by.morph.arranged, width = 6.5, height = 6.5, scale = 1.5)
}

# Plots of adult thorax expression by morph, wing & thorax shapes
{
  volcano.plots.thorax.by.morph.alt <- ggarrange(volcano.plots.by.morph[[1]],volcano.plots.by.wingPC1[[1]],volcano.plots.by.txPC1[[1]],
                                          ncol = 3, nrow = 1, labels="AUTO")
  ggsave("plots/volcano.plots.thorax.by.morph.alt.jpg", plot = volcano.plots.thorax.by.morph.alt, width = 9.75, height = 3.25, scale = 1.5)
  ggsave("plots/volcano.plots.thorax.by.morph.alt.pdf", plot = volcano.plots.thorax.by.morph.alt, width = 9.75, height = 3.25, scale = 1.5)
}

# Plots of adult gonad expression by morph, wing & thorax shapes
{
  volcano.plots.gonad.by.morph.alt <- ggarrange(volcano.plots.by.morph[[2]],volcano.plots.by.wingPC1[[2]],volcano.plots.by.txPC1[[2]],
                                          ncol = 3, nrow = 1, labels="AUTO")
  ggsave("plots/volcano.plots.gonad.by.morph.alt.jpg", plot = volcano.plots.gonad.by.morph.alt, width = 9.75, height = 3.25, scale = 1.5)
  ggsave("plots/volcano.plots.gonad.by.morph.alt.pdf", plot = volcano.plots.gonad.by.morph.alt, width = 9.75, height = 3.25, scale = 1.5)
}

# Plots of juvenile expression by wing pad and thorax shapes -- Not very interesting!
{
  volcano.plots.for.juv.shape <- list()
  volcano.plots.for.juv.shape[[1]] <- 
    volcano.plot(res.L5.thorax.by.wingpadPC1[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in juvenile thorax ~ batch + sex + wing pad shape",
                 NegLFClabel = "narrow", PosLFClabel = "wide", n = 51)
  volcano.plots.for.juv.shape[[2]] <- 
    volcano.plot(res.L5.gonad.by.wingpadPC1[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in juvenile gonad ~ batch + sex + wing pad shape",
                 NegLFClabel = "narrow", PosLFClabel = "wide", n = 60)
  volcano.plots.for.juv.shape[[3]] <- 
    volcano.plot(res.L5.thorax.by.txPC1[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in juvenile thorax ~ batch + sex + thorax shape",
                 NegLFClabel = "narrower mesonotum", PosLFClabel = "wider mesonotum", n = 50)
  volcano.plots.for.juv.shape[[4]] <- 
    volcano.plot(res.L5.gonad.by.txPC1[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in juvenile gonad ~ batch + sex + thorax shape",
                 NegLFClabel = "narrower mesonotum", PosLFClabel = "wider mesonotum", n = 48)
  volcano.plots.for.juv.shape.arranged <- ggarrange(plotlist = volcano.plots.for.juv.shape, ncol = 2, nrow = 2, labels="AUTO")
  ggsave("plots/volcano.plots.for.juv.shape.jpg", plot = volcano.plots.for.juv.shape.arranged, width = 6.5, height = 6.5, scale = 1.5)
  ggsave("plots/volcano.plots.for.juv.shape.pdf", plot = volcano.plots.for.juv.shape.arranged, width = 6.5, height = 6.5, scale = 1.5)
}

# Volcano plots by food density
{
  volcano.plots.by.food_density <- list()
  volcano.plots.by.food_density[[1]] <- 
    volcano.plot(res.L5.thorax.by.food_density[,1:5], pCutoff = 0.05,
                 title = "Expression in juvenile thorax ~ batch + sex + food density",
                 NegLFClabel = "less food", PosLFClabel = "more food", n=51)
  volcano.plots.by.food_density[[2]] <- 
    volcano.plot(res.adult.thorax.by.food_density[,1:5], pCutoff = 0.05,
                 title = "Expression in adult thorax ~ batch + sex + morph + food density",
                 NegLFClabel = "less food", PosLFClabel = "more food", n=91)
  volcano.plots.by.food_density[[3]] <- 
    volcano.plot(res.L5.gonad.by.food_density[,1:5], pCutoff = 0.05,
                 title = "Expression in juvenile gonad ~ batch + sex + food density",
                 NegLFClabel = "less food", PosLFClabel = "more food", n=49)
  volcano.plots.by.food_density[[4]] <- 
    volcano.plot(res.adult.gonad.by.food_density[,1:5], pCutoff = 0.05,
                 title = "Expression in adult gonad ~ batch + sex + morph + food density",
                 NegLFClabel = "less food", PosLFClabel = "more food", n=90)
  volcano.plots.by.food_density.arranged <- ggarrange(plotlist = volcano.plots.by.food_density, ncol = 2, nrow = 2, labels="AUTO")
  ggsave("plots/volcano.plots.by.food_density.jpg", plot = volcano.plots.by.food_density.arranged, width = 6.5, height = 6.5, scale = 1.5)
  ggsave("plots/volcano.plots.by.food_density.pdf", plot = volcano.plots.by.food_density.arranged, width = 6.5, height = 6.5, scale = 1.5)
}

# Volcano plots by food regime
{
  volcano.plots.by.food_regime <- list()
  volcano.plots.by.food_regime[[1]] <- 
    volcano.plot(res.L5.thorax.by.food[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in juvenile thorax ~ batch + sex + food regime",
                 NegLFClabel = "low food", PosLFClabel = "high food", n=51)
  volcano.plots.by.food_regime[[2]] <- 
    volcano.plot(res.adult.thorax.by.food[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in adult thorax ~ batch + sex + morph + food regime",
                 NegLFClabel = "low food", PosLFClabel = "high food", n=91)
  volcano.plots.by.food_regime[[3]] <- 
    volcano.plot(res.L5.gonad.by.food[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in juvenile gonad ~ batch + sex + food regime",
                 NegLFClabel = "low food", PosLFClabel = "high food", n=49)
  volcano.plots.by.food_regime[[4]] <- 
    volcano.plot(res.adult.gonad.by.food[,1:5], pCutoff = 0.05, reverseLFC = TRUE,
                 title = "Expression in adult gonad ~ batch + sex + morph + food regime",
                 NegLFClabel = "low food", PosLFClabel = "high food", n=90)
  volcano.plots.by.food_regime.arranged <- ggarrange(plotlist = volcano.plots.by.food_regime, ncol = 2, nrow = 2, labels="AUTO")
  ggsave("plots/volcano.plots.by.food_regime.jpg", plot = volcano.plots.by.food_regime.arranged, width = 6.5, height = 6.5, scale = 1.5)
  ggsave("plots/volcano.plots.by.food_regime.pdf", plot = volcano.plots.by.food_regime.arranged, width = 6.5, height = 6.5, scale = 1.5)
}

# Volcano plots by stage
{
  volcano.plots.by.stage <- list()
  volcano.plots.by.stage[[1]] <- 
    volcano.plot(res.thorax.by.stage[,1:5], max.overlaps = 40, reverseLFC = TRUE,
                 title = "Expression in thorax ~ batch + sex + stage",
                 NegLFClabel = "juvenile", PosLFClabel = "adult", n=142)
  # volcano.plots.by.stage[[2]] <- ggplot() + theme_void() + geom_blank()
  volcano.plots.by.stage[[2]] <- 
    volcano.plot(res.ovaries.by.stage[,1:5], max.overlaps = 30, reverseLFC = TRUE, 
                 title = "Expression in ovaries ~ batch + stage",
                 NegLFClabel = "juvenile", PosLFClabel = "adult", n=68)
  volcano.plots.by.stage[[3]] <- 
    volcano.plot(res.testes.by.stage[,1:5], max.overlaps = 30, reverseLFC = TRUE,
                 title = "Expression in testes ~ batch + stage",
                 NegLFClabel = "juvenile", PosLFClabel = "adult", n=71)
  volcano.plots.by.stage.arranged <- ggarrange(plotlist = volcano.plots.by.stage, ncol = 3, nrow = 1, labels="AUTO")
  ggsave("plots/volcano.plots.by.stage.jpg", plot = volcano.plots.by.stage.arranged, width = 9.75, height = 3.25, scale = 1.5)
  ggsave("plots/volcano.plots.by.stage.pdf", plot = volcano.plots.by.stage.arranged, width = 9.75, height = 3.25, scale = 1.5)
}

###################
# Enrichment Analysis
###################
# Enrichment Analysis for GO, KEGG and EggNOG Protein Domain terms 
# was performed using the hypergeometric test as implemented in stats::phyper
# 
# source('enrichment.R')

load("enrichment.pvalues.rda",verbose = TRUE)

###################
# The End!
###################

sessionInfo()
