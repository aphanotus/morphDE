rm(list=ls())
gc()

library(tidyverse)
library(viridis)
library(DESeq2)
library(IHW)
library(pheatmap)

###################
# Import the Kallisto gene count matrix
###################
cts <- read.csv("kallisto.counts.csv")
dim(cts)
# 395125    370
head(names(cts))
tail(names(cts))

# Exclude samples of uncertain tissue
x <- c(grep("FC355_06",colnames(cts)),grep("Undetermined",colnames(cts)))
cts <- cts[,-x]

sum(is.na(cts))
# 0

# Read sums by sample
hist(colSums(cts))
barplot(colSums(cts))
abline(h=summary(colSums(cts))[c(2,5)], col="darkred", lty=3)

# What's the range of counts for transcripts?
summary(rowSums(cts))

# Read sums by transcript
hist(log10(rowSums(cts)+1))

# Coefficient of variation (CV) in transcript counts by sample,
# for those with at least 100 counts
x <- cts[which(rowSums(cts)>=100),]
cv <- function (x) { sd(x) / mean (x) }
hist(apply(x, 1, cv))

# Sample metadata
sample.metadata <- read.csv("3seq.metadata.csv", stringsAsFactors = FALSE)
# Metadata rows much match order of samples (columns) in the count matrix
all(colnames(cts) %in% sample.metadata$kallisto_name)

sample.metadata <- sample.metadata %>%
  column_to_rownames(var = "kallisto_name") %>%
  select(class, population, stage, sex, morph, morph_sex, tissue, food_regime,
         cohort, seeds, plate, well, seq_date) %>%
  mutate(plate = as.factor(plate), seq_date = as.factor(seq_date)) %>%
  mutate(tissue = fct_relevel(tissue, "thorax")) %>%
  mutate(morph = fct_relevel(morph, "SW", "LW"))

# Apply filtering criteria
# remove redundant samples (not sure right now if the June '20 or July '20 reads are better)
x <- sample.metadata %>%
  # filter(seq_date==200611) %>%
  filter(seq_date==200701) %>%
  rownames()
x <- which(colnames(cts) %in% x)
cts <- cts[,-x]
sample.metadata <- sample.metadata %>%
  filter(seq_date!=200701)

# Note that DESeq2 recommends only minimal filtering of low-count genes
gene.read.number.cut.off <- 100

###################
# Gene annotations
###################
chr.scafs <- c(
  "Scaffold_8819;HRSCAF=9471" = "Chr4",
  "Scaffold_16598;HRSCAF=18332" = "Chr5",
  "Scaffold_23484;HRSCAF=26600" = "ChrM",
  "Scaffold_49954;HRSCAF=55183" = "Chr1",
  "Scaffold_50151;HRSCAF=55390" = "ChrX",
  "Scaffold_50591;HRSCAF=55866" = "Chr2",
  "Scaffold_50594;HRSCAF=55992" = "Chr3"
)

ann <- read.delim("final_annotations_lvl0_postive.tsv")
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

mapping <- read.delim("JhaeGU.c90.cds.vs.gmap_JhaeHiC.positions.tsv", header = FALSE)
colnames(mapping) <- c("qname","flag","rname","pos")
mapping$qname <- str_split_fixed(mapping$qname, "::",3)[,2]
mapping$rname <- str_replace_all(mapping$rname, chr.scafs)

apply.annotation <- function(df, ann, mapping) {
  i <- match(rownames(df), mapping$qname)
  df$chr <- mapping$rname[i]
  df$pos <- mapping$pos[i]
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
# # Differential expression analysis
###################

# See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

###################
### Null model
###################

#### Sub-setting the data and metadata
cts.i <- as.matrix(cts[which(rowSums(cts) > gene.read.number.cut.off),])
dim(cts.i)

# Why a null model accounting for batch is critical!
uncorrected.pca <- prcomp(t(cts.i))
eigs <- uncorrected.pca$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)

uncorrected.pca.by.biol.grp.plot <- uncorrected.pca$x %>%
  as.data.frame() %>%
  ggplot(aes(x=PC1, y=PC2)) +
  theme_bw() + theme(legend.position="bottom") +
  geom_point(aes(color=as.factor(grp)), size = 3, alpha = 0.85) +
  scale_color_viridis(name = NULL, discrete = TRUE, begin = 0, end = 1) +
  labs(x="PC1 (78.4% variance)", y="PC2 (20.1% variance)") +
  coord_fixed()
uncorrected.pca.by.biol.grp.plot

uncorrected.pca.by.batch.plot <- uncorrected.pca$x %>%
  as.data.frame() %>%
  ggplot(aes(x=PC1, y=PC2)) +
  theme_bw() + theme(legend.position="bottom") +
  geom_point(aes(color=as.factor(as.numeric(sample.metadata$plate))), size = 3, alpha = 0.85) +
  scale_color_viridis(name = "batch", discrete = TRUE, begin = 0.2, end = 0.8, option = "magma") +
  labs(x="PC1 (78.4% variance)", y="PC2 (20.1% variance)") +
  coord_fixed()
uncorrected.pca.by.batch.plot

uncorrected.pca.plots <- ggpubr::ggarrange(
  uncorrected.pca.by.biol.grp.plot, uncorrected.pca.by.batch.plot,
  ncol = 2)

ggsave("images/uncorrected.pca.plots.pdf", uncorrected.pca.plots, width = 14, height = 5, scale = 1)
ggsave("images/uncorrected.pca.plots.jpg", uncorrected.pca.plots, width = 14, height = 5, scale = 1)

# DEseq2
dds.null <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=sample.metadata,
  design = ~ plate )
dds.null <- DESeq(dds.null)
dds.null <- estimateSizeFactors(dds.null)

# Extract and save normalized counts
ctn.null <- counts(dds.null, normalized=TRUE)
write.csv(ctn.null, 'normalized.counts.null.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.null)
res.null <- results(dds.null, contrast=c("plate","3","2"), alpha=0.05, filterFun=ihw)
res.null <- lfcShrink(dds.null, coef="plate_3_vs_2", type="apeglm", res = res.null)
summary(res.null)

# We can order our results table by the smallest p value:
deg100.null23 <- res.null[order(res.null$padj),]
deg100.null23 <- deg100.null23[1:100,]
deg100.null23 <- apply.annotation(deg100.null23, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.null23), "deg100.null23.csv")

res.null <- results(dds.null, contrast=c("plate","4","2"), alpha=0.05, filterFun=ihw)
res.null <- lfcShrink(dds.null, coef="plate_4_vs_2", type="apeglm", res = res.null)
summary(res.null)

# We can order our results table by the smallest p value:
deg100.null24 <- res.null[order(res.null$padj),]
deg100.null24 <- deg100.null24[1:100,]
deg100.null24 <- apply.annotation(deg100.null24, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.null23), "deg100.null23.csv")

vsd.null <- vst(dds.null, blind=FALSE)

### heirarchical clustering
vsd.cor.null <- cor(assay(vsd.null))
colnames(sample.metadata)
sample.metadata %>% select(7,3,4,5,2,8) %>%
  pheatmap(vsd.cor.null, annotation = ., labels_row = " ", labels_col = " ")

# PCA
pca <- prcomp(t(assay(vsd.null)))
eigs <- pca$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)

grp <- with(sample.metadata, paste(stage,sex,tissue))
grp <- sub("abdomen","gonad",grp)
grp <- sub("L5"," juvenile",grp)
grp <- sub(" f "," female ",grp)
grp <- sub(" m "," male ",grp)
unique(grp)

null.pca.plot <- pca$x %>%
  as.data.frame() %>%
  ggplot(aes(x=PC1, y=PC2)) +
  theme_bw() + theme(legend.position="bottom") +
  geom_point(aes(color=as.factor(grp)), size = 3, alpha = 0.85) +
  scale_color_viridis(name = NULL, discrete = TRUE, begin = 0, end = 1) +
  labs(x="PC1 (33.5% variance", y="PC2 (7.5% variance)") +
  coord_fixed()
null.pca.plot
ggsave("images/null.pca.plot.pdf", null.pca.plot, width = 7, height = 5, scale = 1)
ggsave("images/null.pca.plot.jpg", null.pca.plot, width = 7, height = 5, scale = 0.95)

rm(cts.i, dds.null, res.null, deg100.null23, deg100.null24)



###################
### M v. F for adult gonad
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.gonad.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.adult.gonad.by.sex <- DESeq(dds.adult.gonad.by.sex)
dds.adult.gonad.by.sex <- estimateSizeFactors(dds.adult.gonad.by.sex)

# Extract and save normalized counts
ctn.adult.gonad.by.sex <- counts(dds.adult.gonad.by.sex, normalized=TRUE)
write.csv(ctn.adult.gonad.by.sex, 'normalized.counts.adult.gonad.by.sex.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.sex)
res.adult.gonad.by.sex <- results(dds.adult.gonad.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.sex <- lfcShrink(dds.adult.gonad.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.adult.gonad.by.sex)
summary(res.adult.gonad.by.sex)
plotMA(res.adult.gonad.by.sex, main="apeGLM")
# interactively detect the row number
i <- identify(res.adult.gonad.by.sex$baseMean, res.adult.gonad.by.sex$log2FoldChange)
plotCounts(dds.adult.gonad.by.sex, gene=i[1], intgroup="sex", main="adult.gonad.by.sex")
plotCounts(dds.adult.gonad.by.sex, gene=which.min(res.adult.gonad.by.sex$padj), intgroup="sex", main="adult.gonad.by.sex")
# Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.sex <- vst(dds.adult.gonad.by.sex, blind=FALSE)
plotPCA(vsd.adult.gonad.by.sex, intgroup=c("sex"))
# Check the sex on the contrary samples!

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.sex <- res.adult.gonad.by.sex[order(res.adult.gonad.by.sex$padj),]
deg100.adult.gonad.by.sex <- deg100.adult.gonad.by.sex[1:100,]
deg100.adult.gonad.by.sex <- apply.annotation(deg100.adult.gonad.by.sex, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.adult.gonad.by.sex), "deg100.adult.gonad.by.sex.csv")
save(meta.i, cts.i, dds.adult.gonad.by.sex, ctn.adult.gonad.by.sex, res.adult.gonad.by.sex, vsd.adult.gonad.by.sex, deg100.adult.gonad.by.sex,
     file = "adult.gonad.by.sex.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.sex, ctn.adult.gonad.by.sex, res.adult.gonad.by.sex, vsd.adult.gonad.by.sex)
load("adult.gonad.by.sex.rda", verbose = TRUE)

###################
### M v. F for L5 gonad
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "L5", tissue != "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.L5.gonad.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.L5.gonad.by.sex <- DESeq(dds.L5.gonad.by.sex)
dds.L5.gonad.by.sex <- estimateSizeFactors(dds.L5.gonad.by.sex)

# Extract and save normalized counts
ctn.L5.gonad.by.sex <- counts(dds.L5.gonad.by.sex, normalized=TRUE)
write.csv(ctn.L5.gonad.by.sex, 'normalized.counts.L5.gonad.by.sex.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.sex)
res.L5.gonad.by.sex <- results(dds.L5.gonad.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.sex <- lfcShrink(dds.L5.gonad.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.L5.gonad.by.sex)
summary(res.L5.gonad.by.sex)
plotMA(res.L5.gonad.by.sex, main="apeGLM")
# Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.sex <- vst(dds.L5.gonad.by.sex, blind=FALSE)
plotPCA(vsd.L5.gonad.by.sex, intgroup=c("sex"))
# Check the sex on the contrary samples!

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.sex <- res.L5.gonad.by.sex[order(res.L5.gonad.by.sex$padj),]
deg100.L5.gonad.by.sex <- deg100.L5.gonad.by.sex[1:100,]
deg100.L5.gonad.by.sex <- apply.annotation(deg100.L5.gonad.by.sex, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.L5.gonad.by.sex), "deg100.L5.gonad.by.sex.csv")
save(meta.i, cts.i, dds.L5.gonad.by.sex, ctn.L5.gonad.by.sex, res.L5.gonad.by.sex, vsd.L5.gonad.by.sex, deg100.L5.gonad.by.sex,
     file = "L5.gonad.by.sex.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.sex, ctn.L5.gonad.by.sex, res.L5.gonad.by.sex, vsd.L5.gonad.by.sex)
load("L5.gonad.by.sex.rda", verbose = TRUE)

###################
### M v. F for adult thorax
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.thorax.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.adult.thorax.by.sex <- DESeq(dds.adult.thorax.by.sex)
dds.adult.thorax.by.sex <- estimateSizeFactors(dds.adult.thorax.by.sex)

# Extract and save normalized counts
ctn.adult.thorax.by.sex <- counts(dds.adult.thorax.by.sex, normalized=TRUE)
write.csv(ctn.adult.thorax.by.sex, 'normalized.counts.adult.thorax.by.sex.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.sex)
res.adult.thorax.by.sex <- results(dds.adult.thorax.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.sex <- lfcShrink(dds.adult.thorax.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.adult.thorax.by.sex)
summary(res.adult.thorax.by.sex)
plotMA(res.adult.thorax.by.sex, main="apeGLM")
# Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.sex <- vst(dds.adult.thorax.by.sex, blind=FALSE)
plotPCA(vsd.adult.thorax.by.sex, intgroup=c("sex"))
# Check the sex on the contrary samples!

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.sex <- res.adult.thorax.by.sex[order(res.adult.thorax.by.sex$padj),]
deg100.adult.thorax.by.sex <- deg100.adult.thorax.by.sex[1:100,]
deg100.adult.thorax.by.sex <- apply.annotation(deg100.adult.thorax.by.sex, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.adult.thorax.by.sex), "deg100.adult.thorax.by.sex.csv")
save(meta.i, cts.i, dds.adult.thorax.by.sex, ctn.adult.thorax.by.sex, res.adult.thorax.by.sex, vsd.adult.thorax.by.sex, deg100.adult.thorax.by.sex,
     file = "adult.thorax.by.sex.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.sex, ctn.adult.thorax.by.sex, res.adult.thorax.by.sex, vsd.adult.thorax.by.sex)
load("adult.thorax.by.sex.rda", verbose = TRUE)

###################
### M v. F for L5 thorax
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "L5", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.L5.thorax.by.sex <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex )
dds.L5.thorax.by.sex <- DESeq(dds.L5.thorax.by.sex)
dds.L5.thorax.by.sex <- estimateSizeFactors(dds.L5.thorax.by.sex)

# Extract and save normalized counts
ctn.L5.thorax.by.sex <- counts(dds.L5.thorax.by.sex, normalized=TRUE)
write.csv(ctn.L5.thorax.by.sex, 'normalized.counts.L5.thorax.by.sex.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.sex)
res.L5.thorax.by.sex <- results(dds.L5.thorax.by.sex, contrast=c("sex","m","f"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.sex <- lfcShrink(dds.L5.thorax.by.sex, coef="sex_m_vs_f", type="apeglm", res = res.L5.thorax.by.sex)
summary(res.L5.thorax.by.sex)
plotMA(res.L5.thorax.by.sex, main="apeGLM")
# Positive LFC values are associated with male expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.sex <- vst(dds.L5.thorax.by.sex, blind=FALSE)
plotPCA(vsd.L5.thorax.by.sex, intgroup=c("sex"))
# Check the sex on the contrary samples!

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.sex <- res.L5.thorax.by.sex[order(res.L5.thorax.by.sex$padj),]
deg100.L5.thorax.by.sex <- deg100.L5.thorax.by.sex[1:100,]
deg100.L5.thorax.by.sex <- apply.annotation(deg100.L5.thorax.by.sex, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.L5.thorax.by.sex), "deg100.L5.thorax.by.sex.csv")
save(meta.i, cts.i, dds.L5.thorax.by.sex, ctn.L5.thorax.by.sex, res.L5.thorax.by.sex, vsd.L5.thorax.by.sex, deg100.L5.thorax.by.sex,
     file = "L5.thorax.by.sex.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.sex, ctn.L5.thorax.by.sex, res.L5.thorax.by.sex, vsd.L5.thorax.by.sex)
load("L5.thorax.by.sex.rda", verbose = TRUE)

###################
### Adult thorax by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.thorax.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph )
dds.adult.thorax.by.morph <- DESeq(dds.adult.thorax.by.morph)
dds.adult.thorax.by.morph <- estimateSizeFactors(dds.adult.thorax.by.morph)

# Extract and save normalized counts
ctn.adult.thorax.by.morph <- counts(dds.adult.thorax.by.morph, normalized=TRUE)
write.csv(ctn.adult.thorax.by.morph, 'normalized.counts.adult.thorax.by.morph.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.morph)
res.adult.thorax.by.morph <- results(dds.adult.thorax.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.morph <- lfcShrink(dds.adult.thorax.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.adult.thorax.by.morph)
summary(res.adult.thorax.by.morph)
plotMA(res.adult.thorax.by.morph, main="apeGLM")
# interactively detect the row number
i <- identify(res.adult.thorax.by.morph$baseMean, res.adult.thorax.by.morph$log2FoldChange)
res.adult.thorax.by.morph[i[1],]
rownames(ctn.adult.thorax.by.morph)[i[1]]
data.frame(
  counts = ctn.adult.thorax.by.morph[i[2],],
  morph = meta.i$morph) %>%
  ggplot(aes(morph,counts)) +
  geom_violin(fill="gray85", color = NA, alpha = 0.65) +
  geom_point(position=position_jitter(w=0.1,h=0), alpha = 0.65) +
  scale_y_log10()
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.morph <- vst(dds.adult.thorax.by.morph, blind=FALSE)
plotPCA(vsd.adult.thorax.by.morph, intgroup=c("morph"))
# The separation of groups in the NW/SE axis is by plate!

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.morph <- res.adult.thorax.by.morph[order(res.adult.thorax.by.morph$padj),]
deg100.adult.thorax.by.morph <- deg100.adult.thorax.by.morph[1:100,]
deg100.adult.thorax.by.morph <- apply.annotation(deg100.adult.thorax.by.morph, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.adult.thorax.by.morph), "deg100.adult.thorax.by.morph.csv")
save(meta.i, cts.i, dds.adult.thorax.by.morph, ctn.adult.thorax.by.morph, res.adult.thorax.by.morph, vsd.adult.thorax.by.morph, deg100.adult.thorax.by.morph,
     file = "adult.thorax.by.morph.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.morph, ctn.adult.thorax.by.morph, res.adult.thorax.by.morph, vsd.adult.thorax.by.morph)
load("adult.thorax.by.morph.rda", verbose = TRUE)

###################
### Adult gonad by morph
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.gonad.by.morph <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph )
dds.adult.gonad.by.morph <- DESeq(dds.adult.gonad.by.morph)
dds.adult.gonad.by.morph <- estimateSizeFactors(dds.adult.gonad.by.morph)

# Extract and save normalized counts
ctn.adult.gonad.by.morph <- counts(dds.adult.gonad.by.morph, normalized=TRUE)
write.csv(ctn.adult.gonad.by.morph, 'normalized.counts.adult.gonad.by.morph.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.morph)
res.adult.gonad.by.morph <- results(dds.adult.gonad.by.morph, contrast=c("morph","LW","SW"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.morph <- lfcShrink(dds.adult.gonad.by.morph, coef="morph_LW_vs_SW", type="apeglm", res = res.adult.gonad.by.morph)
summary(res.adult.gonad.by.morph)
plotMA(res.adult.gonad.by.morph, main="apeGLM")
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.morph <- vst(dds.adult.gonad.by.morph, blind=FALSE)
plotPCA(vsd.adult.gonad.by.morph, intgroup=c("morph"))
# The separation of groups in the NW/SE axis is by plate!

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.morph <- res.adult.gonad.by.morph[order(res.adult.gonad.by.morph$padj),]
deg100.adult.gonad.by.morph <- deg100.adult.gonad.by.morph[1:100,]
deg100.adult.gonad.by.morph <- apply.annotation(deg100.adult.gonad.by.morph, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.adult.gonad.by.morph), "deg100.adult.gonad.by.morph.csv")
save(meta.i, cts.i, dds.adult.gonad.by.morph, ctn.adult.gonad.by.morph, res.adult.gonad.by.morph, vsd.adult.gonad.by.morph, deg100.adult.gonad.by.morph,
     file = "adult.gonad.by.morph.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.morph, ctn.adult.gonad.by.morph, res.adult.gonad.by.morph, vsd.adult.gonad.by.morph)
load("adult.gonad.by.morph.rda", verbose = TRUE)

###################
### Adult gonad by food regime
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "adult", tissue == "gonad")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.adult.gonad.by.food <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph + food_regime)
dds.adult.gonad.by.food <- DESeq(dds.adult.gonad.by.food)
dds.adult.gonad.by.food <- estimateSizeFactors(dds.adult.gonad.by.food)

# Extract and save normalized counts
ctn.adult.gonad.by.food <- counts(dds.adult.gonad.by.food, normalized=TRUE)
write.csv(ctn.adult.gonad.by.food, 'normalized.counts.adult.gonad.by.food.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.food)
res.adult.gonad.by.food <- results(dds.adult.gonad.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.food <- lfcShrink(dds.adult.gonad.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.adult.gonad.by.food)
summary(res.adult.gonad.by.food)
plotMA(res.adult.gonad.by.food, main="apeGLM")
# Positive LFC values are associated with low food regime expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.food <- vst(dds.adult.gonad.by.food, blind=FALSE)
plotPCA(vsd.adult.gonad.by.food, intgroup=c("food_regime"))
# The separation of groups in the NW/SE axis is by plate!

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.food <- res.adult.gonad.by.food[order(res.adult.gonad.by.food$padj),]
deg100.adult.gonad.by.food <- deg100.adult.gonad.by.food[1:100,]
deg100.adult.gonad.by.food <- apply.annotation(deg100.adult.gonad.by.food, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.adult.gonad.by.food), "deg100.adult.gonad.by.food.csv")
save(meta.i, cts.i, dds.adult.gonad.by.food, ctn.adult.gonad.by.food, res.adult.gonad.by.food, vsd.adult.gonad.by.food, deg100.adult.gonad.by.food,
     file = "adult.gonad.by.food.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.food, ctn.adult.gonad.by.food, res.adult.gonad.by.food, vsd.adult.gonad.by.food)
load("adult.gonad.by.food.rda", verbose = TRUE)

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

# Extract and save normalized counts
ctn.adult.thorax.by.food <- counts(dds.adult.thorax.by.food, normalized=TRUE)
write.csv(ctn.adult.thorax.by.food, 'normalized.counts.adult.thorax.by.food.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.food)
res.adult.thorax.by.food <- results(dds.adult.thorax.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.food <- lfcShrink(dds.adult.thorax.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.adult.thorax.by.food)
summary(res.adult.thorax.by.food)
plotMA(res.adult.thorax.by.food, main="apeGLM")
# Positive LFC values are associated with low food regime expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.food <- vst(dds.adult.thorax.by.food, blind=FALSE)
plotPCA(vsd.adult.thorax.by.food, intgroup=c("food_regime"))
# The separation of groups in the NW/SE axis is by plate!

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.food <- res.adult.thorax.by.food[order(res.adult.thorax.by.food$padj),]
deg100.adult.thorax.by.food <- deg100.adult.thorax.by.food[1:100,]
deg100.adult.thorax.by.food <- apply.annotation(deg100.adult.thorax.by.food, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.adult.thorax.by.food), "deg100.adult.thorax.by.food.csv")
save(meta.i, cts.i, dds.adult.thorax.by.food, ctn.adult.thorax.by.food, res.adult.thorax.by.food, vsd.adult.thorax.by.food, deg100.adult.thorax.by.food,
     file = "adult.thorax.by.food.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.food, ctn.adult.thorax.by.food, res.adult.thorax.by.food, vsd.adult.thorax.by.food)
load("adult.thorax.by.food.rda", verbose = TRUE)

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

# Extract and save normalized counts
ctn.L5.gonad.by.food <- counts(dds.L5.gonad.by.food, normalized=TRUE)
write.csv(ctn.L5.gonad.by.food, 'normalized.counts.L5.gonad.by.food.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.food)
res.L5.gonad.by.food <- results(dds.L5.gonad.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.food <- lfcShrink(dds.L5.gonad.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.L5.gonad.by.food)
summary(res.L5.gonad.by.food)
plotMA(res.L5.gonad.by.food, main="apeGLM")
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.food <- vst(dds.L5.gonad.by.food, blind=FALSE)
plotPCA(vsd.L5.gonad.by.food, intgroup=c("food_regime"))

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.food <- res.L5.gonad.by.food[order(res.L5.gonad.by.food$padj),]
deg100.L5.gonad.by.food <- deg100.L5.gonad.by.food[1:100,]
deg100.L5.gonad.by.food <- apply.annotation(deg100.L5.gonad.by.food, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.L5.gonad.by.food), "deg100.L5.gonad.by.food.csv")
save(meta.i, cts.i, dds.L5.gonad.by.food, ctn.L5.gonad.by.food, res.L5.gonad.by.food, vsd.L5.gonad.by.food, deg100.L5.gonad.by.food,
     file = "L5.gonad.by.food.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.food, ctn.L5.gonad.by.food, res.L5.gonad.by.food, vsd.L5.gonad.by.food)
load("L5.gonad.by.food.rda", verbose = TRUE)

###################
### L5 thorax by food regime
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage == "L5", tissue == "thorax")
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i)
dim(cts.i)
c(with(meta.i, by(food_regime,food_regime,length)))

dds.L5.thorax.by.food <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + food_regime )
dds.L5.thorax.by.food <- DESeq(dds.L5.thorax.by.food)
dds.L5.thorax.by.food <- estimateSizeFactors(dds.L5.thorax.by.food)

# Extract and save normalized counts
ctn.L5.thorax.by.food <- counts(dds.L5.thorax.by.food, normalized=TRUE)
write.csv(ctn.L5.thorax.by.food, 'normalized.counts.L5.thorax.by.food.csv')
# Use normalized counts for downstream visualization

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.food)
res.L5.thorax.by.food <- results(dds.L5.thorax.by.food, contrast=c("food_regime","low","high"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.food <- lfcShrink(dds.L5.thorax.by.food, coef="food_regime_low_vs_high", type="apeglm", res = res.L5.thorax.by.food)
summary(res.L5.thorax.by.food)
plotMA(res.L5.thorax.by.food, main="apeGLM")
# Positive LFC values are associated with LW expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.food <- vst(dds.L5.thorax.by.food, blind=FALSE)
plotPCA(vsd.L5.thorax.by.food, intgroup=c("food_regime"))

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.food <- res.L5.thorax.by.food[order(res.L5.thorax.by.food$padj),]
deg100.L5.thorax.by.food <- deg100.L5.thorax.by.food[1:100,]
deg100.L5.thorax.by.food <- apply.annotation(deg100.L5.thorax.by.food, ann, mapping)
# Save results
write.csv(as.data.frame(deg100.L5.thorax.by.food), "deg100.L5.thorax.by.food.csv")
save(meta.i, cts.i, dds.L5.thorax.by.food, ctn.L5.thorax.by.food, res.L5.thorax.by.food, vsd.L5.thorax.by.food, deg100.L5.thorax.by.food,
     file = "L5.thorax.by.food.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.food, ctn.L5.thorax.by.food, res.L5.thorax.by.food, vsd.L5.thorax.by.food)
load("L5.thorax.by.food.rda", verbose = TRUE)



save(ann, apply.annotation, chr.scafs, cts, deg100.adult.gonad.by.food, deg100.adult.gonad.by.morph, deg100.adult.gonad.by.sex, deg100.adult.thorax.by.food, deg100.adult.thorax.by.morph, deg100.adult.thorax.by.sex, deg100.L5.gonad.by.sex, deg100.L5.thorax.by.sex, deg100.L5.gonad.by.food, deg100.L5.thorax.by.food, gene.read.number.cut.off, mapping, sample.metadata,
     file = "dge.analysis.rda")
save(deg100.adult.gonad.by.food, deg100.adult.gonad.by.morph, deg100.adult.gonad.by.sex, deg100.adult.thorax.by.food, deg100.adult.thorax.by.morph, deg100.adult.thorax.by.sex, deg100.L5.gonad.by.sex, deg100.L5.thorax.by.sex, deg100.L5.gonad.by.food, deg100.L5.thorax.by.food,
     file = "deg100s.rda")
