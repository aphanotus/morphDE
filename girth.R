# Does girth correlate with food density?
sex.stage <- paste(sample.metadata$sex,sample.metadata$stage)
sex.stage <- sub("L5","juvenile",sex.stage)
sex.stage <- sub("^m","male",sex.stage)
sex.stage <- sub("^f","female",sex.stage)
plot.abdominal.width <- 
  ggplot(data = sample.metadata, aes(x = log10(seeds/cohort), y = girth, color = sex.stage)) + 
  theme_bw() +
  geom_point(alpha=0.65) +
  scale_color_viridis(name="sex & stage", discrete = TRUE, end = 0.925) +
  geom_smooth(method=lm, fill="gray90") +
  labs(x="food density (log10 seeds / bug)", y = "abdominal width (mm)")
plot.abdominal.width

ggsave("plots/abdominal.width.pdf",plot.abdominal.width, width = 6, height = 4.5, scale = 1)

m <- lm(girth~stage*sex, data=sample.metadata)
summary(lm(m))

sample.metadata %>% 
  mutate(sex.stage=sex.stage) %>% 
  filter(!is.na(sample.metadata$girth)) %>% 
  mutate(
    food_density = log10(seeds/cohort),
    resids = resid(m)
  ) %>% 
ggplot(aes(x = food_density, y = resids, color = sex.stage)) + 
  theme_bw() +
  geom_point(alpha=0.85) +
  scale_color_viridis(name="sex & stage", discrete = TRUE) +
  geom_smooth(method=lm, fill="gray90") +
  labs(x="food density (log10 seeds / bug)", y = "residuals(abdominal width ~ stage * sex)")

###################
### L5 thorax by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue == "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(sex,sex,length)))

dds.L5.thorax.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + girth )
dds.L5.thorax.by.girth <- DESeq(dds.L5.thorax.by.girth)
dds.L5.thorax.by.girth <- estimateSizeFactors(dds.L5.thorax.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.thorax.by.girth)
res.L5.thorax.by.girth <- results(dds.L5.thorax.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.L5.thorax.by.girth <- lfcShrink(dds.L5.thorax.by.girth, coef="girth", type="apeglm", res = res.L5.thorax.by.girth)
summary(res.L5.thorax.by.girth)
# out of 24958 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 11, 0.044%
# LFC < 0 (down)     : 29, 0.12%
# outliers [1]       : 0, 0%

plotMA(res.L5.thorax.by.girth, main="apeGLM", ylim=c(-3,4))
plot(plotCounts(dds.L5.thorax.by.girth, gene=order(res.L5.thorax.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.L5.thorax.by.girth$log2FoldChange[order(res.L5.thorax.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.thorax.by.girth <- vst(dds.L5.thorax.by.girth, blind=FALSE)
plotPCA(vsd.L5.thorax.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.L5.thorax.by.girth <- res.L5.thorax.by.girth[order(res.L5.thorax.by.girth$padj),]
deg100.L5.thorax.by.girth <- deg100.L5.thorax.by.girth[1:100,]
deg100.L5.thorax.by.girth <- apply.annotation(deg100.L5.thorax.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.L5.thorax.by.girth), "deg100.L5.thorax.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.thorax.by.girth, res.L5.thorax.by.girth, vsd.L5.thorax.by.girth, deg100.L5.thorax.by.girth,
     file = "L5.thorax.by.girth.rda")
rm(meta.i, cts.i, dds.L5.thorax.by.girth, res.L5.thorax.by.girth, vsd.L5.thorax.by.girth)
# load("L5.thorax.by.girth.rda", verbose = TRUE)


###################
### L5 gonad by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", tissue != "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

dds.L5.gonad.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + girth )
dds.L5.gonad.by.girth <- DESeq(dds.L5.gonad.by.girth)
dds.L5.gonad.by.girth <- estimateSizeFactors(dds.L5.gonad.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.gonad.by.girth)
res.L5.gonad.by.girth <- results(dds.L5.gonad.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.L5.gonad.by.girth <- lfcShrink(dds.L5.gonad.by.girth, coef="girth", type="apeglm", res = res.L5.gonad.by.girth)
summary(res.L5.gonad.by.girth)
# out of 32603 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 196, 0.6%
# LFC < 0 (down)     : 106, 0.33%
# outliers [1]       : 0, 0%

plotMA(res.L5.gonad.by.girth, main="apeGLM", ylim=c(-4,6))
plot(plotCounts(dds.L5.gonad.by.girth, gene=order(res.L5.gonad.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.L5.gonad.by.girth$log2FoldChange[order(res.L5.gonad.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.gonad.by.girth <- vst(dds.L5.gonad.by.girth, blind=FALSE)
plotPCA(vsd.L5.gonad.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.L5.gonad.by.girth <- res.L5.gonad.by.girth[order(res.L5.gonad.by.girth$padj),]
deg100.L5.gonad.by.girth <- deg100.L5.gonad.by.girth[1:100,]
deg100.L5.gonad.by.girth <- apply.annotation(deg100.L5.gonad.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.L5.gonad.by.girth), "deg100.L5.gonad.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.gonad.by.girth, res.L5.gonad.by.girth, vsd.L5.gonad.by.girth, deg100.L5.gonad.by.girth,
     file = "L5.gonad.by.girth.rda")
rm(meta.i, cts.i, dds.L5.gonad.by.girth, res.L5.gonad.by.girth, vsd.L5.gonad.by.girth)
# load("L5.gonad.by.girth.rda", verbose = TRUE)


###################
### L5 ovaries by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", sex == "f", tissue != "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

dds.L5.ovaries.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + girth )
dds.L5.ovaries.by.girth <- DESeq(dds.L5.ovaries.by.girth)
dds.L5.ovaries.by.girth <- estimateSizeFactors(dds.L5.ovaries.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.ovaries.by.girth)
res.L5.ovaries.by.girth <- results(dds.L5.ovaries.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.L5.ovaries.by.girth <- lfcShrink(dds.L5.ovaries.by.girth, coef="girth", type="apeglm", res = res.L5.ovaries.by.girth)
summary(res.L5.ovaries.by.girth)
# out of 17635 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 170, 0.96%
# LFC < 0 (down)     : 74, 0.42%
# outliers [1]       : 0, 0%

plotMA(res.L5.ovaries.by.girth, main="apeGLM", ylim=c(-3,5))
plot(plotCounts(dds.L5.ovaries.by.girth, gene=order(res.L5.ovaries.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.L5.ovaries.by.girth$log2FoldChange[order(res.L5.ovaries.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.ovaries.by.girth <- vst(dds.L5.ovaries.by.girth, blind=FALSE)
plotPCA(vsd.L5.ovaries.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.L5.ovaries.by.girth <- res.L5.ovaries.by.girth[order(res.L5.ovaries.by.girth$padj),]
deg100.L5.ovaries.by.girth <- deg100.L5.ovaries.by.girth[1:100,]
deg100.L5.ovaries.by.girth <- apply.annotation(deg100.L5.ovaries.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.L5.ovaries.by.girth), "deg100.L5.ovaries.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.ovaries.by.girth, res.L5.ovaries.by.girth, vsd.L5.ovaries.by.girth, deg100.L5.ovaries.by.girth,
     file = "L5.ovaries.by.girth.rda")
rm(meta.i, cts.i, dds.L5.ovaries.by.girth, res.L5.ovaries.by.girth, vsd.L5.ovaries.by.girth)
# load("L5.ovaries.by.girth.rda", verbose = TRUE)

###################
### L5 testes by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="L5", sex == "m", tissue != "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

dds.L5.testes.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + girth )
dds.L5.testes.by.girth <- DESeq(dds.L5.testes.by.girth)
dds.L5.testes.by.girth <- estimateSizeFactors(dds.L5.testes.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.L5.testes.by.girth)
res.L5.testes.by.girth <- results(dds.L5.testes.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.L5.testes.by.girth <- lfcShrink(dds.L5.testes.by.girth, coef="girth", type="apeglm", res = res.L5.testes.by.girth)
summary(res.L5.testes.by.girth)
# out of 24508 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 5, 0.02%
# LFC < 0 (down)     : 7, 0.029%
# outliers [1]       : 0, 0%

plotMA(res.L5.testes.by.girth, main="apeGLM", ylim=c(-3,4))
plot(plotCounts(dds.L5.testes.by.girth, gene=order(res.L5.testes.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.L5.testes.by.girth$log2FoldChange[order(res.L5.testes.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.L5.testes.by.girth <- vst(dds.L5.testes.by.girth, blind=FALSE)
plotPCA(vsd.L5.testes.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.L5.testes.by.girth <- res.L5.testes.by.girth[order(res.L5.testes.by.girth$padj),]
deg100.L5.testes.by.girth <- deg100.L5.testes.by.girth[1:100,]
deg100.L5.testes.by.girth <- apply.annotation(deg100.L5.testes.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.L5.testes.by.girth), "deg100.L5.testes.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.L5.testes.by.girth, res.L5.testes.by.girth, vsd.L5.testes.by.girth, deg100.L5.testes.by.girth,
     file = "L5.testes.by.girth.rda")
rm(meta.i, cts.i, dds.L5.testes.by.girth, res.L5.testes.by.girth, vsd.L5.testes.by.girth)
# load("L5.testes.by.girth.rda", verbose = TRUE)

###################
### Adult thorax by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue == "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

dds.adult.thorax.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph + girth )
dds.adult.thorax.by.girth <- DESeq(dds.adult.thorax.by.girth)
dds.adult.thorax.by.girth <- estimateSizeFactors(dds.adult.thorax.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.thorax.by.girth)
res.adult.thorax.by.girth <- results(dds.adult.thorax.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.adult.thorax.by.girth <- lfcShrink(dds.adult.thorax.by.girth, coef="girth", type="apeglm", res = res.adult.thorax.by.girth)
summary(res.adult.thorax.by.girth)
# out of 29858 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 46, 0.15%
# LFC < 0 (down)     : 90, 0.3%
# outliers [1]       : 0, 0%

plotMA(res.adult.thorax.by.girth, main="apeGLM", ylim=c(-5,5))
plot(plotCounts(dds.adult.thorax.by.girth, gene=order(res.adult.thorax.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.adult.thorax.by.girth$log2FoldChange[order(res.adult.thorax.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.thorax.by.girth <- vst(dds.adult.thorax.by.girth, blind=FALSE)
plotPCA(vsd.adult.thorax.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.adult.thorax.by.girth <- res.adult.thorax.by.girth[order(res.adult.thorax.by.girth$padj),]
deg100.adult.thorax.by.girth <- deg100.adult.thorax.by.girth[1:100,]
deg100.adult.thorax.by.girth <- apply.annotation(deg100.adult.thorax.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.adult.thorax.by.girth), "deg100.adult.thorax.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.thorax.by.girth, res.adult.thorax.by.girth, vsd.adult.thorax.by.girth, deg100.adult.thorax.by.girth,
     file = "adult.thorax.by.girth.rda")
rm(meta.i, cts.i, dds.adult.thorax.by.girth, res.adult.thorax.by.girth, vsd.adult.thorax.by.girth)
# load("adult.thorax.by.girth.rda", verbose = TRUE)


###################
### Adult gonad by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", tissue != "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

dds.adult.gonad.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + sex + morph + girth )
dds.adult.gonad.by.girth <- DESeq(dds.adult.gonad.by.girth)
dds.adult.gonad.by.girth <- estimateSizeFactors(dds.adult.gonad.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.gonad.by.girth)
res.adult.gonad.by.girth <- results(dds.adult.gonad.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.adult.gonad.by.girth <- lfcShrink(dds.adult.gonad.by.girth, coef="girth", type="apeglm", res = res.adult.gonad.by.girth)
summary(res.adult.gonad.by.girth)
# out of 49719 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 32, 0.064%
# LFC < 0 (down)     : 270, 0.54%
# outliers [1]       : 0, 0%

plotMA(res.adult.gonad.by.girth, main="apeGLM", ylim=c(-6,4))
plot(plotCounts(dds.adult.gonad.by.girth, gene=order(res.adult.gonad.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.adult.gonad.by.girth$log2FoldChange[order(res.adult.gonad.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.gonad.by.girth <- vst(dds.adult.gonad.by.girth, blind=FALSE)
plotPCA(vsd.adult.gonad.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.adult.gonad.by.girth <- res.adult.gonad.by.girth[order(res.adult.gonad.by.girth$padj),]
deg100.adult.gonad.by.girth <- deg100.adult.gonad.by.girth[1:100,]
deg100.adult.gonad.by.girth <- apply.annotation(deg100.adult.gonad.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.adult.gonad.by.girth), "deg100.adult.gonad.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.gonad.by.girth, res.adult.gonad.by.girth, vsd.adult.gonad.by.girth, deg100.adult.gonad.by.girth,
     file = "adult.gonad.by.girth.rda")
rm(meta.i, cts.i, dds.adult.gonad.by.girth, res.adult.gonad.by.girth, vsd.adult.gonad.by.girth)
# load("adult.gonad.by.girth.rda", verbose = TRUE)


###################
### Adult ovaries by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", sex == "f", tissue != "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

dds.adult.ovaries.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph + girth )
dds.adult.ovaries.by.girth <- DESeq(dds.adult.ovaries.by.girth)
dds.adult.ovaries.by.girth <- estimateSizeFactors(dds.adult.ovaries.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.ovaries.by.girth)
res.adult.ovaries.by.girth <- results(dds.adult.ovaries.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.adult.ovaries.by.girth <- lfcShrink(dds.adult.ovaries.by.girth, coef="girth", type="apeglm", res = res.adult.ovaries.by.girth)
summary(res.adult.ovaries.by.girth)
# out of 25956 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 47, 0.18%
# LFC < 0 (down)     : 404, 1.6%
# outliers [1]       : 0, 0%

plotMA(res.adult.ovaries.by.girth, main="apeGLM", ylim=c(-4,2))
plot(plotCounts(dds.adult.ovaries.by.girth, gene=order(res.adult.ovaries.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.adult.ovaries.by.girth$log2FoldChange[order(res.adult.ovaries.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.ovaries.by.girth <- vst(dds.adult.ovaries.by.girth, blind=FALSE)
plotPCA(vsd.adult.ovaries.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.adult.ovaries.by.girth <- res.adult.ovaries.by.girth[order(res.adult.ovaries.by.girth$padj),]
deg100.adult.ovaries.by.girth <- deg100.adult.ovaries.by.girth[1:100,]
deg100.adult.ovaries.by.girth <- apply.annotation(deg100.adult.ovaries.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.adult.ovaries.by.girth), "deg100.adult.ovaries.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.ovaries.by.girth, res.adult.ovaries.by.girth, vsd.adult.ovaries.by.girth, deg100.adult.ovaries.by.girth,
     file = "adult.ovaries.by.girth.rda")
rm(meta.i, cts.i, dds.adult.ovaries.by.girth, res.adult.ovaries.by.girth, vsd.adult.ovaries.by.girth)
# load("adult.ovaries.by.girth.rda", verbose = TRUE)

###################
### Adult testes by girth
###################

#### Sub-setting the data and metadata
meta.i <- sample.metadata %>%
  filter(stage=="adult", sex == "m", tissue != "thorax", !is.na(girth))
cts.i <- cts[,which(colnames(cts) %in% rownames(meta.i))]
cts.i <- as.matrix(cts.i[which(rowSums(cts.i) > gene.read.number.cut.off),])
dim(meta.i); dim(cts.i)
c(with(meta.i, by(morph_sex,morph_sex,length)))

dds.adult.testes.by.girth <- DESeqDataSetFromMatrix(
  countData=round(cts.i), colData=meta.i,
  design = ~ plate + morph + girth )
dds.adult.testes.by.girth <- DESeq(dds.adult.testes.by.girth)
dds.adult.testes.by.girth <- estimateSizeFactors(dds.adult.testes.by.girth)

# Results table with log2 fold changes, p values and adjusted p values.
resultsNames(dds.adult.testes.by.girth)
res.adult.testes.by.girth <- results(dds.adult.testes.by.girth, name=c("girth"), alpha=0.05, filterFun=ihw)
res.adult.testes.by.girth <- lfcShrink(dds.adult.testes.by.girth, coef="girth", type="apeglm", res = res.adult.testes.by.girth)
summary(res.adult.testes.by.girth)
# out of 35528 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1, 0.0028%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%

plotMA(res.adult.testes.by.girth, main="apeGLM", ylim=c(-1,4.5))
plot(plotCounts(dds.adult.testes.by.girth, gene=order(res.adult.testes.by.girth$padj)[1], intgroup="girth", returnData=TRUE))
res.adult.testes.by.girth$log2FoldChange[order(res.adult.testes.by.girth$padj)[1]]
# Positive LFC values are associated with higher girth expression bias

# Use variance stabilizing transformations only for ML applications
vsd.adult.testes.by.girth <- vst(dds.adult.testes.by.girth, blind=FALSE)
plotPCA(vsd.adult.testes.by.girth, intgroup=c("girth"))

# We can order our results table by the smallest p value:
deg100.adult.testes.by.girth <- res.adult.testes.by.girth[order(res.adult.testes.by.girth$padj),]
deg100.adult.testes.by.girth <- deg100.adult.testes.by.girth[1:100,]
deg100.adult.testes.by.girth <- apply.annotation(deg100.adult.testes.by.girth, ann)

# Save results
write.table(as.data.frame(deg100.adult.testes.by.girth), "deg100.adult.testes.by.girth.tsv", quote = FALSE, sep = "\t")
save(meta.i, cts.i, dds.adult.testes.by.girth, res.adult.testes.by.girth, vsd.adult.testes.by.girth, deg100.adult.testes.by.girth,
     file = "adult.testes.by.girth.rda")
rm(meta.i, cts.i, dds.adult.testes.by.girth, res.adult.testes.by.girth, vsd.adult.testes.by.girth)
# load("adult.testes.by.girth.rda", verbose = TRUE)


save(ann, apply.annotation, gene.read.number.cut.off, 
     res.adult.gonad.by.girth,res.adult.ovaries.by.girth,res.adult.testes.by.girth,res.adult.thorax.by.girth,res.L5.gonad.by.girth,res.L5.ovaries.by.girth,res.L5.testes.by.girth,res.L5.thorax.by.girth,
     file = "dge.analysis.res.girth.rda")

###################
# Volcano plots by girth
###################

{
  volcano.plots.by.girth <- list()
  volcano.plots.by.girth[[1]] <- 
    volcano.plot(res.L5.thorax.by.girth[,1:5], pCutoff = 0.05,
                 title = "Expression in juvenile thorax ~ batch + sex + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth[[2]] <- 
    volcano.plot(res.L5.gonad.by.girth[,1:5], pCutoff = 0.05,  
                 title = "Expression in juvenile gonads ~ batch + sex + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth[[3]] <- 
    volcano.plot(res.L5.ovaries.by.girth[,1:5], pCutoff = 0.05,  
                 title = "Expression in juvenile ovaries ~ batch + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth[[4]] <- 
    volcano.plot(res.L5.testes.by.girth[,1:5], pCutoff = 0.05, 
                 title = "Expression in juvenile testes ~ batch + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth[[5]] <- 
    volcano.plot(res.adult.thorax.by.girth[,1:5], pCutoff = 0.05, 
                 title = "Expression in adult thorax ~ batch + sex + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth[[6]] <- 
    volcano.plot(res.adult.gonad.by.girth[,1:5], pCutoff = 0.05,  
                 title = "Expression in adult gonads ~ batch + sex + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth[[7]] <- 
    volcano.plot(res.adult.ovaries.by.girth[,1:5], pCutoff = 0.05,  
                 title = "Expression in adult ovaries ~ batch + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth[[8]] <- 
    volcano.plot(res.adult.testes.by.girth[,1:5], pCutoff = 0.05, 
                 title = "Expression in adult testes ~ batch + girth",
                 NegLFClabel = "narrower", PosLFClabel = "wider", n=1)
  volcano.plots.by.girth.arranged <- ggarrange(plotlist = volcano.plots.by.girth, ncol = 4, nrow = 2, labels="AUTO")
  ggsave("plots/volcano.plots.by.girth.jpg", plot = volcano.plots.by.girth.arranged, width = 13, height = 6.5, scale = 1.5)
  ggsave("plots/volcano.plots.by.girth.pdf", plot = volcano.plots.by.girth.arranged, width = 13, height = 6.5, scale = 1.5)
}

###################
# Enrichment
###################

enrichment.go.adult.gonad.by.girth <- enrichment.go(res.adult.gonad.by.girth, ann)
enrichment.go.adult.ovaries.by.girth <- enrichment.go(res.adult.ovaries.by.girth, ann)
enrichment.go.adult.testes.by.girth <- enrichment.go(res.adult.testes.by.girth, ann)
enrichment.go.adult.thorax.by.girth <- enrichment.go(res.adult.thorax.by.girth, ann)
enrichment.go.L5.gonad.by.girth <- enrichment.go(res.L5.gonad.by.girth, ann)
enrichment.go.L5.ovaries.by.girth <- enrichment.go(res.L5.ovaries.by.girth, ann)
enrichment.go.L5.testes.by.girth <- enrichment.go(res.L5.testes.by.girth, ann)
enrichment.go.L5.thorax.by.girth <- enrichment.go(res.L5.thorax.by.girth, ann)

enrichment.kegg.adult.gonad.by.girth <- enrichment.kegg(res.adult.gonad.by.girth, ann)
enrichment.kegg.adult.ovaries.by.girth <- enrichment.kegg(res.adult.ovaries.by.girth, ann)
enrichment.kegg.adult.testes.by.girth <- enrichment.kegg(res.adult.testes.by.girth, ann)
enrichment.kegg.adult.thorax.by.girth <- enrichment.kegg(res.adult.thorax.by.girth, ann)
enrichment.kegg.L5.gonad.by.girth <- enrichment.kegg(res.L5.gonad.by.girth, ann)
enrichment.kegg.L5.ovaries.by.girth <- enrichment.kegg(res.L5.ovaries.by.girth, ann)
enrichment.kegg.L5.testes.by.girth <- enrichment.kegg(res.L5.testes.by.girth, ann)
enrichment.kegg.L5.thorax.by.girth <- enrichment.kegg(res.L5.thorax.by.girth, ann)

enrichment.pd.adult.gonad.by.girth <- enrichment.pd(res.adult.gonad.by.girth, ann)
enrichment.pd.adult.ovaries.by.girth <- enrichment.pd(res.adult.ovaries.by.girth, ann)
enrichment.pd.adult.testes.by.girth <- enrichment.pd(res.adult.testes.by.girth, ann)
enrichment.pd.adult.thorax.by.girth <- enrichment.pd(res.adult.thorax.by.girth, ann)
enrichment.pd.L5.gonad.by.girth <- enrichment.pd(res.L5.gonad.by.girth, ann)
enrichment.pd.L5.ovaries.by.girth <- enrichment.pd(res.L5.ovaries.by.girth, ann)
enrichment.pd.L5.testes.by.girth <- enrichment.pd(res.L5.testes.by.girth, ann)
enrichment.pd.L5.thorax.by.girth <- enrichment.pd(res.L5.thorax.by.girth, ann)

enrichment.xenic.adult.gonad.by.girth <- enrichment.xenic(res.adult.gonad.by.girth, ann)
enrichment.xenic.adult.ovaries.by.girth <- enrichment.xenic(res.adult.ovaries.by.girth, ann)
enrichment.xenic.adult.testes.by.girth <- enrichment.xenic(res.adult.testes.by.girth, ann)
enrichment.xenic.adult.thorax.by.girth <- enrichment.xenic(res.adult.thorax.by.girth, ann)
enrichment.xenic.L5.gonad.by.girth <- enrichment.xenic(res.L5.gonad.by.girth, ann)
enrichment.xenic.L5.ovaries.by.girth <- enrichment.xenic(res.L5.ovaries.by.girth, ann)
enrichment.xenic.L5.testes.by.girth <- enrichment.xenic(res.L5.testes.by.girth, ann)
enrichment.xenic.L5.thorax.by.girth <- enrichment.xenic(res.L5.thorax.by.girth, ann)

save( # Error with enrichment.go.adult.testes.by.girth
  enrichment.go.adult.gonad.by.girth, enrichment.go.adult.ovaries.by.girth, enrichment.go.adult.thorax.by.girth, enrichment.go.L5.gonad.by.girth, enrichment.go.L5.ovaries.by.girth, enrichment.go.L5.testes.by.girth, enrichment.go.L5.thorax.by.girth,
  enrichment.kegg.adult.gonad.by.girth, enrichment.kegg.adult.ovaries.by.girth, enrichment.kegg.adult.testes.by.girth, enrichment.kegg.adult.thorax.by.girth, enrichment.kegg.L5.gonad.by.girth, enrichment.kegg.L5.ovaries.by.girth, enrichment.kegg.L5.testes.by.girth, enrichment.kegg.L5.thorax.by.girth,
  enrichment.pd.adult.gonad.by.girth, enrichment.pd.adult.ovaries.by.girth, enrichment.pd.adult.testes.by.girth, enrichment.pd.adult.thorax.by.girth, enrichment.pd.L5.gonad.by.girth, enrichment.pd.L5.ovaries.by.girth, enrichment.pd.L5.testes.by.girth, enrichment.pd.L5.thorax.by.girth,
  enrichment.xenic.adult.gonad.by.girth, enrichment.xenic.adult.ovaries.by.girth, enrichment.xenic.adult.testes.by.girth, enrichment.xenic.adult.thorax.by.girth, enrichment.xenic.L5.gonad.by.girth, enrichment.xenic.L5.ovaries.by.girth, enrichment.xenic.L5.testes.by.girth, enrichment.xenic.L5.thorax.by.girth,
  file="enrichment.pvalues.girth.rda"
)

objs <- ls()[which(grepl("^enrichment\\.go\\.",ls()) & grepl("girth$",ls()))]
x <- unlist(lapply(objs, function (x) {
  meta.s <- sub("enrichment\\.","",x)
  meta.s <- sub("by\\.","by ",meta.s)
  meta.s <- paste0(unlist(str_split(meta.s,"\\.")), collapse = "\t")
  fdr.x <- get(x)
  s <- paste0(meta.s,"\t",top.go.terms(fdr.x) )
  s <- sub("\\\nnegative",paste0("\n",meta.s,"\tnegative"),s)
  return(s)
}) )
# Output to a text file
write.table(x,"enrichment.descriptions.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE,
            append = TRUE)

objs <- ls()[which(grepl("^enrichment\\.kegg\\.",ls()) & grepl("girth$",ls()))]
x <- unlist(lapply(objs, function (x) {
  meta.s <- sub("enrichment\\.","",x)
  meta.s <- sub("by\\.","by ",meta.s)
  meta.s <- paste0(unlist(str_split(meta.s,"\\.")), collapse = "\t")
  fdr.x <- get(x)
  if (is.null(fdr.x)) { 
    s <- paste0(meta.s,"\tpositive\t0\n",meta.s,"\tnegative\t0" ) 
  } else {
    s <- paste0(meta.s,"\t",top.kegg.terms(fdr.x) )
    s <- sub("\\\nnegative",paste0("\n",meta.s,"\tnegative"),s)
  }
  return(s)
}) )
# Output to a text file
write.table(x,"enrichment.descriptions.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE,
            append = TRUE)

objs <- ls()[which(grepl("^enrichment\\.pd\\.",ls()) & grepl("girth$",ls()))]
x <- unlist(lapply(objs, function (x) {
  meta.s <- sub("enrichment\\.","",x)
  meta.s <- sub("by\\.","by ",meta.s)
  meta.s <- paste0(unlist(str_split(meta.s,"\\.")), collapse = "\t")
  fdr.x <- get(x)
  if (is.null(fdr.x)) { 
    s <- paste0(meta.s,"\tpositive\t0\n",meta.s,"\tnegative\t0" ) 
  } else {
    s <- paste0(meta.s,"\t",top.pd.terms(fdr.x) )
    s <- sub("\\\nnegative",paste0("\n",meta.s,"\tnegative"),s)
  }
  return(s)
}) )
# Output to a text file
write.table(x,"enrichment.descriptions.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE,
            append = TRUE)

# Are contaminant or xenic transcripts enriched in any datasets? -- No
objs <- ls()[which(grepl("^enrichment\\.xenic\\.",ls()) & grepl("girth$",ls()))]
any(unlist(lapply(objs, function (x) {
  fdr.x <- get(x)
  return(any(unlist(fdr.x[,4:5]) < 0.05))
}) ))

# There's a significant enrichment of xenic transcripts in
# adult.gonad.by.girth *
# adult.ovaries.by.girth *
# adult.testes.by.girth *
# all associated with negative abdominal girth
# I spot check one of these transcripts with the greatest baseMean
# TR75328|c0_g2_i1 is a perfect BLASTn hit for AP013028 Wolbachia endosymbiont of Cimex lectularius 
# and it appears to encode the 30S ribosomal protein S11 
