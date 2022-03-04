rm(list=ls())
gc(verbose = TRUE)

library(dplyr)
library(magrittr)
require(abind)

###################
# Import the Kallisto gene count matrix
###################
cts <- read.csv("kallisto.counts.csv")
dim(cts)
# 395125    370
head(colnames(cts))
tail(colnames(cts))

sum(is.na(cts))
# 0

# Range of read sums by sample
summary(colSums(cts))
barplot(colSums(cts))
abline(h=summary(colSums(cts))[c(2,5)], col="darkred", lty=3)

# Range of counts for transcripts
summary(rowSums(cts))

# Read sums by transcript
hist(log10(rowSums(cts)+1))

# Coefficient of variation (CV) in transcript counts by sample,
# for those with at least 100 counts
x <- cts[which(rowSums(cts)>=100),]
cv <- function (x) { sd(x) / mean(x) }
hist(apply(x, 1, cv))

# Sample metadata
sample.metadata <- read.csv("3seq.metadata.csv", stringsAsFactors = FALSE)
# Metadata rows much match order of samples (columns) in the count matrix
length(unique(colnames(cts))); length(unique(sample.metadata$kallisto_name))
all(colnames(cts) %in% sample.metadata$kallisto_name)

sample.metadata <- sample.metadata %>%
  column_to_rownames(var = "kallisto_name") %>%
  select(class, population, stage, sex, morph, morph_sex, tissue, food_regime,
         seeds, cohort, plate, well, seq_date, note) %>%
  mutate(plate = as.factor(plate), seq_date = as.factor(seq_date)) %>%
  mutate(tissue = fct_relevel(tissue, "thorax")) %>%
  mutate(morph = fct_relevel(morph, "SW", "LW"))

# Examine read totals by sample, grouped by sequencing run
boxplot(colSums(cts) ~ sample.metadata$seq_date)

###################
# Combine the two sequencing runs (200611 and 200701) for plate 4 samples
###################

# Subset the count and metadata tables for the duplicate samples
x <- sample.metadata %>%
  filter(seq_date %in% c(200611,200701)) %>%
  rownames()
x <- which(colnames(cts) %in% x)
cts.i <- cts[,x]
metadata.i <- sample.metadata %>%
  filter(seq_date %in% c(200611,200701))

# Inspect the results
dim(cts.i); dim(metadata.i)

head(metadata.i)
head(colnames(cts.i))

# Sample one area of the matrix
cts.i[,which(metadata.i$seq_date==200611)][300000:300010,1:3]
cts.i[,which(metadata.i$seq_date==200701)][300000:300010,1:3]

# Reorganize the count matrix into an array with the duplicates in the thrid dimension
a <- abind::abind(
  cts.i[,which(metadata.i$seq_date==200611)],
  cts.i[,which(metadata.i$seq_date==200701)],
  along = 3
)
dim(a)

# Sum the array along its third (duplicate) dimension
# Warning: This takes a few minutes!
cts.j <- apply(a,c(1,2),sum)
dim(cts.j)
# Inspect that same area of the matrix and confirm that the values are summed.
cts.j[300000:300010,1:3]

# Subset the original count table to include the samples that were not sequenced twice
x <- sample.metadata %>%
  filter(seq_date %in% c(200611,200701)) %>%
  rownames()
x <- which(colnames(cts) %in% x)
cts.i <- cts[,-x]

# Overwrite the count table with a new one, combining counts for all samples
cts <- cbind(cts.i,cts.j)
colnames(cts) <- sub("_UMI_S\\d+","",colnames(cts))
cts <- cts[,order(colnames(cts))]
dim(cts.i); dim(cts.j); dim(cts)

# Exclude the redundant entries from the metadata
sample.metadata <- sample.metadata %>% 
  filter(seq_date!=200701) %>% 
  mutate(seq_date = droplevels(seq_date))
dim(sample.metadata)
rownames(sample.metadata) <- sub("_UMI_S\\d+","",rownames(sample.metadata))
sample.metadata <- sample.metadata[order(rownames(sample.metadata)),]
all(rownames(sample.metadata) == colnames(cts))

# Save the resulting count table and metadata
save(cts, sample.metadata, file = "combined.counts.rda")
# load("combined.counts.rda", verbose = TRUE)

# Clean up memory
rm(cts.i, cts.j, a, x)
gc(verbose = TRUE)

