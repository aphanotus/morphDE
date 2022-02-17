# Merge counts from the abundance.tsv files produced by Kallisto for
# individual specimens into a single tabular (CSV) file.
#
# Produces two files based on the "est_count" and "TPM" columns.
# Run from the command line as...
# bash$  R -e "source('merge.kallisto.results.R')"

start.time <- Sys.time()

# Set paths
# Be sure folders end in '/'
output.folder <- "./"
cts.file <- "kallisto.counts.csv"
tpm.file <- "kallisto.tpm.csv"
table.sep <- ","

pathbase <- getwd()
abundance.file.base.name <- ".abundance.tsv"

# Load packages
require("dplyr", quietly = TRUE, warn.conflicts = FALSE)
require("stringr", quietly = TRUE, warn.conflicts = FALSE)
require("magrittr", quietly = TRUE, warn.conflicts = FALSE)

# Get specimen IDs
abundance.files <- list.files(pathbase, pattern = abundance.file.base.name)
spec.names <- sub(abundance.file.base.name,"",abundance.files)

cat("Processing:\t")

# Initialize the abundance tables
k <- read.delim(abundance.files[1], stringsAsFactors = FALSE) %>%
  select(target_id, est_counts, tpm) %>%
  arrange(target_id)
gene.list <- k$target_id
counts <- data.frame(
  matrix(nrow = length(gene.list),
         ncol = length(spec.names)),
  row.names = gene.list)
colnames(counts) <- spec.names
tpm <- counts

number.of.samples <- dim(tpm)[2]

cat(number.of.samples,"samples\t",dim(tpm)[1],"features\n")

cat("1 of",number.of.samples,"\t",spec.names[1],"\n")
counts[,spec.names[1]] <- k$est_counts
tpm[,spec.names[1]] <- k$tpm

# Import gene-level abundance
for (i in 2:length(spec.names)) {
  cat(i,"of",number.of.samples,"\t",spec.names[i],"\n")
  k <- read.delim(abundance.files[i], stringsAsFactors = FALSE) %>%
    select(target_id, est_counts, tpm) %>%
    arrange(target_id)
  counts[,spec.names[i]] <- k$est_counts
  tpm[,spec.names[i]] <- k$tpm
}

if (any(is.na(counts))) { cat("WARNING: NA values in the count table!\n") }
if (any(is.na(tpm))) { cat("WARNING: NA values in the TPM table!\n") }

# Write out the resulting matrix as a CSV
cts.file <- paste0(output.folder, cts.file)
tpm.file <- paste0(output.folder, tpm.file)
write.table(counts, file=cts.file, quote=FALSE, sep=table.sep, row.names=TRUE, col.names = TRUE)
write.table(tpm, file=tpm.file, quote=FALSE, sep=table.sep, row.names=TRUE, col.names = TRUE)

cat("Done.\n",dim(tpm)[2],"samples\n",dim(tpm)[1],"features\n",sum(counts, na.rm = TRUE),"reads mapped\n")

x <- file.info(cts.file)
x <- utils:::format.object_size(x$size, "auto")
cat(cts.file,"\t",x,"\n")
x <- file.info(tpm.file)
x <- utils:::format.object_size(x$size, "auto")
cat(tpm.file,"\t",x,"\n")

elapsed.time <- Sys.time()-start.time
print(elapsed.time)

# End of R Script
