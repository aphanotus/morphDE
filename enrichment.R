library(stringr)

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

# See:
# https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/

enrichment.go <- function(res, ann, padj.cutoff = 0.05) {
  require(dplyr)
  require(tibble)
  # Annotate the results table
  res <- as.data.frame(res)
  res <- apply.annotation(res, ann)
  # Internal function to perform the hypergeometric test
  HGT <- function (res, up.regulated = TRUE) {
    if (up.regulated) { res <- res[which(res$log2FoldChange > 0),]}
    else {            res <- res[which(res$log2FoldChange < 0),]}
    # number of DEGs
    k <- sum(res$padj < padj.cutoff, na.rm = TRUE)       
    # Find list of GO terms applied to DEGs
    go.cols <- grep("GO",colnames(res))
    s <- paste0(unlist(res[which(res$padj<padj.cutoff),go.cols]), collapse = ",")
    s <- unlist(str_split(s,","))
    go.terms <- sort(unique(s))
    go.terms <- go.terms[which(nchar(go.terms)==15 | nchar(go.terms)==16)]
    # Mark the characters searchable
    go.terms <- sub("\\:","\\\\\\:",go.terms)
    go.terms <- sub("\\(","\\\\\\(",go.terms)
    go.terms <- sub("\\)","\\\\\\)",go.terms)
    # Iterate for each GO term
    p.values <- unlist(lapply(go.terms, function (go.i) {
      # DEGs with the GO term
      x <- sum(apply(res[which(res$padj<padj.cutoff),go.cols],1, function(x) any(grepl(go.i,x)) ), na.rm = TRUE) 
      # Genes with the GO term
      m <- sum(apply(res[,go.cols],1, function(x) any(grepl(go.i,x)) ), na.rm = TRUE) 
      # Genes without the GO term
      n <- dim(res)[1] - m    
      # Calculate the probability of enrichment (one-sided, lower-tail p-value)
      p <- phyper(x-1, m, n, k, lower.tail= FALSE)
      cat(ifelse(up.regulated,"up: ","down: "),gsub("\\\\","",go.i),"\t",p,"\n")
      return(p)
    })) # End lapply
    names(p.values) <- go.terms
    return(p.values)
  } # end function HGT
  enr.up.p <- HGT(res, up.regulated = TRUE)
  enr.down.p <- HGT(res, up.regulated = FALSE)
  df <- full_join(enframe(enr.up.p),enframe(enr.down.p), by="name")
  colnames(df) <- c("term","positive.DEGs.p","negative.DEGs.p") 
  df$positive.DEGs.fdr <- p.adjust(df$positive.DEGs.p, method = "fdr")
  df$negative.DEGs.fdr <- p.adjust(df$negative.DEGs.p, method = "fdr")
  df$term <- gsub("\\\\","",df$term)
  return(df)
} # end function

enrichment.kegg <- function(res, ann, padj.cutoff = 0.05) {
  require(dplyr)
  require(tibble)
  # Annotate the results table
  res <- as.data.frame(res)
  res <- apply.annotation(res, ann)
  # Internal function to perform the hypergeometric test
  HGT <- function (res, up.regulated = TRUE) {
    if (up.regulated) { res <- res[which(res$log2FoldChange > 0),]}
    else {            res <- res[which(res$log2FoldChange < 0),]}
    # number of DEGs
    k <- sum(res$padj < padj.cutoff, na.rm = TRUE)       
    # Find list of KEGG terms applied to DEGs
    kegg.col <- grep("EggNOG.KEGG.Terms",colnames(res))
    s <- paste0(unlist(res[which(res$padj<padj.cutoff),kegg.col]), collapse = ",")
    s <- gsub("NA\\,","",s)
    s <- unlist(str_split(s,","))
    s <- s[which(nchar(s)>0)]
    kegg.terms <- sort(unique(s))
    # Iterate for each KEGG term
    p.values <- unlist(lapply(kegg.terms, function (kegg.i) {
      s.i <- paste0("map",str_pad(kegg.i,width=5,side="left","0"))
      kegg.i <- paste0("[\\^\\,]",kegg.i,"[\\,\\$]")
      # DEGs with the KEGG term
      x <- sum((res$padj < padj.cutoff) & grepl(kegg.i,res[,kegg.col]), na.rm = TRUE) 
      # Genes with the KEGG term
      m <- sum(grepl(kegg.i,res[,kegg.col]), na.rm = TRUE) 
      # Genes without the KEGG term
      n <- dim(res)[1] - m    
      # Calculate the probability of enrichment (one-sided, lower-tail p-value)
      p <- phyper(x-1, m, n, k, lower.tail= FALSE)
      cat(ifelse(up.regulated,"up: ","down: "),s.i,"\t",p,"\n")
      return(p)
    })) # End lapply
    if (!is.null(kegg.terms) & !is.null(p.values)) { names(p.values) <- kegg.terms }
    return(p.values)
  } # end function HGT
  enr.up.p <- HGT(res, up.regulated = TRUE)
  enr.down.p <- HGT(res, up.regulated = FALSE)
  if (is.null(enr.up.p) & is.null(enr.down.p)) {
    warning("No DEGs - Returning NULL")
    return(NULL)
  }
  if (!is.null(enr.up.p) & !is.null(enr.down.p)) {
    df <- full_join(enframe(enr.up.p),enframe(enr.down.p), by="name")
    colnames(df) <- c("term","positive.DEGs.p","negative.DEGs.p") 
    df$positive.DEGs.fdr <- p.adjust(df$positive.DEGs.p, method = "fdr")
    df$negative.DEGs.fdr <- p.adjust(df$negative.DEGs.p, method = "fdr")
    df$term <- paste0("map",str_pad(df$term,width=5,side="left","0"))
    return(df)
  }
  if (!is.null(enr.up.p)) {
    df <- data_frame(
      term = paste0("map",str_pad(names(enr.up.p),width=5,side="left","0")),
      positive.DEGs.p = enr.up.p,
      negative.DEGs.p = rep(NA, length(enr.up.p)),
      positive.DEGs.fdr = p.adjust(enr.up.p, method = "fdr"),
      negative.DEGs.fdr = rep(NA, length(enr.up.p)),
    ) 
  } else {
    df <- data_frame(
      term = paste0("map",str_pad(names(enr.down.p),width=5,side="left","0")),
      positive.DEGs.p = rep(NA, length(enr.down.p)),
      negative.DEGs.p = enr.down.p,
      positive.DEGs.fdr = rep(NA, length(enr.down.p)),
      negative.DEGs.fdr = p.adjust(enr.down.p, method = "fdr"),
    ) 
  }
  return(df)
} # end function

enrichment.pd <- function(res, ann, padj.cutoff = 0.05) {
  require(dplyr)
  require(tibble)
  # Annotate the results table
  res <- as.data.frame(res)
  res <- apply.annotation(res, ann)
  clean.string <- function(s) {
    s <- gsub("PFAM \\(","",s)
    s <- gsub("SMART \\(","",s)
    s <- gsub("\\), SMART \\(",",",s)
    s <- gsub("\\)","",s)
    s <- gsub("\\, ",",",s)
    s <- gsub("\\,\\,",",",s)
    s <- gsub("\\,\\,",",",s)
    s <- gsub("\\,\\,",",",s)
    s <- gsub("\\,\\,",",",s)
    s <- unlist(str_split(s,","))
    s <- s[which(nchar(s)>0)]
    return(s)
  }
  # Internal function to perform the hypergeometric test
  HGT <- function (res, up.regulated = TRUE) {
    if (up.regulated) { res <- res[which(res$log2FoldChange > 0),]}
    else {            res <- res[which(res$log2FoldChange < 0),]}
    # number of DEGs
    k <- sum(res$padj < padj.cutoff, na.rm = TRUE)       
    # Find list of protein domains applied to DEGs
    pd.col <- grep("EggNOG.Protein.Domains",colnames(res))
    s <- clean.string(paste0(unlist(res[which(res$padj<padj.cutoff),pd.col]), collapse = ","))
    pd.terms <- sort(unique(s))
    # Iterate for each KEGG term
    p.values <- unlist(lapply(pd.terms, function (pd.i) {
      s.i <- pd.i
      pd.i <- paste0("[ \\(]",pd.i,"[\\,\\)]")
      # DEGs with the KEGG term
      x <- sum((res$padj < padj.cutoff) & grepl(pd.i,res[,pd.col]), na.rm = TRUE) 
      # Genes with the EggNOG Protein Domain
      m <- sum(grepl(pd.i,res[,pd.col]), na.rm = TRUE) 
      # Genes without the KEGG term
      n <- dim(res)[1] - m    
      # Calculate the probability of enrichment (one-sided, lower-tail p-value)
      p <- phyper(x-1, m, n, k, lower.tail= FALSE)
      cat(ifelse(up.regulated,"up: ","down: "),s.i,"\t",p,"\n")
      return(p)
    })) # End lapply
    if (!is.null(pd.terms) & !is.null(p.values)) { names(p.values) <- pd.terms }
    return(p.values)
  } # end function HGT
  enr.up.p <- HGT(res, up.regulated = TRUE)
  enr.down.p <- HGT(res, up.regulated = FALSE)
  if (is.null(enr.up.p) & is.null(enr.down.p)) {
    warning("No DEGs - Returning NULL")
    return(NULL)
  }
  if (!is.null(enr.up.p) & !is.null(enr.down.p)) {
    df <- full_join(enframe(enr.up.p),enframe(enr.down.p), by="name")
    colnames(df) <- c("term","positive.DEGs.p","negative.DEGs.p") 
    df$positive.DEGs.fdr <- p.adjust(df$positive.DEGs.p, method = "fdr")
    df$negative.DEGs.fdr <- p.adjust(df$negative.DEGs.p, method = "fdr")
    return(df)
  }
  if (!is.null(enr.up.p)) {
    df <- data_frame(
      term = paste0("map",str_pad(names(enr.up.p),width=5,side="left","0")),
      positive.DEGs.p = enr.up.p,
      negative.DEGs.p = rep(NA, length(enr.up.p)),
      positive.DEGs.fdr = p.adjust(enr.up.p, method = "fdr"),
      negative.DEGs.fdr = rep(NA, length(enr.up.p)),
    ) 
  } else {
    df <- data_frame(
      term = paste0("map",str_pad(names(enr.down.p),width=5,side="left","0")),
      positive.DEGs.p = rep(NA, length(enr.down.p)),
      negative.DEGs.p = enr.down.p,
      positive.DEGs.fdr = rep(NA, length(enr.down.p)),
      negative.DEGs.fdr = p.adjust(enr.down.p, method = "fdr"),
    ) 
  }
  return(df)
} # end function

enrichment.xenic <- function(res, ann, padj.cutoff = 0.05) {
  require(dplyr)
  require(tibble)
  # Annotate the results table
  res <- as.data.frame(res)
  res <- apply.annotation(res, ann)
  
  df <- data_frame(
    term = c("Contaminant","xenic"),
    positive.DEGs.p = c(1,1),
    negative.DEGs.p = c(1,1)
  )
  
  # Positive DEGs
  res.i <- res[which(res$log2FoldChange > 0),]
  k <- sum(res.i$padj < padj.cutoff, na.rm = TRUE)       
  p.values <- c(NA,NA)
  names(p.values) <- c("Contaminant","xenic")
  con <- res.i$Contaminant=="Yes"
  # Contaminant DEGs 
  x <- sum((res.i$padj < padj.cutoff) & (con), na.rm = TRUE)  
  # Contaminant Genes 
  m <- sum(con, na.rm = TRUE) 
  # Non-Contaminant Genes
  n <- dim(res.i)[1] - m    
  p.values[1] <- phyper(x-1, m, n, k, lower.tail= FALSE)
  cat("up: Contaminant\t",p.values[1],"\n")
  # Xenic DEGs 
  x <- sum((res.i$padj < padj.cutoff) & (res.i$xenic), na.rm = TRUE) 
  # Contaminant Genes 
  m <- sum(res.i$xenic, na.rm = TRUE) 
  # Non-Contaminant Genes
  n <- dim(res.i)[1] - m    
  p.values[2] <- phyper(x-1, m, n, k, lower.tail= FALSE)
  cat("up: xenic\t",p.values[2],"\n")
  df$positive.DEGs.p <- p.values
  
  # Negative DEGs
  res.i <- res[which(res$log2FoldChange < 0),]
  k <- sum(res.i$padj < padj.cutoff, na.rm = TRUE)       
  p.values <- c(NA,NA)
  names(p.values) <- c("Contaminant","xenic")
  con <- res.i$Contaminant=="Yes"
  # Contaminant DEGs 
  x <- sum((res.i$padj < padj.cutoff) & (con), na.rm = TRUE)  
  # Contaminant Genes 
  m <- sum(con, na.rm = TRUE) 
  # Non-Contaminant Genes
  n <- dim(res.i)[1] - m    
  p.values[1] <- phyper(x-1, m, n, k, lower.tail= FALSE)
  cat("down: Contaminant\t",p.values[1],"\n")
  # Xenic DEGs 
  x <- sum((res.i$padj < padj.cutoff) & (res.i$xenic), na.rm = TRUE) 
  # Contaminant Genes 
  m <- sum(res.i$xenic, na.rm = TRUE) 
  # Non-Contaminant Genes
  n <- dim(res.i)[1] - m    
  p.values[2] <- phyper(x-1, m, n, k, lower.tail= FALSE)
  cat("down: xenic\t",p.values[2],"\n")
  df$negative.DEGs.p <- p.values
  
  df$positive.DEGs.fdr <- p.adjust(df$positive.DEGs.p, method = "fdr")
  df$negative.DEGs.fdr <- p.adjust(df$negative.DEGs.p, method = "fdr")
  return(df)
} # end function

load("dge.results.res.rda", verbose = TRUE)

{
  enrichment.go.adult.thorax.by.sex <- enrichment.go(res.adult.thorax.by.sex, ann)
  enrichment.go.adult.gonad.by.sex <- enrichment.go(res.adult.gonad.by.sex, ann)
  enrichment.go.adult.thorax.by.morph <- enrichment.go(res.adult.thorax.by.morph, ann)
  enrichment.go.adult.gonad.by.morph <- enrichment.go(res.adult.gonad.by.morph, ann)
  enrichment.go.adult.thorax.by.food <- enrichment.go(res.adult.thorax.by.food, ann)
  enrichment.go.adult.gonad.by.food <- enrichment.go(res.adult.gonad.by.food, ann)
  enrichment.go.L5.gonad.by.sex <- enrichment.go(res.L5.gonad.by.sex, ann)
  enrichment.go.L5.gonad.by.food <- enrichment.go(res.L5.gonad.by.food, ann)
  enrichment.go.L5.thorax.by.sex <- enrichment.go(res.L5.thorax.by.sex, ann)
  enrichment.go.L5.thorax.by.food <- enrichment.go(res.L5.thorax.by.food, ann)
  enrichment.go.ovaries.by.morph <- enrichment.go(res.ovaries.by.morph, ann)
  enrichment.go.ovaries.by.stage <- enrichment.go(res.ovaries.by.stage, ann)
  enrichment.go.testes.by.morph <- enrichment.go(res.testes.by.morph, ann)
  enrichment.go.testes.by.stage <- enrichment.go(res.testes.by.stage, ann)
  enrichment.go.thorax.by.stage <- enrichment.go(res.thorax.by.stage, ann)
  enrichment.go.adult.thorax.by.food_density <- enrichment.go(res.adult.thorax.by.food_density, ann)
  enrichment.go.adult.gonad.by.food_density <- enrichment.go(res.adult.gonad.by.food_density, ann)
  enrichment.go.adult.ovaries.by.food_density <- enrichment.go(res.adult.ovaries.by.food_density, ann)
  enrichment.go.adult.testes.by.food_density <- enrichment.go(res.adult.testes.by.food_density, ann)
  enrichment.go.L5.thorax.by.food_density <- enrichment.go(res.L5.thorax.by.food_density, ann)
  enrichment.go.L5.gonad.by.food_density <- enrichment.go(res.L5.gonad.by.food_density, ann)
  enrichment.go.L5.ovaries.by.food_density <- enrichment.go(res.L5.ovaries.by.food_density, ann)
  enrichment.go.L5.testes.by.food_density <- enrichment.go(res.L5.testes.by.food_density, ann)
  enrichment.go.adult.thorax.by.wingPC1 <- enrichment.go(res.adult.thorax.by.wingPC1, ann)
  enrichment.go.adult.gonad.by.wingPC1 <- enrichment.go(res.adult.gonad.by.wingPC1, ann)
  enrichment.go.ovaries.by.wingPC1 <- enrichment.go(res.ovaries.by.wingPC1, ann)
  enrichment.go.testes.by.wingPC1 <- enrichment.go(res.testes.by.wingPC1, ann)
  enrichment.go.adult.thorax.by.txPC1 <- enrichment.go(res.adult.thorax.by.txPC1, ann)
  enrichment.go.L5.gonad.by.wingpadPC1 <- enrichment.go(res.L5.gonad.by.wingpadPC1, ann)
}

{
  enrichment.kegg.adult.thorax.by.sex <- enrichment.kegg(res.adult.thorax.by.sex, ann)
  enrichment.kegg.adult.gonad.by.sex <- enrichment.kegg(res.adult.gonad.by.sex, ann)
  enrichment.kegg.adult.thorax.by.morph <- enrichment.kegg(res.adult.thorax.by.morph, ann)
  enrichment.kegg.adult.gonad.by.morph <- enrichment.kegg(res.adult.gonad.by.morph, ann)
  enrichment.kegg.adult.thorax.by.food <- enrichment.kegg(res.adult.thorax.by.food, ann)
  enrichment.kegg.adult.gonad.by.food <- enrichment.kegg(res.adult.gonad.by.food, ann)
  enrichment.kegg.L5.gonad.by.sex <- enrichment.kegg(res.L5.gonad.by.sex, ann)
  enrichment.kegg.L5.gonad.by.food <- enrichment.kegg(res.L5.gonad.by.food, ann)
  enrichment.kegg.L5.thorax.by.sex <- enrichment.kegg(res.L5.thorax.by.sex, ann)
  enrichment.kegg.L5.thorax.by.food <- enrichment.kegg(res.L5.thorax.by.food, ann)
  enrichment.kegg.ovaries.by.morph <- enrichment.kegg(res.ovaries.by.morph, ann)
  enrichment.kegg.ovaries.by.stage <- enrichment.kegg(res.ovaries.by.stage, ann)
  enrichment.kegg.testes.by.morph <- enrichment.kegg(res.testes.by.morph, ann)
  enrichment.kegg.testes.by.stage <- enrichment.kegg(res.testes.by.stage, ann)
  enrichment.kegg.thorax.by.stage <- enrichment.kegg(res.thorax.by.stage, ann)
  enrichment.kegg.adult.thorax.by.food_density <- enrichment.kegg(res.adult.thorax.by.food_density, ann)
  enrichment.kegg.adult.gonad.by.food_density <- enrichment.kegg(res.adult.gonad.by.food_density, ann)
  enrichment.kegg.adult.ovaries.by.food_density <- enrichment.kegg(res.adult.ovaries.by.food_density, ann)
  enrichment.kegg.adult.testes.by.food_density <- enrichment.kegg(res.adult.testes.by.food_density, ann)
  enrichment.kegg.L5.thorax.by.food_density <- enrichment.kegg(res.L5.thorax.by.food_density, ann)
  enrichment.kegg.L5.gonad.by.food_density <- enrichment.kegg(res.L5.gonad.by.food_density, ann)
  enrichment.kegg.L5.ovaries.by.food_density <- enrichment.kegg(res.L5.ovaries.by.food_density, ann)
  enrichment.kegg.L5.testes.by.food_density <- enrichment.kegg(res.L5.testes.by.food_density, ann)
  enrichment.kegg.adult.thorax.by.wingPC1 <- enrichment.kegg(res.adult.thorax.by.wingPC1, ann)
  enrichment.kegg.adult.gonad.by.wingPC1 <- enrichment.kegg(res.adult.gonad.by.wingPC1, ann)
  enrichment.kegg.ovaries.by.wingPC1 <- enrichment.kegg(res.ovaries.by.wingPC1, ann)
  enrichment.kegg.testes.by.wingPC1 <- enrichment.kegg(res.testes.by.wingPC1, ann)
  enrichment.kegg.adult.thorax.by.txPC1 <- enrichment.kegg(res.adult.thorax.by.txPC1, ann)
  enrichment.kegg.L5.gonad.by.wingpadPC1 <- enrichment.kegg(res.L5.gonad.by.wingpadPC1, ann)
}

{
  enrichment.pd.adult.thorax.by.sex <- enrichment.pd(res.adult.thorax.by.sex, ann)
  enrichment.pd.adult.gonad.by.sex <- enrichment.pd(res.adult.gonad.by.sex, ann)
  enrichment.pd.adult.thorax.by.morph <- enrichment.pd(res.adult.thorax.by.morph, ann)
  enrichment.pd.adult.gonad.by.morph <- enrichment.pd(res.adult.gonad.by.morph, ann)
  enrichment.pd.adult.thorax.by.food <- enrichment.pd(res.adult.thorax.by.food, ann)
  enrichment.pd.adult.gonad.by.food <- enrichment.pd(res.adult.gonad.by.food, ann)
  enrichment.pd.L5.gonad.by.sex <- enrichment.pd(res.L5.gonad.by.sex, ann)
  enrichment.pd.L5.gonad.by.food <- enrichment.pd(res.L5.gonad.by.food, ann)
  enrichment.pd.L5.thorax.by.sex <- enrichment.pd(res.L5.thorax.by.sex, ann)
  enrichment.pd.L5.thorax.by.food <- enrichment.pd(res.L5.thorax.by.food, ann)
  enrichment.pd.ovaries.by.morph <- enrichment.pd(res.ovaries.by.morph, ann)
  enrichment.pd.ovaries.by.stage <- enrichment.pd(res.ovaries.by.stage, ann)
  enrichment.pd.testes.by.morph <- enrichment.pd(res.testes.by.morph, ann)
  enrichment.pd.testes.by.stage <- enrichment.pd(res.testes.by.stage, ann)
  enrichment.pd.thorax.by.stage <- enrichment.pd(res.thorax.by.stage, ann)
  enrichment.pd.adult.thorax.by.food_density <- enrichment.pd(res.adult.thorax.by.food_density, ann)
  enrichment.pd.adult.gonad.by.food_density <- enrichment.pd(res.adult.gonad.by.food_density, ann)
  enrichment.pd.adult.ovaries.by.food_density <- enrichment.pd(res.adult.ovaries.by.food_density, ann)
  enrichment.pd.adult.testes.by.food_density <- enrichment.pd(res.adult.testes.by.food_density, ann)
  enrichment.pd.L5.thorax.by.food_density <- enrichment.pd(res.L5.thorax.by.food_density, ann)
  enrichment.pd.L5.gonad.by.food_density <- enrichment.pd(res.L5.gonad.by.food_density, ann)
  enrichment.pd.L5.ovaries.by.food_density <- enrichment.pd(res.L5.ovaries.by.food_density, ann)
  enrichment.pd.L5.testes.by.food_density <- enrichment.pd(res.L5.testes.by.food_density, ann)
  enrichment.pd.adult.thorax.by.wingPC1 <- enrichment.pd(res.adult.thorax.by.wingPC1, ann)
  enrichment.pd.adult.gonad.by.wingPC1 <- enrichment.pd(res.adult.gonad.by.wingPC1, ann)
  enrichment.pd.ovaries.by.wingPC1 <- enrichment.pd(res.ovaries.by.wingPC1, ann)
  enrichment.pd.testes.by.wingPC1 <- enrichment.pd(res.testes.by.wingPC1, ann)
  enrichment.pd.adult.thorax.by.txPC1 <- enrichment.pd(res.adult.thorax.by.txPC1, ann)
  enrichment.pd.L5.gonad.by.wingpadPC1 <- enrichment.pd(res.L5.gonad.by.wingpadPC1, ann)
}

{
  enrichment.xenic.adult.thorax.by.sex <- enrichment.xenic(res.adult.thorax.by.sex, ann)
  enrichment.xenic.adult.gonad.by.sex <- enrichment.xenic(res.adult.gonad.by.sex, ann)
  enrichment.xenic.adult.thorax.by.morph <- enrichment.xenic(res.adult.thorax.by.morph, ann)
  enrichment.xenic.adult.gonad.by.morph <- enrichment.xenic(res.adult.gonad.by.morph, ann)
  enrichment.xenic.adult.thorax.by.food <- enrichment.xenic(res.adult.thorax.by.food, ann)
  enrichment.xenic.adult.gonad.by.food <- enrichment.xenic(res.adult.gonad.by.food, ann)
  enrichment.xenic.L5.gonad.by.sex <- enrichment.xenic(res.L5.gonad.by.sex, ann)
  enrichment.xenic.L5.gonad.by.food <- enrichment.xenic(res.L5.gonad.by.food, ann)
  enrichment.xenic.L5.thorax.by.sex <- enrichment.xenic(res.L5.thorax.by.sex, ann)
  enrichment.xenic.L5.thorax.by.food <- enrichment.xenic(res.L5.thorax.by.food, ann)
  enrichment.xenic.ovaries.by.morph <- enrichment.xenic(res.ovaries.by.morph, ann)
  enrichment.xenic.ovaries.by.stage <- enrichment.xenic(res.ovaries.by.stage, ann)
  enrichment.xenic.testes.by.morph <- enrichment.xenic(res.testes.by.morph, ann)
  enrichment.xenic.testes.by.stage <- enrichment.xenic(res.testes.by.stage, ann)
  enrichment.xenic.thorax.by.stage <- enrichment.xenic(res.thorax.by.stage, ann)
  enrichment.xenic.adult.thorax.by.food_density <- enrichment.xenic(res.adult.thorax.by.food_density, ann)
  enrichment.xenic.adult.gonad.by.food_density <- enrichment.xenic(res.adult.gonad.by.food_density, ann)
  enrichment.xenic.adult.ovaries.by.food_density <- enrichment.xenic(res.adult.ovaries.by.food_density, ann)
  enrichment.xenic.adult.testes.by.food_density <- enrichment.xenic(res.adult.testes.by.food_density, ann)
  enrichment.xenic.L5.thorax.by.food_density <- enrichment.xenic(res.L5.thorax.by.food_density, ann)
  enrichment.xenic.L5.gonad.by.food_density <- enrichment.xenic(res.L5.gonad.by.food_density, ann)
  enrichment.xenic.L5.ovaries.by.food_density <- enrichment.xenic(res.L5.ovaries.by.food_density, ann)
  enrichment.xenic.L5.testes.by.food_density <- enrichment.xenic(res.L5.testes.by.food_density, ann)
  enrichment.xenic.adult.thorax.by.wingPC1 <- enrichment.xenic(res.adult.thorax.by.wingPC1, ann)
  enrichment.xenic.adult.gonad.by.wingPC1 <- enrichment.xenic(res.adult.gonad.by.wingPC1, ann)
  enrichment.xenic.ovaries.by.wingPC1 <- enrichment.xenic(res.ovaries.by.wingPC1, ann)
  enrichment.xenic.testes.by.wingPC1 <- enrichment.xenic(res.testes.by.wingPC1, ann)
  enrichment.xenic.adult.thorax.by.txPC1 <- enrichment.xenic(res.adult.thorax.by.txPC1, ann)
  enrichment.xenic.L5.gonad.by.wingpadPC1 <- enrichment.xenic(res.L5.gonad.by.wingpadPC1, ann)
}

{
  save(
    enrichment.go, 
    enrichment.go.adult.gonad.by.food, 
    enrichment.go.adult.gonad.by.food_density, 
    enrichment.go.adult.gonad.by.morph, 
    enrichment.go.adult.gonad.by.sex, 
    enrichment.go.adult.gonad.by.wingPC1, 
    enrichment.go.adult.ovaries.by.food_density, 
    enrichment.go.adult.testes.by.food_density, 
    enrichment.go.adult.thorax.by.food, 
    enrichment.go.adult.thorax.by.food_density, 
    enrichment.go.adult.thorax.by.morph, 
    enrichment.go.adult.thorax.by.sex, 
    enrichment.go.adult.thorax.by.txPC1, 
    enrichment.go.adult.thorax.by.wingPC1, 
    enrichment.go.L5.gonad.by.food, 
    enrichment.go.L5.gonad.by.food_density, 
    enrichment.go.L5.gonad.by.sex, 
    enrichment.go.L5.gonad.by.wingpadPC1, 
    enrichment.go.L5.ovaries.by.food_density, 
    enrichment.go.L5.testes.by.food_density, 
    enrichment.go.L5.thorax.by.food, 
    enrichment.go.L5.thorax.by.food_density, 
    enrichment.go.L5.thorax.by.sex, 
    enrichment.go.ovaries.by.morph, 
    enrichment.go.ovaries.by.stage, 
    enrichment.go.ovaries.by.wingPC1, 
    enrichment.go.testes.by.morph, 
    enrichment.go.testes.by.stage, 
    enrichment.go.testes.by.wingPC1, 
    enrichment.go.thorax.by.stage, 
    enrichment.kegg, 
    enrichment.kegg.adult.gonad.by.food, 
    enrichment.kegg.adult.gonad.by.food_density, 
    enrichment.kegg.adult.gonad.by.morph, 
    enrichment.kegg.adult.gonad.by.sex, 
    enrichment.kegg.adult.gonad.by.wingPC1, 
    enrichment.kegg.adult.ovaries.by.food_density, 
    enrichment.kegg.adult.testes.by.food_density, 
    enrichment.kegg.adult.thorax.by.food, 
    enrichment.kegg.adult.thorax.by.food_density, 
    enrichment.kegg.adult.thorax.by.morph, 
    enrichment.kegg.adult.thorax.by.sex, 
    enrichment.kegg.adult.thorax.by.txPC1, 
    enrichment.kegg.adult.thorax.by.wingPC1, 
    enrichment.kegg.L5.gonad.by.food, 
    enrichment.kegg.L5.gonad.by.food_density, 
    enrichment.kegg.L5.gonad.by.sex, 
    enrichment.kegg.L5.gonad.by.wingpadPC1, 
    enrichment.kegg.L5.ovaries.by.food_density, 
    enrichment.kegg.L5.testes.by.food_density, 
    enrichment.kegg.L5.thorax.by.food, 
    enrichment.kegg.L5.thorax.by.food_density, 
    enrichment.kegg.L5.thorax.by.sex, 
    enrichment.kegg.ovaries.by.morph, 
    enrichment.kegg.ovaries.by.stage, 
    enrichment.kegg.ovaries.by.wingPC1, 
    enrichment.kegg.testes.by.morph, 
    enrichment.kegg.testes.by.stage, 
    enrichment.kegg.testes.by.wingPC1, 
    enrichment.kegg.thorax.by.stage, 
    enrichment.pd, 
    enrichment.pd.adult.gonad.by.food, 
    enrichment.pd.adult.gonad.by.food_density, 
    enrichment.pd.adult.gonad.by.morph, 
    enrichment.pd.adult.gonad.by.sex, 
    enrichment.pd.adult.gonad.by.wingPC1, 
    enrichment.pd.adult.ovaries.by.food_density, 
    enrichment.pd.adult.testes.by.food_density, 
    enrichment.pd.adult.thorax.by.food, 
    enrichment.pd.adult.thorax.by.food_density, 
    enrichment.pd.adult.thorax.by.morph, 
    enrichment.pd.adult.thorax.by.sex, 
    enrichment.pd.adult.thorax.by.txPC1, 
    enrichment.pd.adult.thorax.by.wingPC1, 
    enrichment.pd.L5.gonad.by.food, 
    enrichment.pd.L5.gonad.by.food_density, 
    enrichment.pd.L5.gonad.by.sex, 
    enrichment.pd.L5.gonad.by.wingpadPC1, 
    enrichment.pd.L5.ovaries.by.food_density, 
    enrichment.pd.L5.testes.by.food_density, 
    enrichment.pd.L5.thorax.by.food, 
    enrichment.pd.L5.thorax.by.food_density, 
    enrichment.pd.L5.thorax.by.sex, 
    enrichment.pd.ovaries.by.morph, 
    enrichment.pd.ovaries.by.stage, 
    enrichment.pd.ovaries.by.wingPC1, 
    enrichment.pd.testes.by.morph, 
    enrichment.pd.testes.by.stage, 
    enrichment.pd.testes.by.wingPC1, 
    enrichment.pd.thorax.by.stage, 
    enrichment.xenic, 
    enrichment.xenic.adult.gonad.by.food, 
    enrichment.xenic.adult.gonad.by.food_density, 
    enrichment.xenic.adult.gonad.by.morph, 
    enrichment.xenic.adult.gonad.by.sex, 
    enrichment.xenic.adult.gonad.by.wingPC1, 
    enrichment.xenic.adult.ovaries.by.food_density, 
    enrichment.xenic.adult.testes.by.food_density, 
    enrichment.xenic.adult.thorax.by.food, 
    enrichment.xenic.adult.thorax.by.food_density, 
    enrichment.xenic.adult.thorax.by.morph, 
    enrichment.xenic.adult.thorax.by.sex, 
    enrichment.xenic.adult.thorax.by.txPC1, 
    enrichment.xenic.adult.thorax.by.wingPC1, 
    enrichment.xenic.L5.gonad.by.food, 
    enrichment.xenic.L5.gonad.by.food_density, 
    enrichment.xenic.L5.gonad.by.sex, 
    enrichment.xenic.L5.gonad.by.wingpadPC1, 
    enrichment.xenic.L5.ovaries.by.food_density, 
    enrichment.xenic.L5.testes.by.food_density, 
    enrichment.xenic.L5.thorax.by.food, 
    enrichment.xenic.L5.thorax.by.food_density, 
    enrichment.xenic.L5.thorax.by.sex, 
    enrichment.xenic.ovaries.by.morph, 
    enrichment.xenic.ovaries.by.stage, 
    enrichment.xenic.ovaries.by.wingPC1, 
    enrichment.xenic.testes.by.morph, 
    enrichment.xenic.testes.by.stage, 
    enrichment.xenic.testes.by.wingPC1, 
    enrichment.xenic.thorax.by.stage, 
    file="enrichment.pvalues.rda"
  )
}

load("enrichment.pvalues.rda", verbose = TRUE)

# Find the top GO term descriptions via API 
top.go.terms <- function(x, top = 10) {
  require(GO.db)
  require(stringr)
  n.sig.pos <- sum(x$positive.DEGs.fdr < 0.05, na.rm = TRUE)
  n.sig.neg <- sum(x$negative.DEGs.fdr < 0.05, na.rm = TRUE)
  # Positive DEG enrichment
  x.i <- x[order(x$positive.DEGs.fdr),c(1,4)]
  # Exclude level 0 and level 1 terms
  x.i <- x.i[which(!grepl("L=[01]\\)$",x.i$term)),]
  x.i$Levels <- sub("\\)","",str_split_fixed(x.i$term,"=",2)[,2])
  x.i$term <- strtrim(x.i$term,10)
  x.i <- x.i[which(!is.na(GOID(x.i$term))),]
  x.i <- x.i[c(1:top),]
  x.i <- x.i[which(x.i$positive.DEGs.fdr<0.05),]
  x.i$description <- Term(x.i$term)
  # Record output string
  s <- with(x.i, paste0(term, " ", description, ", ", signif(positive.DEGs.fdr,2), collapse = "; ") )
  # Negative DEG enrichment
  x.i <- x[order(x$negative.DEGs.fdr),c(1,5)]
  # Exclude level 0 and level 1 terms
  x.i <- x.i[which(!grepl("L=[01]\\)$",x.i$term)),]
  x.i$Levels <- sub("\\)","",str_split_fixed(x.i$term,"=",2)[,2])
  x.i$term <- strtrim(x.i$term,10)
  x.i <- x.i[which(!is.na(GOID(x.i$term))),]
  x.i <- x.i[c(1:top),]
  x.i <- x.i[which(x.i$negative.DEGs.fdr<0.05),]
  x.i$description <- Term(x.i$term)
  # Record output string
  s <- paste0("postive\t",n.sig.pos,"\t",s,"\nnegative\t",n.sig.neg,"\t")
  s <- paste0(s,with(x.i, paste0(term, " ", description, ", ", signif(negative.DEGs.fdr,2), collapse = "; ") ) )
  return(s)
} # End function

x <- unlist(lapply(ls()[grep("^enrichment\\.go\\.",ls())], function (x) {
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
            append = FALSE)

# Find the top KEGG term descriptions via API 
top.kegg.terms <- function(x, top = 10) {
  require(RCurl)
  require(stringr)
  n.sig.pos <- sum(x$positive.DEGs.fdr < 0.05, na.rm = TRUE)
  n.sig.neg <- sum(x$negative.DEGs.fdr < 0.05, na.rm = TRUE)
  # Positive DEG enrichment
  x.i <- x[order(x$positive.DEGs.fdr),c(1,4)]
  x.i <- x.i[c(1:top),]
  x.i <- x.i[which(x.i$positive.DEGs.fdr<0.05),]
  x.i$description <- getURL(paste0("http://togows.dbcls.jp/entry/pathway/",x.i$term,"/name"))
  x.i$description <- sub("\n$","",x.i$description)
  # Record output string
  s <- with(x.i, paste0(term, " ", description, ", ", signif(positive.DEGs.fdr,2), collapse = "; ") )
  # Negative DEG enrichment
  x.i <- x[order(x$negative.DEGs.fdr),c(1,5)]
  x.i <- x.i[c(1:top),]
  x.i <- x.i[which(x.i$negative.DEGs.fdr<0.05),]
  x.i$description <- getURL(paste0("http://togows.dbcls.jp/entry/pathway/",x.i$term,"/name"))
  x.i$description <- sub("\n$","",x.i$description)
  # Record output string
  s <- paste0("postive\t",n.sig.pos,"\t",s,"\nnegative\t",n.sig.neg,"\t")
  s <- paste0(s,with(x.i, paste0(term, " ", description, ", ", signif(negative.DEGs.fdr,2), collapse = "; ") ) )
  return(s)
} # End function

x <- unlist(lapply(ls()[grep("^enrichment\\.kegg\\.",ls())], function (x) {
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


top.pd.terms <- function(x, top = 10) {
  require(stringr)
  x$term <- sub("^map0000","",x$term)
  x$term <- sub("^map000","",x$term)
  x$term <- sub("^map00","",x$term)
  x$term <- sub("^map0","",x$term)
  x$term <- sub("^map","",x$term)
  n.sig.pos <- sum(x$positive.DEGs.fdr < 0.05, na.rm = TRUE)
  n.sig.neg <- sum(x$negative.DEGs.fdr < 0.05, na.rm = TRUE)
  # Positive DEG enrichment
  x.i <- x[order(x$positive.DEGs.fdr),c(1,4)]
  x.i <- x.i[c(1:top),]
  x.i <- x.i[which(x.i$positive.DEGs.fdr<0.05),]
  # Record output string
  s <- with(x.i, paste0(term, ", ", signif(positive.DEGs.fdr,2), collapse = "; ") )
  # Negative DEG enrichment
  x.i <- x[order(x$negative.DEGs.fdr),c(1,5)]
  x.i <- x.i[c(1:top),]
  x.i <- x.i[which(x.i$negative.DEGs.fdr<0.05),]
  # Record output string
  s <- paste0("postive\t",n.sig.pos,"\t",s,"\nnegative\t",n.sig.neg,"\t")
  s <- paste0(s,with(x.i, paste0(term, ", ", signif(negative.DEGs.fdr,2), collapse = "; ") ) )
  return(s)
} # End function

x <- unlist(lapply(ls()[grep("^enrichment\\.pd\\.",ls())], function (x) {
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
any(unlist(lapply(ls()[grep("^enrichment\\.xenic\\.",ls())], function (x) {
  fdr.x <- get(x)
  return(any(unlist(fdr.x[,4:5]) < 0.05))
}) ))

