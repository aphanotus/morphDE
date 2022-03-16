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

enrichment.go <- function(res, ann, padj.cutoff = 0.05) {
  # See:
  # https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
  # https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
  res <- as.data.frame(res)
  # Annotate table
  res <- apply.annotation(res, ann)
  # number of DEGs
  k <- sum(res$padj < padj.cutoff, na.rm = TRUE)       
  
  if (k<1) {
    warning(paste0("No genes are differentially expressed at a padj cutoff of ",padj.cutoff,"\n"))
    return(NA)
  }
  
  # Find list of GO terms
  # Filter to terms applied to DEGs
  go.cols <- grep("GO",colnames(res))
  s <- paste0(unlist(res[which(res$padj<padj.cutoff),go.cols]), collapse = ",")
  s <- unlist(str_split(s,","))
  # hist(log10(by(s,s,length)))
  go.terms <- sort(unique(s))
  go.terms <- go.terms[which(nchar(go.terms)==15 | nchar(go.terms)==16)]
  # length(go.terms)
  
  p.values <- rep(NA,length(go.terms))
  names(p.values) <- go.terms
  for (i in 1:length(go.terms)) {
    go.i <- go.terms[i]
    go.i <- sub("\\:","\\\\\\:",go.i)
    go.i <- sub("\\(","\\\\\\(",go.i)
    go.i <- sub("\\)","\\\\\\)",go.i)
    
    # DEGs with the GO term
    x <- sum(apply(res[which(res$padj<padj.cutoff),go.cols],1, function(x) any(grepl(go.i,x)) ), na.rm = TRUE) 
    # Genes with the GO term
    m <- sum(apply(res[,go.cols],1, function(x) any(grepl(go.i,x)) ), na.rm = TRUE) 
    # Genes without the GO term
    n <- dim(res)[1] - m    
    
    # Calculate the probability of enrichment (one-sided, lower-tail p-value)
    p.values[i] <- phyper(x-1, m, n, k, lower.tail= FALSE)
    # phyper is faster than fisher.test
    # fisher.test(matrix(c(x, m-x, k-x, n-(k-x)), 2, 2), alternative='greater')$p.value
    print(p.values[i])
  } # end for loop
  return(p.values)
} # end function

enrichment.kegg <- function(res, ann, padj.cutoff = 0.05) {
  res <- as.data.frame(res)
  # Annotate table
  res <- apply.annotation(res, ann)
  # number of DEGs
  k <- sum(res$padj < padj.cutoff, na.rm = TRUE)       
  
  if (k<1) {
    warning(paste0("No genes are differentially expressed at a padj cutoff of ",padj.cutoff,"\n"))
    return(NA)
  }
  
  # Find list of KEGG terms
  # Filter to terms applied to DEGs
  kegg.col <- grep("EggNOG.KEGG.Terms",colnames(res))
  s <- paste0(unlist(res[which(res$padj<padj.cutoff),kegg.col]), collapse = ",")
  s <- gsub("NA\\,","",s)
  s <- unlist(str_split(s,","))
  s <- s[which(nchar(s)>0)]
  kegg.terms <- sort(unique(s))
  
  p.values <- rep(NA,length(kegg.terms))
  names(p.values) <- kegg.terms
  for (i in 1:length(kegg.terms)) {
    kegg.i <- kegg.terms[i]
    kegg.i <- paste0("[\\^\\,]",kegg.i,"[\\,\\$]")
    # DEGs with the KEGG term
    x <- sum((res$padj < padj.cutoff) & grepl(kegg.i,res[,kegg.col]), na.rm = TRUE) 
    # Genes with the KEGG term
    m <- sum(grepl(kegg.i,res[,kegg.col]), na.rm = TRUE) 
    # Genes without the KEGG term
    n <- dim(res)[1] - m    
    
    # Calculate the probability of enrichment (one-sided, lower-tail p-value)
    p.values[i] <- phyper(x-1, m, n, k, lower.tail= FALSE)
    print(p.values[i])
  } # end for loop
  return(p.values)
} # end function

enrichment.pd <- function(res, ann, padj.cutoff = 0.05) {
  res <- as.data.frame(res)
  # Annotate table
  res <- apply.annotation(res, ann)
  # number of DEGs
  k <- sum(res$padj < padj.cutoff, na.rm = TRUE)       
  
  if (k<1) {
    warning(paste0("No genes are differentially expressed at a padj cutoff of ",padj.cutoff,"\n"))
    return(NA)
  }
  
  clean.string <- function(s) {
    s <- gsub("PFAM \\(","",s)
    s <- gsub("SMART \\(","",s)
    s <- gsub("\\), SMART \\(",",",s)
    s <- gsub("\\)","",s)
    s <- gsub("\\, ",",",s)
    s <- gsub("\\,\\,",",",s)
    s <- gsub("\\,\\,",",",s)
    s <- gsub("\\,\\,",",",s)
    s <- unlist(str_split(s,","))
    s <- s[which(nchar(s)>0)]
    return(s)
  }
  
  # Find list of EggNOG Protein Domains
  # Filter to terms applied to DEGs
  pd.col <- grep("EggNOG.Protein.Domains",colnames(res))
  s <- clean.string(paste0(unlist(res[which(res$padj<padj.cutoff),pd.col]), collapse = ","))
  pd.terms <- sort(unique(s))
  
  p.values <- rep(NA,length(pd.terms))
  names(p.values) <- pd.terms
  for (i in 1:length(pd.terms)) {
    pd.i <- pd.terms[i]
    pd.i <- paste0("[ \\(]",pd.i,"[\\,\\)]")
    # DEGs with the EggNOG Protein Domain
    x <- sum((res$padj < padj.cutoff) & grepl(pd.i,res[,pd.col]), na.rm = TRUE) 
    # Genes with the EggNOG Protein Domain
    m <- sum(grepl(pd.i,res[,pd.col]), na.rm = TRUE) 
    # Genes without the EggNOG Protein Domain
    n <- dim(res)[1] - m    
    
    # Calculate the probability of enrichment (one-sided, lower-tail p-value)
    p.values[i] <- phyper(x-1, m, n, k, lower.tail= FALSE)
    print(p.values[i])
  } # end for loop
  return(p.values)
} # end function

enrichment.xenic <- function(res, ann, padj.cutoff = 0.05) {
  res <- as.data.frame(res)
  # Annotate table
  res <- apply.annotation(res, ann)
  # number of DEGs
  k <- sum(res$padj < padj.cutoff, na.rm = TRUE)       
  
  if (k<1) {
    warning(paste0("No genes are differentially expressed at a padj cutoff of ",padj.cutoff,"\n"))
    return(NA)
  }
  
  p.values <- c(NA,NA)
  names(p.values) <- c("Contaminant","xenic")
  con <- res$Contaminant=="Yes"
  
  # Contaminant DEGs 
  x <- sum((res$padj < padj.cutoff) & (con), na.rm = TRUE) 
  # Contaminant Genes 
  m <- sum(con, na.rm = TRUE) 
  # Non-Contaminant Genes
  n <- dim(res)[1] - m    
  p.values[1] <- phyper(x-1, m, n, k, lower.tail= FALSE)
  print(p.values[1])
  
  # Xenic DEGs 
  x <- sum((res$padj < padj.cutoff) & (res$xenic), na.rm = TRUE) 
  # Contaminant Genes 
  m <- sum(res$xenic, na.rm = TRUE) 
  # Non-Contaminant Genes
  n <- dim(res)[1] - m    
  p.values[2] <- phyper(x-1, m, n, k, lower.tail= FALSE)
  print(p.values[2])
  
  return(p.values)
} # end function

load("dge.results.dfs.rda", verbose = TRUE)

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

save(
  enrichment.go.adult.thorax.by.sex, enrichment.go.adult.gonad.by.sex, enrichment.go.adult.thorax.by.morph, enrichment.go.adult.gonad.by.morph, enrichment.go.adult.thorax.by.food, enrichment.go.adult.gonad.by.food, enrichment.go.L5.gonad.by.sex, enrichment.go.L5.gonad.by.food, enrichment.go.L5.thorax.by.sex, enrichment.go.L5.thorax.by.food,
  enrichment.kegg.adult.thorax.by.sex, enrichment.kegg.adult.gonad.by.sex, enrichment.kegg.adult.thorax.by.morph, enrichment.kegg.adult.gonad.by.morph, enrichment.kegg.adult.thorax.by.food, enrichment.kegg.adult.gonad.by.food, enrichment.kegg.L5.gonad.by.sex, enrichment.kegg.L5.gonad.by.food, enrichment.kegg.L5.thorax.by.sex, enrichment.kegg.L5.thorax.by.food,
  enrichment.pd.adult.thorax.by.sex, enrichment.pd.adult.gonad.by.sex, enrichment.pd.adult.thorax.by.morph, enrichment.pd.adult.gonad.by.morph, enrichment.pd.adult.thorax.by.food, enrichment.pd.adult.gonad.by.food, enrichment.pd.L5.gonad.by.sex, enrichment.pd.L5.gonad.by.food, enrichment.pd.L5.thorax.by.sex, enrichment.pd.L5.thorax.by.food,
  enrichment.xenic.adult.thorax.by.sex, enrichment.xenic.adult.gonad.by.sex, enrichment.xenic.adult.thorax.by.morph, enrichment.xenic.adult.gonad.by.morph, enrichment.xenic.adult.thorax.by.food, enrichment.xenic.adult.gonad.by.food, enrichment.xenic.L5.gonad.by.sex, enrichment.xenic.L5.gonad.by.food, enrichment.xenic.L5.thorax.by.sex, enrichment.xenic.L5.thorax.by.food,
  enrichment.go.ovaries.by.morph, enrichment.go.ovaries.by.stage, enrichment.go.testes.by.morph, enrichment.go.testes.by.stage, enrichment.go.thorax.by.stage, enrichment.kegg.ovaries.by.morph, enrichment.kegg.ovaries.by.stage, enrichment.kegg.testes.by.morph, enrichment.kegg.testes.by.stage, enrichment.kegg.thorax.by.stage, 
  enrichment.pd.ovaries.by.morph, enrichment.pd.ovaries.by.stage, enrichment.pd.testes.by.morph, enrichment.pd.testes.by.stage, enrichment.pd.thorax.by.stage, enrichment.xenic.ovaries.by.morph, enrichment.xenic.ovaries.by.stage, enrichment.xenic.testes.by.morph, enrichment.xenic.testes.by.stage, enrichment.xenic.thorax.by.stage,
  file="enrichment.pvalues.rda"
)

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

save(
  enrichment.go.adult.thorax.by.food_density, enrichment.go.adult.gonad.by.food_density, enrichment.go.adult.ovaries.by.food_density, enrichment.go.adult.testes.by.food_density, enrichment.go.L5.thorax.by.food_density, enrichment.go.L5.gonad.by.food_density, enrichment.go.L5.ovaries.by.food_density, enrichment.go.L5.testes.by.food_density, enrichment.go.adult.thorax.by.wingPC1, enrichment.go.adult.gonad.by.wingPC1, enrichment.go.ovaries.by.wingPC1, enrichment.go.testes.by.wingPC1, enrichment.go.adult.thorax.by.txPC1, enrichment.go.L5.gonad.by.wingpadPC1,
  enrichment.kegg.adult.thorax.by.food_density, enrichment.kegg.adult.gonad.by.food_density, enrichment.kegg.adult.ovaries.by.food_density, enrichment.kegg.adult.testes.by.food_density, enrichment.kegg.L5.thorax.by.food_density, enrichment.kegg.L5.gonad.by.food_density, enrichment.kegg.L5.ovaries.by.food_density, enrichment.kegg.L5.testes.by.food_density, enrichment.kegg.adult.thorax.by.wingPC1, enrichment.kegg.adult.gonad.by.wingPC1, enrichment.kegg.ovaries.by.wingPC1, enrichment.kegg.testes.by.wingPC1, enrichment.kegg.adult.thorax.by.txPC1, enrichment.kegg.L5.gonad.by.wingpadPC1,
  enrichment.pd.adult.thorax.by.food_density, enrichment.pd.adult.gonad.by.food_density, enrichment.pd.adult.ovaries.by.food_density, enrichment.pd.adult.testes.by.food_density, enrichment.pd.L5.thorax.by.food_density, enrichment.pd.L5.gonad.by.food_density, enrichment.pd.L5.ovaries.by.food_density, enrichment.pd.L5.testes.by.food_density, enrichment.pd.adult.thorax.by.wingPC1, enrichment.pd.adult.gonad.by.wingPC1, enrichment.pd.ovaries.by.wingPC1, enrichment.pd.testes.by.wingPC1, enrichment.pd.adult.thorax.by.txPC1, enrichment.pd.L5.gonad.by.wingpadPC1,
  enrichment.xenic.adult.thorax.by.food_density, enrichment.xenic.adult.gonad.by.food_density, enrichment.xenic.adult.ovaries.by.food_density, enrichment.xenic.adult.testes.by.food_density, enrichment.xenic.L5.thorax.by.food_density, enrichment.xenic.L5.gonad.by.food_density, enrichment.xenic.L5.ovaries.by.food_density, enrichment.xenic.L5.testes.by.food_density, enrichment.xenic.adult.thorax.by.wingPC1, enrichment.xenic.adult.gonad.by.wingPC1, enrichment.xenic.ovaries.by.wingPC1, enrichment.xenic.testes.by.wingPC1, enrichment.xenic.adult.thorax.by.txPC1, enrichment.xenic.L5.gonad.by.wingpadPC1, 
  file="enrichment.pvalues.part2.rda"
)

# Number of significant terms
x <- unlist(lapply(ls()[grep("^enrichment\\.go",ls())], function (x) {
  p.x <- get(x)
  p.x <- p.x[!grepl("L=[01]",names(p.x))]
  sig <- sum(p.adjust(p.x, method = "fdr") < 0.05, na.rm = TRUE)
  names(sig) <- x
  return(sig)
}) )
cat(names(x))
cat(unname(x))

# Find the top GO term descriptions via API 
top.go.terms <- function(x, top = 10, return.string = TRUE) {
  require(GO.db)
  require(stringr)
  
  ids <- names(sort(p.adjust(x, method = "fdr")))
  
  # Exclude level 0 and level 1 terms
  ids <- ids[!grepl("L=0",ids)]
  ids <- ids[!grepl("L=1",ids)]
  
  ids <- strtrim(ids,10)
  ids <- ids[which(!is.na(GOID(ids)))]
  ids <- ids[1:top]
  
  long.ids <- unlist(lapply(ids, function(i) { names(x)[grep(i,names(x))] }))
  Levels <- sub("\\)","",str_split_fixed(long.ids,"=",2)[,2])
  p.values <- unlist(lapply(ids, function(i) { x[grep(i,names(x))] }))
  FDR <- unlist(lapply(ids, function(i) { p.adjust(x, method = "fdr")[grep(i,names(x))] }))
  GOTerms <- Term(ids)
  
  if (return.string) {
    s <- paste0(long.ids, " ", GOTerms, ", ", signif(FDR,4), collapse = "; ")
    s <- gsub("\\de-","×10^-",s)
    return(s)
  } else {
    df <- data.frame(
      id = ids,
      level = Levels,
      p = p.values,
      FDR,
      term = GOTerms
    )
    rownames(df) <- 1:top
    return(df)
  }
} # End function

ls()[grep("enrichment\\.go",ls())]
top.go.terms(enrichment.go.adult.gonad.by.morph)
top.go.terms(enrichment.go.adult.thorax.by.morph)
top.go.terms(enrichment.go.adult.thorax.by.sex)
top.go.terms(enrichment.go.L5.gonad.by.sex)
top.go.terms(enrichment.go.L5.thorax.by.food)
top.go.terms(enrichment.go.ovaries.by.morph)
top.go.terms(enrichment.go.ovaries.by.stage)
top.go.terms(enrichment.go.testes.by.morph)
top.go.terms(enrichment.go.thorax.by.stage)
top.go.terms(enrichment.go.adult.thorax.by.food_density)
top.go.terms(enrichment.go.adult.gonad.by.food_density)
top.go.terms(enrichment.go.adult.gonad.by.wingPC1)
top.go.terms(enrichment.go.adult.ovaries.by.food_density)
top.go.terms(enrichment.go.adult.testes.by.food_density)
top.go.terms(enrichment.go.adult.thorax.by.wingPC1)
top.go.terms(enrichment.go.L5.gonad.by.food_density)
top.go.terms(enrichment.go.L5.ovaries.by.food_density)
top.go.terms(enrichment.go.ovaries.by.wingPC1)
top.go.terms(enrichment.go.testes.by.wingPC1)

# Finding descriptions of the KEGG and EggNOG terms must be done manually
ls()[grep("enrichment\\.kegg",ls())]
# https://www.genome.jp/kegg/kegg2.html

x <- unlist(lapply(ls()[grep("^enrichment\\.kegg",ls())], function (x) {
  sig <- sum(p.adjust(get(x), method = "fdr") < 0.05, na.rm = TRUE)
  names(sig) <- x
  return(sig)
}) )
cat(names(x))
cat(unname(x))

head(sort(p.adjust(enrichment.kegg.adult.gonad.by.food, method = "fdr")),2)
head(sort(p.adjust(enrichment.kegg.adult.gonad.by.morph, method = "fdr")),10)
head(sort(p.adjust(enrichment.kegg.adult.thorax.by.morph, method = "fdr")),10)
head(sort(p.adjust(enrichment.kegg.adult.thorax.by.sex, method = "fdr")),3)
head(sort(p.adjust(enrichment.kegg.L5.gonad.by.sex, method = "fdr")),10)
head(sort(p.adjust(enrichment.kegg.ovaries.by.morph, method = "fdr")),10)
head(sort(p.adjust(enrichment.kegg.testes.by.morph, method = "fdr")),7)
head(sort(p.adjust(enrichment.kegg.thorax.by.stage, method = "fdr")),10)
head(sort(p.adjust(enrichment.kegg.adult.gonad.by.food_density, method = "fdr")),2)
head(sort(p.adjust(enrichment.kegg.adult.gonad.by.wingPC1, method = "fdr")),5)
head(sort(p.adjust(enrichment.kegg.adult.testes.by.food_density, method = "fdr")),4)
head(sort(p.adjust(enrichment.kegg.adult.thorax.by.wingPC1, method = "fdr")),10)
head(sort(p.adjust(enrichment.kegg.L5.gonad.by.food_density, method = "fdr")),10)
head(sort(p.adjust(enrichment.kegg.L5.gonad.by.wingpadPC1, method = "fdr")),3)

ls()[grep("enrichment\\.pd",ls())]
# https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml

x <- unlist(lapply(ls()[grep("^enrichment\\.pd",ls())], function (x) {
  sig <- sum(p.adjust(get(x), method = "fdr") < 0.05, na.rm = TRUE)
  names(sig) <- x
  return(sig)
}) )
cat(names(x))
cat(unname(x))

head(sort(p.adjust(enrichment.pd.adult.gonad.by.food, method = "fdr")),8)
head(sort(p.adjust(enrichment.pd.adult.gonad.by.morph, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.thorax.by.food, method = "fdr")),5)
head(sort(p.adjust(enrichment.pd.adult.thorax.by.morph, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.thorax.by.sex, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.L5.gonad.by.food, method = "fdr")),4)
head(sort(p.adjust(enrichment.pd.L5.gonad.by.sex, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.L5.thorax.by.food, method = "fdr")),7)
head(sort(p.adjust(enrichment.pd.L5.thorax.by.sex, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.ovaries.by.morph, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.ovaries.by.stage, method = "fdr")),5)
head(sort(p.adjust(enrichment.pd.testes.by.morph, method = "fdr")),7)
head(sort(p.adjust(enrichment.pd.thorax.by.stage, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.gonad.by.food_density, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.gonad.by.wingPC1, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.ovaries.by.food_density, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.testes.by.food_density, method = "fdr")),9)
head(sort(p.adjust(enrichment.pd.adult.thorax.by.food_density, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.thorax.by.txPC1, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.adult.thorax.by.wingPC1, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.L5.gonad.by.food_density, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.L5.gonad.by.wingpadPC1, method = "fdr")),1)
head(sort(p.adjust(enrichment.pd.L5.ovaries.by.food_density, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.L5.thorax.by.food_density, method = "fdr")),10)
head(sort(p.adjust(enrichment.pd.ovaries.by.wingPC1, method = "fdr")),6)
head(sort(p.adjust(enrichment.pd.testes.by.wingPC1, method = "fdr")),10)

unlist(lapply(ls()[grep("^enrichment\\.xenic",ls())], function (x) {
  sig <- sum(p.adjust(get(x), method = "fdr") < 0.05, na.rm = TRUE)
  names(sig) <- x
  return(sig)
}) )

