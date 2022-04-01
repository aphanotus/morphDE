load("model.output/adult.thorax.by.morph.rda", verbose = TRUE)
meta.adult.thorax.by.morph <- meta.i
load("model.output/adult.gonad.by.morph.rda", verbose = TRUE)
meta.adult.gonad.by.morph <- meta.i
load("model.output/ovaries.by.morph.rda", verbose = TRUE)
meta.adult.ovaries.by.morph <- meta.i
load("model.output/testes.by.morph.rda", verbose = TRUE)
meta.adult.testes.by.morph <- meta.i

ctn.adult.thorax.by.morph <- counts(dds.adult.thorax.by.morph, normalized=TRUE)
ctn.adult.gonad.by.morph <- counts(dds.adult.gonad.by.morph, normalized=TRUE)
ctn.adult.ovaries.by.morph <- counts(dds.adult.ovaries.by.morph, normalized=TRUE)
ctn.adult.testes.by.morph <- counts(dds.adult.testes.by.morph, normalized=TRUE)

save(ctn.adult.thorax.by.morph, ctn.adult.gonad.by.morph, ctn.adult.ovaries.by.morph, ctn.adult.testes.by.morph,
     meta.adult.thorax.by.morph, meta.adult.gonad.by.morph, meta.adult.ovaries.by.morph, meta.adult.testes.by.morph,
     file = "model.output/ctn.adult.tissues.by.morph.rda")

rm(dds.adult.thorax.by.morph, dds.adult.gonad.by.morph, dds.adult.ovaries.by.morph, dds.adult.testes.by.morph)

goi.by.morph.plots <- function (transcript.name, transcript.id, panel.labels=LETTERS[1:4]) {
  plot.list <- list()
  # Thorax
  title.i <- paste0(transcript.name," in thorax")
  ids.i <- transcript.id[which(transcript.id %in% rownames(ctn.adult.thorax.by.morph))]
  if (length(ids.i)<1) {
    title.i <- paste0("      ",title.i)
    plot.list[[1]] <- data.frame(x = c(1:3)-2, y = c(1:3)-2) %>% 
      ggplot() + theme_void() +
      ggtitle(title.i) + 
      annotate(geom="text", x=0, y=0, label="No expression") 
  } else {
    if (length(ids.i)>1) {
      cts <- colSums(ctn.adult.thorax.by.morph[ids.i,])
    } else {
      cts <- ctn.adult.thorax.by.morph[ids.i,]
    }
    p.val <- signif(t.test(cts ~ meta.adult.thorax.by.morph$morph, var.equal = FALSE)$p.value,2)
    p.val <- p.adjust(p.val, method = "fdr", n=4)
    if (p.val < 0.05) { title.i <- paste0(title.i," (p = ",p.val,")") }
    plot.list[[1]] <- data.frame(
      x = meta.adult.thorax.by.morph$morph_sex,
      y = cts
    ) %>% 
      mutate(x = stringr::str_replace_all(x, c(
        "LWf" = "long-winged\nfemale", "LWm" = "long-winged\nmale",
        "SWf" = "short-winged\nfemale", "SWm" = "short-winged\nmale"))) %>% 
      mutate(x = as.factor(x)) %>% 
      ggplot(aes(x,y)) +
      theme_bw() +
      geom_violin(color = "gray50") + 
      geom_boxplot(color = "gray50", fill = "gray90", width = 0.2, outlier.shape = NA, coef = 0) +
      geom_jitter(width = 0.2, alpha = 0.75) +
      scale_y_continuous(limits=c(0,NA)) +
      labs(x=NULL,y="normalized transcript counts", title = title.i)
  }
  
  # Gonads
  title.i <- paste0("gonads")
  ids.i <- transcript.id[which(transcript.id %in% rownames(ctn.adult.gonad.by.morph))]
  if (length(ids.i)<1) {
    title.i <- paste0("      ",title.i)
    plot.list[[2]] <- data.frame(x = c(1:3)-2, y = c(1:3)-2) %>% 
      ggplot() + theme_void() +
      ggtitle(title.i) + 
      annotate(geom="text", x=0, y=0, label="No expression") 
  } else {
    if (length(ids.i)>1) {
      cts <- colSums(ctn.adult.gonad.by.morph[ids.i,])
    } else {
      cts <- ctn.adult.gonad.by.morph[ids.i,]
    }
    p.val <- signif(t.test(cts ~ meta.adult.gonad.by.morph$morph, var.equal = FALSE)$p.value,2)
    p.val <- p.adjust(p.val, method = "fdr", n=4)
    if (p.val < 0.05) { title.i <- paste0(title.i," (p = ",p.val,")") }
    plot.list[[2]] <- data.frame(
      x = meta.adult.gonad.by.morph$morph_sex,
      y = cts
    ) %>% 
      mutate(x = stringr::str_replace_all(x, c(
        "LWf" = "long-winged\nfemale", "LWm" = "long-winged\nmale",
        "SWf" = "short-winged\nfemale", "SWm" = "short-winged\nmale"))) %>% 
      mutate(x = as.factor(x)) %>% 
      ggplot(aes(x,y)) +
      theme_bw() +
      geom_violin(color = "gray50") + 
      geom_boxplot(color = "gray50", fill = "gray90", width = 0.2, outlier.shape = NA, coef = 0) +
      geom_jitter(width = 0.2, alpha = 0.75) +
      scale_y_continuous(limits=c(0,NA)) +
      labs(x=NULL,y="normalized transcript counts", title = title.i)
  }
  
  # Ovaries
  title.i <- paste0("ovaries")
  ids.i <- transcript.id[which(transcript.id %in% rownames(ctn.adult.ovaries.by.morph))]
  if (length(ids.i)<1) {
    title.i <- paste0("      ",title.i)
    plot.list[[3]] <- data.frame(x = c(1:3)-2, y = c(1:3)-2) %>% 
      ggplot() + theme_void() +
      ggtitle(title.i) + 
      annotate(geom="text", x=0, y=0, label="No expression") 
  } else {
    if (length(ids.i)>1) {
      cts <- colSums(ctn.adult.ovaries.by.morph[ids.i,])
    } else {
      cts <- ctn.adult.ovaries.by.morph[ids.i,]
    }
    p.val <- signif(t.test(cts ~ meta.adult.ovaries.by.morph$morph, var.equal = FALSE)$p.value,2)
    p.val <- p.adjust(p.val, method = "fdr", n=4)
    if (p.val < 0.05) { title.i <- paste0(title.i," (p = ",p.val,")") }
    plot.list[[3]] <-data.frame(
      x = meta.adult.ovaries.by.morph$morph,
      y = cts
    ) %>% 
      mutate(x = stringr::str_replace_all(x, c(
        "LW" = "long-winged\nfemale", "SW" = "short-winged\nfemale"))) %>% 
      mutate(x = as.factor(x)) %>% 
      ggplot(aes(x,y)) +
      theme_bw() +
      geom_violin(color = "gray50") + 
      geom_boxplot(color = "gray50", fill = "gray90", width = 0.2, outlier.shape = NA, coef = 0) +
      geom_jitter(width = 0.2, alpha = 0.75) +
      scale_y_continuous(limits=c(0,NA)) +
      # scale_y_log10(limits=c(1,NA)) +
      labs(x=NULL,y="normalized transcript counts", title = title.i)
  }
  
  # Testes
  title.i <- paste0("testes")
  ids.i <- transcript.id[which(transcript.id %in% rownames(ctn.adult.testes.by.morph))]
  if (length(ids.i)<1) {
    title.i <- paste0("      ",title.i)
    plot.list[[4]] <- data.frame(x = c(1:3)-2, y = c(1:3)-2) %>% 
      ggplot() + theme_void() +
      ggtitle(title.i) + 
      annotate(geom="text", x=0, y=0, label="No expression") 
  } else {
    if (length(ids.i)>1) {
      cts <- colSums(ctn.adult.testes.by.morph[ids.i,])
    } else {
      cts <- ctn.adult.testes.by.morph[ids.i,]
    }
    p.val <- signif(t.test(cts ~ meta.adult.testes.by.morph$morph, var.equal = FALSE)$p.value,2)
    p.val <- p.adjust(p.val, method = "fdr", n=4)
    if (p.val < 0.05) { title.i <- paste0(title.i," (p = ",p.val,")") }
    plot.list[[4]] <-data.frame(
      x = meta.adult.testes.by.morph$morph,
      y = cts
    ) %>% 
      mutate(x = stringr::str_replace_all(x, c(
        "LW" = "long-winged\nmale", "SW" = "short-winged\nmale"))) %>% 
      mutate(x = as.factor(x)) %>% 
      ggplot(aes(x,y)) +
      theme_bw() +
      geom_violin(color = "gray50") + 
      geom_boxplot(color = "gray50", fill = "gray90", width = 0.2, outlier.shape = NA, coef = 0) +
      geom_jitter(width = 0.2, alpha = 0.75) +
      scale_y_continuous(limits=c(0,NA)) +
      # scale_y_log10(limits=c(1,NA)) +
      labs(x=NULL,y="normalized transcript counts", title = title.i)
  }
  
  output.plot <- ggarrange(plotlist = plot.list, ncol = 4, nrow = 1, 
                           labels = panel.labels, widths = c(1,1,0.6,0.6))
  return(output.plot)
} # End of function

# Novel transcription factors
goi.schnurri <- goi.by.morph.plots("schnurri", "TR198894|c1_g1_i2")
goi.mohawk <- goi.by.morph.plots("Mohawk", "TR281024|c0_g1_i1", LETTERS[5:8])
goi.castor <- goi.by.morph.plots("castor", "TR201825|c0_g7_i4", LETTERS[9:12])
goi.sox14 <- goi.by.morph.plots("Sox14", "TR252663|c0_g1_i2", c("","M","N",""))
goi.morph.txf.plot <- ggarrange(goi.schnurri, goi.mohawk, goi.castor, goi.sox14, ncol = 1, nrow = 4)
ggsave("figures/FigS6.GoI.TXF.plots.jpg", goi.morph.txf.plot, width = 13, height = 15, scale = 1)
ggsave("figures/FigS6.GoI.TXF.plots.pdf", goi.morph.txf.plot, width = 13, height = 15, scale = 1)

# Insulin signaling components
goi.inr1 <- goi.by.morph.plots("InR1", "TR202697|c1_g1_i11")
goi.inr2 <- goi.by.morph.plots("InR2", c("TR230867|c1_g4_i10","TR230867|c1_g4_i11","TR230867|c1_g4_i12","TR230867|c1_g4_i14","TR230867|c1_g4_i2","TR230867|c1_g4_i3","TR230867|c1_g4_i4","TR230867|c1_g4_i5","TR230867|c1_g4_i9"), LETTERS[1:4])
goi.akt  <- goi.by.morph.plots("Akt", c("TR275531|c0_g1_i1","TR278472|c6_g6_i2","TR278472|c6_g6_i3"), LETTERS[5:8])
goi.foxo <- goi.by.morph.plots("FoxO", "TR228786|c0_g1_i1", LETTERS[9:12])
ins.plot <- ggarrange(goi.inr2, goi.akt, goi.foxo,
                      ncol = 1, nrow = 3)
ggsave("figures/FigS7.GoI.insulin.plots.png", ins.plot, width = 13, height = 9, scale = 1)
ggsave("figures/FigS7.GoI.insulin.plots.pdf", ins.plot, width = 13, height = 9, scale = 1)

goi.slim <- goi.by.morph.plots("slimfast", c("TR208893|c1_g4_i1"), LETTERS[1:4])
goi.rheb <- goi.by.morph.plots("rheb", c("TR238584|c0_g1_i1"), LETTERS[5:8])
goi.char <- goi.by.morph.plots("charybde", c("TR95086|c0_g1_i1"), LETTERS[9:12])
goi.tor <- goi.by.morph.plots("TOR", c("TR283202|c0_g1_i1"), LETTERS[13:16])
goi.s6k <- goi.by.morph.plots("S6K", c("TR205756|c0_g10_i1","TR205756|c0_g15_i1","TR205756|c0_g2_i1","TR205756|c0_g7_i1"), LETTERS[17:20])
tor.plot <- ggarrange(goi.slim, goi.rheb, goi.char, goi.tor, goi.s6k,
                      ncol = 1, nrow = 5)
ggsave("figures/FigS8.GoI.tor.plots.png", tor.plot, width = 13, height = 15, scale = 1)
ggsave("figures/FigS8.GoI.tor.plots.pdf", tor.plot, width = 13, height = 15, scale = 1)


# Ras: while a t-test indicates differences in the ovaries, DEseq2's model does not.
goi.ras <- goi.by.morph.plots("Ras", c("TR198111|c0_g1_i1"))
ggsave("plots/GoI.Ras.plots.pdf", goi.ras, width = 13, height = 3, scale = 1)


# Others
goi.by.morph.plots("Erk", c("TR239719|c0_g1_i1","TR239719|c0_g1_i2","TR239719|c0_g1_i3"))
goi.by.morph.plots("Mek", c("TR238560|c9_g1_i1","TR238560|c9_g1_i2"))
goi.by.morph.plots("ILPs", c("TR230631|c0_g1_i7","TR261620|c0_g1_i1","TR9168|c0_g1_i1"))
goi.by.morph.plots("Chico", c("TR243146|c2_g1_i1","TR243146|c2_g1_i2","TR243146|c2_g1_i6"))
goi.by.morph.plots("Pi3K", c("TR199615|c0_g1_i1"))
goi.by.morph.plots("PTEN", c("TR201861|c0_g2_i1"))

goi.by.morph.plots("raptor", c("TR203576|c2_g2_i1"))
goi.by.morph.plots("LAMTOR", c("TR110166|c0_g1_i1"))
goi.by.morph.plots("RicTOR", c("TR224123|c0_g3_i3"))
goi.by.morph.plots("Tsc1", c("TR208921|c0_g1_i1"))
goi.by.morph.plots("GSK3", c("TR225601|c1_g3_i2","TR225601|c1_g3_i4"))


goi.by.morph.plots("fln", "TR89791|c0_g1_i1")
goi.by.morph.plots("ErR", "TR215852|c0_g1_i2")
goi.by.morph.plots("sals", "TR252894|c0_g1_i1")
goi.by.morph.plots("ECSIT", "TR135402|c0_g1_i1")
goi.by.morph.plots("FGF", "TR209940|c0_g1_i1")
goi.by.morph.plots("dsxF", "TR199998|c0_g1_i3")
goi.by.morph.plots("dsxM", "TR199998|c0_g1_i2")
goi.by.morph.plots("l(2)tid", "TR282931|c0_g1_i3")
goi.by.morph.plots("Hippo", "TR262253|c0_g1_i2")
goi.by.morph.plots("iap1", "TR211096|c0_g1_i1")
goi.by.morph.plots("sox14", "TR252663|c0_g1_i2")
goi.by.morph.plots("LCMII", "TR237619|c8_g3_i1")
goi.by.morph.plots("H2.0/CG11607", c("TR263459|c0_g1_i4","TR274032|c5_g2_i6"))


