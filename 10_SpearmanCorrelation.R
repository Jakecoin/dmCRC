library(corrplot)
library(dplyr)
require(readxl)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(data.table)
library(tibble)
library(psych)
library(impute)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

target <- c("species", "met", "ko")
group.test <- c("CTRL", "T2DM", "CRC", "T2DM_CRC")
meta <- read.table("dm_meta.txt")
meta$Group <- factor(meta$Group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))

index <- 1

#
otu <- read.table(paste("dm_adj_", target[index], ".txt", sep = ""))
otusig <- read.table("species_log2fc_1_0717.txt")
otusig[grep(otusig[,1], pattern = "CAG|sp\\.|unclassified|bacterium"),] <- NA
otusig <- na.omit(otusig[,1])
table(rownames(otu) %in% otusig)
rownames(otu) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(otu))
otu <- otu[otusig, rownames(meta)]

#met
ev <- read.table("dm_adj_met.txt")
evsig <- read.table("met_sig.txt")
ev <- ev[unique(evsig$kegg),]
metinfo <- read.csv("dm_metinfo.csv")
rownames(ev) <- metinfo$Metabolite[match(rownames(ev), metinfo$KEGG)]

#
ko <- read.table("dm_adj_ko.txt")
kosig <- read.table("ko_q0.05_log2fc_0.5_ko_anno_0716.txt")
ko <- ko[kosig$ko,]

# start
rvalue <- 0.4
pvalue <- 0.05
point.del <- 0

lapply(1:4, function(num){
  gotu <- otu[, meta$Group == group.test[num]]
  gev <- ev[, meta$Group == group.test[num]]
  gko <- ko[, meta$Group == group.test[num]]
  
  # mat <- as.matrix(rbind(gotu, gev))
  cormat <- corr.test(t(gotu), t(gev), method="spearman")
  # write.csv(cormat$r, file=paste("result/9_correlation/OTUvsEV_rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  # write.csv(cormat$p, file=paste("result/9_correlation/OTUvsEV_padjmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r <- -cormat$r
  pdf(paste("result/9_correlation/OTUvsEV_fdr", pvalue, "_", group.test[num],".pdf", sep=""), 
      width = round(ncol(r)/8), 
      height = round(nrow(r)/8))
  corrplot(r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = cormat$p.adj, sig.level = pvalue, insig='label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu')
  )
  dev.off()
  
  # mat <- as.matrix(rbind(gotu, gko))
  cormat <- corr.test(t(gotu), t(gko), method="spearman")
  # write.csv(cormat$r, file=paste("result/9_correlation/OTUvsKO_rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  # write.csv(cormat$p, file=paste("result/9_correlation/OTUvsKO_padjmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r <- -cormat$r
  pdf(paste("result/9_correlation/OTUvsKO_fdr", pvalue, "_", group.test[num],".pdf", sep=""), 
      width = round(ncol(r)/8), 
      height = round(nrow(r)/8))
  corrplot(r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = cormat$p.adj, sig.level = pvalue, insig = 'label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu')
  )
  dev.off()
  
  cormat <- corr.test(t(gev), t(gko), method="spearman")
  # write.csv(cormat$r, file=paste("result/9_correlation/EVvsKO_rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  # write.csv(cormat$p, file=paste("result/9_correlation/EVvsKO_padjmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r <- -cormat$r
  pdf(paste("result/9_correlation/EVvsKO_fdr", pvalue, "_", group.test[num], ".pdf", sep=""), 
      width = round(ncol(r)/8), 
      height = round(nrow(r)/8))
  corrplot(r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = cormat$p.adj, sig.level = pvalue, insig = 'label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu')
  )
  dev.off()
})
