library(psych)
library(reshape2)
library(WGCNA)
library(ggrepel)
library(ggClusterNet)
library(sna)
library(network)

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

temp <- rbind(otu, ev) %>% rbind(ko)

# start
rvalue <- 0.4
pvalue <- 0.0005
point.del <- 1
withname <- FALSE
pdf(paste("result/10_corr_poly/", "gg_r_", rvalue, "_fdr_", pvalue, "_",
          "_network_mean_0802.pdf", sep = ""),
    width = 10, height = 9)
# par(mfrow=c(3, 2)) #, mar=c(1,1,1,1)

group.test <- c("CTRL", "T2DM", "CRC", "T2DM_CRC")

lapply(1:4, function(i){
  idy <- meta$Group %in% group.test[i]
  table(idy)
  mat <- t(temp[,idy])
  mat <- na.omit(mat)
  
  occor <- corr.test(mat, mat,
                     use="pairwise",
                     method="spearman", # 可选pearson/kendall
                     adjust="fdr",
                     alpha=0.05)
  
  write.csv(occor$r, file=paste("result/10_corr_poly/rmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  write.csv(occor$p, file=paste("result/10_corr_poly/pmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  write.csv(occor$p.adj, file=paste("result/10_corr_poly/padjmatrix_", group.test[num],".csv", sep=""), quote=F, row.names=T)
  
  r_matrix <- occor$r
  p_matrix <- occor$p.adj
  
  table(p_matrix < pvalue & abs(r_matrix) > rvalue)
  idx <- p_matrix < pvalue & abs(r_matrix) > rvalue
  r_matrix[!idx]=0
  idx1 <- 1:nrow(otu)
  idx2 <- (nrow(otu)+1):(nrow(otu)+nrow(ev)) #ev
  idx3 <- (nrow(otu)+nrow(ev)+1):nrow(temp) #ko
  r_matrix[idx1, idx1] = 0
  r_matrix[idx2, idx2] = 0
  r_matrix[idx3, idx3] = 0
  
  netClu = data.frame(ID = colnames(mat),
                      group = c(rep("Species", nrow(otu)), 
                                rep("Metabolite", nrow(ev)),
                                rep("KO", nrow(ko))))
  # group1 <- taxinfo$domain[match(netClu$ID, taxinfo$species)]
  # group1 <- ifelse(is.na(group1), netClu$group, group1)
  # netClu$group = as.factor(group1)
  netClu$group = as.factor(netClu$group)
  
  set.seed(12)
  
  result2 = PolygonClusterG(cor = r_matrix, nodeGroup = netClu, zoom = 0.8, zoom2 = 0.8) #PolygonRrClusterG
  node = result2[[1]]
  head(node)
  # ---node节点注释
  # nodes = nodeadd(plotcord = node, otu_table = otu_table, tax_table = tax_table)
  nodes <- node
  nodes$group <- netClu$group[match(node$elements, netClu$ID)]
  if(i == 1 | i == 3){
    foldchange = apply(temp,1,function(x){log2(mean(na.omit(x[meta$Group == "CRC"]))/(mean(na.omit(x[meta$Group == "CTRL"]))))})
  } else if (i == 2) {
    foldchange = apply(temp,1,function(x){log2(mean(na.omit(x[meta$Group == "T2DM"]))/(mean(na.omit(x[meta$Group == "CTRL"]))))})
  } else if (i == 4) {
    foldchange = apply(temp,1,function(x){log2(mean(na.omit(x[meta$Group == "T2DM_CRC"]))/(mean(na.omit(x[meta$Group == "CTRL"]))))})
  }
  # names(foldchange[rownames(nodes)]) == rownames(nodes)
  foldchange <- foldchange[rownames(nodes)]
  nodes$updown <- ifelse(foldchange > 0, "Increased", "Decreased")
  
  # colnames(nodes)
  #-----计算边
  edge = edgeBuild(cor = r_matrix, node = node)
  head(edge)
  # edge$weight <- abseedge
  # colnames(edge)[8] = "cor"
  p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor),
                                data = edge, size = 0.2, alpha = 0.6) +
    geom_point(aes(X1, X2, fill = updown), pch = 21, data = nodes) + #, size = mean
    geom_text_repel(aes(X1, X2, label = elements), size = 2, nudge_y = -0.2, data = nodes,
                    max.overlaps = 20) +
    # geom_text(aes(X1, X2, label = elements), size = 2, nudge_y = -0.2, data = nodes) +
    scale_colour_manual(values = c("+" = "#D00000", "-" = "#009FFD")) +
    scale_fill_manual(values = c("Increased" = "#D00000", "Decreased" = "#009FFD")) +
    scale_size(range = c(4, 14)) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    ggtitle(group.test[i])
  p1
})

dev.off()
.