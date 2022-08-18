library(lattice)
library(Formula)
library(readxl)
library(ggpubr)
library(grid)
library(vcd)
library(tibble)
library(RColorBrewer)

library(caret)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(devtools)

library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(convertid)

library(ComplexHeatmap)
library(circlize)
# library(pheatmap)

getwd()
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

mat <- read.table("data/dm_ko.txt") # dm_adj_species.txt
meta <- read.table("data/dm_meta.txt")

# process.mat <- function(target, testgroup) {
#   testgroupname <- ifelse(length(testgroup)>1, "TOTAL", ifelse(testgroup == 1, "CTRL", "T2DM"))
#   mat <- read.table(paste("data/dm_adj_", target, ".txt", sep=""))
#   diff <- read.csv(paste("result/1_MMUPHin/MMUPHin_", target, "/MMUPHin_maaslin_result.csv", sep = ""))
#   diff <- na.omit(diff)
#   rownames(diff) <- diff[,2]
#   diff <- diff[,c(8,18)]
#   table((diff[,1] < 0.05 & diff[,2] > 0.05))
#   table(diff[,1] > 0.05 & diff[,2] < 0.05)
#   if(testgroup == 1){
#     idx <- (diff[,1] < 0.05 & diff[,2] > 0.05)
#   } else {
#     idx <- (diff[,1] > 0.05 & diff[,2] < 0.05)
#   }
#   s.diff <- diff[idx,]
#   s.diff <- rownames(s.diff)
#   s.diff
#   return(mat[s.diff,])
# }

# mat <- process.mat("ko", 2)
meta$Group <- factor(meta$Group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))
diff <- read.table("result/ko_maaslin2/3group_maaslin2.txt")
diff <- diff[diff$qval<0.05,]
View(diff)
s.mat <- mat[unique(diff$feature),]
# s.mat <- 1000000 * mat
# s.mat <- apply(s.mat, 1, scale)
# s.mat <- t(s.mat)
# colnames(s.mat) <- colnames(mat)

# heatmap_out <- data.frame(item = rownames(s.mat),
#                           CTRL = apply(s.mat[, meta$Group == "CTRL"], 1, function(x){median(x)}),
#                           T2DM = apply(s.mat[, meta$Group == "T2DM"], 1, function(x){median(x)}),
#                           CRC = apply(s.mat[, meta$Group == "CRC"], 1, function(x){median(x)}),
#                           T2DM_CRC = apply(s.mat[, meta$Group == "T2DM_CRC"], 1, function(x){median(x)})
#                           )
heatmap_out <- data.frame(
  item = rownames(s.mat),
  T2DM = apply(s.mat,1,function(x){log2(mean(x[meta$Group=="T2DM"] + 1)/mean(x[meta$Group=="CTRL"] + 1) )}),
  CRC = apply(s.mat,1,function(x){log2(mean(x[meta$Group=="CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1) )}),
  T2DM_CRC = apply(s.mat,1,function(x){log2(mean(x[meta$Group=="T2DM_CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1) )})
)
heatmap_out[1:5,]
heatmap_out <- heatmap_out[,-1]
# heatmap_out <- heatmap_out[abs(heatmap_out[,1]) != Inf & abs(heatmap_out[,2]) != Inf & abs(heatmap_out[,3]) != Inf,]

rownames(heatmap_out) <- gsub(heatmap_out$item, pattern = "s__", replacement = "")
manual_color_vector <- c("#16CAB2", "#16697A", "#FF990A", "#FF0022")
col <- colorRamp2(c(-2, 0, 2), c("#16697A", "white", "#FF0022"))
# simple heatmap
# pdf("heatmap.pdf", width = 5, height = 6)
Heatmap(heatmap_out, name = "log2fc", col = col,
          # show_column_names = TRUE, show_row_names = TRUE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = TRUE,
          row_names_gp = gpar(fontsize = 10),
          column_names_rot = 45, column_names_centered = FALSE,
          column_order = c("T2DM", "CRC", "T2DM_CRC")
)
# dev.off()

 # 以下还没用
###############################################################################################
# signdiff <- data.frame(item = rownames(mat),
#                        T2DM = rep(NA,dim(heatmap_out)[1]),
#                        CRC = rep(NA,dim(heatmap_out)[1]),
#                        T2DM_CRC = rep(NA,dim(heatmap_out)[1]))

sig_out <- data.frame(
  item = rownames(heatmap_out),
  T2DM = ifelse(rownames(heatmap_out) %in% diff[diff$value=="T2DM",]$feature, TRUE, FALSE),
  CRC = ifelse(rownames(heatmap_out) %in% diff[diff$value=="CRC",]$feature, TRUE, FALSE),
  T2DM_CRC = ifelse(rownames(heatmap_out) %in% diff[diff$value=="T2DM_CRC",]$feature, TRUE, FALSE)
)

rownames(sig_out) <- sig_out[,1]
sig_out <- sig_out[,-1]
head(sig_out)

heatmap <- scale(heatmap_out)

pdf(paste("result/5_heatmap/", target[index], "_heatmap.pdf", sep = ""), height = 12, width = 2)
print(
  Heatmap(heatmap, name = "log2FC",
          # left_annotation = ha,
          show_column_names = TRUE,# show_row_names = FALSE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 0, column_names_centered = TRUE,
          # col = colorRamp2(c(-max, -7, 0, 7, max), c("#063737", "#109393", "white", "#A8201A", "#A01F18")),
          col = colorRamp2(c(-2, 0, 2), c("#0F8B8D", "white", "#A8201A")),
          column_order  = c("T2DM", "CRC", "T2DM_CRC"),
          row_order = rownames(heatmap), #rownames(heatmap_out[unique(idx),]),
          # width = unit(10, "cm"), height = unit(40, "cm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(sig_out[i,j]) {
              sig <- ifelse(heatmap[i,j] > 0, "+", "-")
              grid.text(sig, x, y, gp = gpar(fontsize = 5))
            }
          }
  )
)
dev.off()

write.table(rownames(heatmap), "ko_sig.txt")
#
kotext <- read.table("ko_database.txt")
colnames(kotext) <- kotext[1,]
kotext <- kotext[-1,]
kotext <- kotext[,-1]

koname <- read.table("ko_sig.txt", sep = "") #ko_sig_p0.05.txt
koname[,1]
koname <- kotext[kotext$ko %in% koname[,1],]
table(koname$type)
ko.duplicated <- duplicated(koname$ko)
koname <- koname[!ko.duplicated,]
write.table(koname, "ko_sig_anno.txt")

heatmap <- heatmap[koname$ko,]
# rownames(heatmap) <- paste(koname[,1], koname[,2], sep = " ")
koname$type <- factor(koname$type, levels = unique(koname$type))
table(koname$type)
colist <- c("#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6", "#086788", "#E9D2F4", "#5B5941", "#340068", "#D88C9A", "#000000",
            "#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6", "#086788")
length(colist)
length(unique(koname$type))

names(colist) <- unique(koname$type)

ha <- rowAnnotation(pathway = koname$type, 
                    col = list(pathway = colist)
)

split = data.frame(type = koname$type)
row_title_rot = 0

pdf(paste("result/5_heatmap/", target[index], "_heatmap_anno.pdf", sep = ""), height = 14, width = 7)
print(
  Heatmap(heatmap, name = "log2FC",
          left_annotation = ha,
          row_split = split, 
          row_title = unique(koname$type),
          row_title_rot = 0,
          row_title_gp = gpar(fontsize = 8),
          show_column_names = TRUE,# show_row_names = FALSE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 0, column_names_centered = TRUE,
          col = colorRamp2(c(-2, 0, 2), c("#0F8B8D", "white", "#A8201A")),
          column_order  = c("T2DM", "CRC", "T2DM_CRC"),
          row_order = rownames(heatmap),
          # width = unit(10, "cm"), height = unit(40, "cm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(sig_out[i,j]) {
              sig <- ifelse(heatmap[i,j] > 0, "+", "-")
              grid.text(sig, x, y, gp = gpar(fontsize = 5))
            }
          }
  )
)
dev.off()

#

dim(s.mat)
ggplot_input <- as.data.frame(reshape::melt(as.matrix(s.mat)) %>% dplyr::select(item=1,sample=2,expr=3) %>% 
                                mutate(group=meta$Group[match(sample, rownames(meta))], expr=as.numeric(expr)))
ggplot_input$group <- factor(ggplot_input$group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))
# ggplot_input$item <- factor(ggplot_input$item, levels = korder)

library(ggplot2)
col <- c("#8EF9F3", "#593C8F", "#FFD9CE", "#DB5461")
# ggplot_input$status <- factor(ggplot_input$status, levels = c("Health", "CRC"))
p <- ggplot(ggplot_input, aes(x=group, y=expr, fill=group))+
  geom_boxplot(width=0.8, outlier.size = 0.5) +
  scale_fill_manual(values = col)+
  # geom_point(position = position_jitterdodge(),size=0.1)+
  # stat_compare_means(aes(group = status),
  #                    label = "p.signif",
  #                    method = "wilcox.test", 
  #                    hide.ns = T)+
  # theme_bw() + 
  theme(
    legend.title = element_blank(),
    legend.text = element_text(color = 'black',size = 6, face = 'plain'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = "transparent"),
    #panel.grid = element_blank(),
    axis.text = element_text(color = 'black',size = 6, face = 'plain'),
    axis.title = element_text(color = 'black',size = 6, face = 'plain'),
    #axis.ticks = element_line(color = 'black')
  )+
  labs(x="", y="Relative Abundance)") +  
  facet_wrap(.~item, ncol = 5, nrow = 29, scales="free_y")

pdf("ko_sig_boxplot.pdf", onefile = T, width = 12, height = 50)
p
dev.off()

#

#pathview
koname <- read.table("ko_sig_anno.txt") #ko_sig_p0.05_0424.txt
colnames(koname) <- koname[1,]
koname <- koname[-1,]
koname <- koname[,-1]

k.obdiff <- koname$ko
ko.pathway.id <- unique(koname$subtypecode)
dim(heatmap_out)
ko.data <- heatmap_out$Obesityq
names(ko.data) <- rownames(heatmap_out)
# ko.data=sim.mol.data(mol.type="gene.ko", nmol=5000)

index = 3
group.index = 1
# obdiff <- read.table(paste("Result/", save_path[index], "/", group.test[group.index], "_maaslin2.txt", sep = ""), 
#                      header = T, row.names = 1)
# num <- 6 # 6 pvalue / 8 fdr
# threshold = 0.05
# table(obdiff[,num] < threshold)
# k.obdiff <- obdiff[obdiff[,num] < threshold,]$feature

# cpd.data
# cpd <- read.table("met.txt")
# cpd.fc <- data.frame(
#   item = rownames(cpd),
#   Normal = apply(cpd,1,function(x){log2(mean(x[meta$multigroup=="CN"])/mean(x[meta$multigroup=="HN"]))}),
#   Overweight = apply(cpd,1,function(x){log2(mean(x[meta$multigroup=="COv"])/mean(x[meta$multigroup=="HOv"]))}),
#   Obesity = apply(cpd,1,function(x){log2(mean(x[meta$multigroup=="COb"])/mean(x[meta$multigroup=="HOb"]))})
# )
# 
# cpd.data <- cpd.fc$Obesity
# names(cpd.data) <- cpd.fc$item

k.obdiff <- read.table("")
lapply(1:length(ko.pathway.id), function(i) { # diff / k.obdiff / m.obdiff / koname$ko
  pv.out <- pathview(gene.data = ko.data[k.obdiff],  pathway.id = ko.pathway.id[i],# cpd.data = cpd.data,
                     species = "ko", out.suffix = "ko.data", kegg.native = T)
})














# column_od <- hclust(dist(t(diff)))$order

col2 <- colorRamp2(c(-floor(max(diff)+1), 0, floor(max(diff)+1)), c("#0F8B8D", "white", "#A8201A"))

circos.par(gap.after = c(10))

split = NA
split[(diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr > 0.05 & diff$ObC_H_fdr < 0.05 & abs(diff$ObC_H_log2fc) > 1)] = 1
split[(diff$NC_H_fdr > 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr < 0.05 & 
         abs(diff$ObC_H_log2fc) > 1 & abs(diff$OvC_H_log2fc) > 1)] = 2
split[(diff$NC_H_fdr < 0.05 & diff$OvC_H_fdr < 0.05 & diff$ObC_H_fdr > 0.05 & 
         abs(diff$NC_H_log2fc) > 2 & abs(diff$OvC_H_log2fc) > 2)] = 3
split = factor(split, levels = 1:3)

circos.heatmap(diff[,1:3], col = col2, split = split,
               cluster = TRUE, # dend.side = "inside",
               rownames.side = "outside", rownames.cex = 0.2,
               track.height = 0.1)

circos.track(
  track.index = get.current.track.index(), 
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 1) {
      # 在最后一个扇形中
      cn = c("N", "Ov", "Ob")
      n = length(cn)
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
        1:n - 0.5, cn, cex = 0.5, adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }, 
  bg.border = NA
)

lgd <- Legend(title = "log2FC", col_fun = col2)
grid.draw(lgd)

circos.clear()
