library(lattice)
library(Formula)
library(readxl)
library(ggpubr)
library(grid)
library(vcd)
library(tibble)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(devtools)
library(biomaRt)
library(ComplexHeatmap)

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

mat <- read.table(paste("dm_adj_", target[index], ".txt", sep = ""))
mat <- mat[, rownames(meta)]

crc.feat <- read.table(paste("result/", target[index], "_maaslin2_CRC.txt", sep = ""),
                       header = T, row.names = 1)
t2dm.feat <- read.table(paste("result/", target[index], "_maaslin2_T2DM.txt", sep = ""),
                        header = T, row.names = 1)
t2dmcrc.feat <- read.table(paste("result/", target[index], "_maaslin2_T2DM_CRC.txt", sep = ""),
                           header = T, row.names = 1)

p.threshold = 0.05

# table(diff[diff$qval < p.threshold,]$value)
table(crc.feat$qval < p.threshold)
table(t2dm.feat$qval < p.threshold)
table(t2dmcrc.feat$qval < p.threshold)

crc.feat <- convert_feature(crc.feat[crc.feat$qval < p.threshold,]$feature)
t2dm.feat <- convert_feature(t2dm.feat[t2dm.feat$qval < p.threshold,]$feature)
t2dmcrc.feat <- convert_feature(t2dmcrc.feat[t2dmcrc.feat$qval < p.threshold,]$feature)

heatmap_out <- data.frame(
  item = rownames(mat),
  T2DM = apply(mat,1,function(x){log2(median(x[meta$Group=="T2DM"] + 1)/median(x[meta$Group=="CTRL"] + 1))}),
  CRC = apply(mat,1,function(x){log2(median(x[meta$Group=="CRC"] + 1)/median(x[meta$Group=="CTRL"] + 1))}),
  T2DM_CRC = apply(mat,1,function(x){log2(median(x[meta$Group=="T2DM_CRC"] + 1)/median(x[meta$Group=="CTRL"] + 1))})
)

# heatmap_out$Abundance = apply(mat,1,median)

heatmap_out <- heatmap_out[,-1]
rownames(heatmap_out) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(heatmap_out))

table(rownames(heatmap_out) %in% t2dmcrc.feat)
rownames(heatmap_out)[!(rownames(heatmap_out) %in% t2dmcrc.feat)]
sig_out <- data.frame(
  item = rownames(heatmap_out),
  T2DM = ifelse(rownames(heatmap_out) %in% t2dm.feat, TRUE, FALSE),
  CRC = ifelse(rownames(heatmap_out) %in% crc.feat, TRUE, FALSE),
  T2DM_CRC = ifelse(rownames(heatmap_out) %in% t2dmcrc.feat, TRUE, FALSE)
)

# sig_out$Abundance = apply(mat,1,median)

rownames(sig_out) <- sig_out[,1]
sig_out <- sig_out[,-1]

library(UpSetR)
# input <- c(
#   T2DM = as.numeric(table(sig_out[,1])[2]),
#   CRC = as.numeric(table(sig_out[,2])[2]),
#   T2DM_CRC = as.numeric(table(sig_out[,3])[2]),
#   "T2DM&CRC" = as.numeric(table(sig_out[,1] & sig_out[,2])[2]),
#   "T2DM&T2DM_CRC" = as.numeric(table(sig_out[,1] & sig_out[,3])[2]),
#   "CRC&T2DM_CRC" = as.numeric(table(sig_out[,2] & sig_out[,3])[2]),
#   "T2DM&CRC&T2DM_CRC" = as.numeric(table(sig_out[,1] & sig_out[,2] & sig_out[,3])[2])
# )

# Plot
# pdf(paste("result/", target[index], "_q", p.threshold, "_upset.pdf", sep = ""))
# upset(fromExpression(input),
#       nintersects = 40,
#       nsets = 3,
#       order.by = "freq",
#       decreasing = T,
#       mb.ratio = c(0.6, 0.4),
#       number.angles = 0,
#       text.scale = 1.1,
#       point.size = 2.8,
#       line.size = 1
# )
# dev.off()

library(ggvenn)
v1 = as.numeric(table(sig_out[,1])[2])
v2 = as.numeric(table(sig_out[,2])[2])
v3 = as.numeric(table(sig_out[,3])[2])
v12 = as.numeric(table(sig_out[,1] & sig_out[,2])[2])
v13 = as.numeric(table(sig_out[,1] & sig_out[,3])[2])
v23 = as.numeric(table(sig_out[,2] & sig_out[,3])[2])
v123 = as.numeric(table(sig_out[,1] & sig_out[,2] & sig_out[,3])[2])

H <-list('T2DM'=c(v1 - v12 - v13 + v123, v12 -v123, v13 - v123, v123),
         'CRC'=c(v2 - v12 - v23 + v123, v12 - v123, v23 - v123, v123),
         'T2DM_CRC'=c(v3 - v13 - v23 + v123, v23 - v123, v13 - v123, v123))

manual_color_vector = c("#009FFD", "#FFA400", "#D00000")

pdf(paste("result/", target[index], "_q", p.threshold, "_venn.pdf", sep = ""), width = 4, height = 4)
ggvenn(H, show_elements=TRUE,
       stroke_color="black",
       fill_color = manual_color_vector,
       fill_alpha = 0.6,
       # stroke_size = 1,
       set_name_size = 6,
       text_color = "black",
       text_size = 6,
       stroke_linetype="solid")
dev.off()

############################################################# choose feature
log2fc.threshold <- 1

# sig 1 dmCRC unique character
x1 <- rownames(mat)[rownames(mat) %in% rownames(sig_out[sig_out[,3] & !sig_out[,2] & !sig_out[,1] & abs(heatmap_out[,3]) > log2fc.threshold,])]
# sig 2 dmCRC enhancer
x2 <- rownames(mat)[rownames(mat) %in% rownames(sig_out[(sig_out[,2] & sig_out[,3] & !sig_out[,1] & 
                                                           (heatmap_out[,3] > heatmap_out[,2] + log2fc.threshold & heatmap_out[,2] > 0)) | (
                                                             sig_out[,2] & sig_out[,3] & !sig_out[,1] & 
                                                               (heatmap_out[,3] < heatmap_out[,2] - log2fc.threshold & heatmap_out[,2] < 0)
                                                           ),])]
# sig 3 dmCRC risk factor
(idx3.1 <- which(sig_out[,1] & sig_out[,3] & !sig_out[,2] &
                   (heatmap_out[,3] > heatmap_out[,1] + log2fc.threshold & heatmap_out[,1] > 0)
))
(idx3.2 <- which(sig_out[,1] & sig_out[,3] & !sig_out[,2] &
                   (heatmap_out[,3] < heatmap_out[,1] - log2fc.threshold & heatmap_out[,1] < 0)
))
(idx3 <- c(idx3.1, idx3.2))
x3 <- rownames(mat)[rownames(mat) %in% rownames(sig_out[idx3,])]
filter.spe <- function(x) {
  return(x[grep(x, pattern = "CAG|sp\\.|unclassified|_bacterium", invert = TRUE)])
}
x1 %>% filter.spe()
x2 %>% filter.spe()
x3 %>% filter.spe()

common <- data.frame(species = filter.spe(x1),
                     group = "dmCRC_Character")
# common <- common[-2,]
non_ob <- data.frame(species = filter.spe(x2),
                     group = "dmCRC_Enhancer")
ob <- data.frame(species = filter.spe(x3),
                 group = "dm_CRC_risk_factor")
# ob <- ob[c(1:5,7,10:18),]
x4 <- rbind(common, non_ob) %>% rbind(ob)
x4
write.table(x4, glue("met_selected_fdr0.05_log2fc_{log2fc.threshold}.txt"))
