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
library(circlize)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

target <- c("species", "met", "ko")
group.test <- c("CTRL", "T2DM", "CRC", "T2DM_CRC")
meta <- read.table("dm_meta.txt")
meta$Group <- factor(meta$Group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))

index <- 3

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

crc.feat <- crc.feat[crc.feat$qval < p.threshold,]$feature
t2dm.feat <- t2dm.feat[t2dm.feat$qval < p.threshold,]$feature
t2dmcrc.feat <- t2dmcrc.feat[t2dmcrc.feat$qval < p.threshold,]$feature

heatmap_out <- data.frame(
  item = rownames(mat),
  T2DM = apply(mat,1,function(x){log2(mean(x[meta$Group=="T2DM"] + 1)/mean(x[meta$Group=="CTRL"] + 1))}),
  CRC = apply(mat,1,function(x){log2(mean(x[meta$Group=="CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1))}),
  T2DM_CRC = apply(mat,1,function(x){log2(mean(x[meta$Group=="T2DM_CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1))})
)

heatmap_out <- heatmap_out[,-1]

sig_out <- data.frame(
  item = rownames(heatmap_out),
  T2DM = ifelse(rownames(heatmap_out) %in% t2dm.feat, TRUE, FALSE),
  CRC = ifelse(rownames(heatmap_out) %in% crc.feat, TRUE, FALSE),
  T2DM_CRC = ifelse(rownames(heatmap_out) %in% t2dmcrc.feat, TRUE, FALSE)
)

rownames(sig_out) <- sig_out[,1]
sig_out <- sig_out[,-1]

library(UpSetR)
input <- c(
  T2DM = as.numeric(table(sig_out[,1])[2]),
  CRC = as.numeric(table(sig_out[,2])[2]),
  T2DM_CRC = as.numeric(table(sig_out[,3])[2]),
  "T2DM&CRC" = as.numeric(table(sig_out[,1] & sig_out[,2])[2]),
  "T2DM&T2DM_CRC" = as.numeric(table(sig_out[,1] & sig_out[,3])[2]),
  "CRC&T2DM_CRC" = as.numeric(table(sig_out[,2] & sig_out[,3])[2]),
  "T2DM&CRC&T2DM_CRC" = as.numeric(table(sig_out[,1] & sig_out[,2] & sig_out[,3])[2])
)

# Plot
pdf(paste("result/", target[index], "_q", p.threshold, "_upset.pdf", sep = ""))
upset(fromExpression(input),
      nintersects = 40,
      nsets = 3,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1
)
dev.off()

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

# heatmap
# s.mat <- apply(mat, 1, scale)
# s.mat <- t(s.mat)
# colnames(s.mat) <- colnames(mat)
# meta$MS <- ifelse(meta$Mutation == "Mutant", "MSI-H", ifelse(meta$Mutation == "WildType", "MSS", NA))
# # meta$Location <- ifelse(meta$Location == "/", NA, meta$Location)
# # col_2 <- colorRamp2(c(0, 1), c("#2FBF71", "#ED7D3A"))
# column_ha = HeatmapAnnotation(Status = meta$Cancer, 
#                               BMI = factor(meta$Group, levels = c("Normal", "Overweight", "Obesity")),
#                               Age = factor(meta$Age, levels = c("Young", "Elder")),
#                               Gender = factor(meta$Gender, levels = c("Male", "Female")),
#                               # Location = factor(meta$Location, levels = c("Left", "Right", "Rectal")),
#                               Differentiation = factor(meta$Differentiation, levels = c("Well-Moderate", "Poor")),
#                               Lymphatic_Metastasis = meta$Lymphatic_Invasion,
#                               # Vascular_Invasion = meta$Vascular_Invasion,
#                               Stage = factor(meta$T),
#                               MSI = meta$MS,
#                               col = list(Status = c("CRC" = "#1F1F1F", "Health" = "#F0F3BD"),
#                                          BMI = c("Normal" = "#009FFD", "Overweight" = "#FFA400", "Obesity" = "#D00000"), 
#                                          Age = c("Young" = "#1F1F1F", "Elder" = "#F0F3BD"),
#                                          Gender = c("Male" = "#1F1F1F", "Female" = "#F0F3BD"),
#                                          # Location = c("Left" = "#C64191", "Right" = "#62BFED", "Rectal" = "#297373"),
#                                          Differentiation = c("Well-Moderate" = "#1F1F1F", "Poor" = "#F0F3BD"),
#                                          Lymphatic_Metastasis = c("Positive" = "#1F1F1F", "Negative" = "#F0F3BD"),
#                                          # Vascular_Invasion = c("Positive" = "#1F1F1F", "Negative" = "#F0F3BD"),
#                                          Stage = c("0" = "#F0F3BD", "1" = "#DEBAC0", "2" = "#C64191", "3" = "#62BFED",
#                                                    "4" = "#297373"),
#                                          MSI = c("MSI-H" = "#1F1F1F", "MSS" = "#F0F3BD")
#                               )
# )
# 
# # col = list(Status = c("CRC" = "#028090", "Health" = "#F0F3BD"),
# #            BMI = c("Normal" = "#009FFD", "Overweight" = "#FFA400", "Obesity" = "#D00000"), 
# #            Age = c("Young" = "#028090", "Elder" = "#F0F3BD"),
# #            Gender = c("Male" = "#028090", "Female" = "#F0F3BD"),
# #            # Location = c("Left" = "#C64191", "Right" = "#62BFED", "Rectal" = "#297373"),
# #            Differentiation = c("Well-Moderate" = "#028090", "Poor" = "#F0F3BD"),
# #            Lymphatic_Metastasis = c("Positive" = "#028090", "Negative" = "#F0F3BD"),
# #            # Vascular_Invasion = c("Positive" = "#028090", "Negative" = "#F0F3BD"),
# #            Stage = c("0" = "#F0F3BD", "1" = "#DEBAC0", "2" = "#C64191", "3" = "#62BFED",
# #                      "4" = "#297373"),
# #            MSI = c("MSI-H" = "#028090", "MSS" = "#F0F3BD")
# 
# p <- Heatmap(s.mat[rownames(s.mat) %in% OCTRL.feat,], name = "Z-score", border = TRUE,
#         show_column_names = FALSE, show_row_names = TRUE,
#         top_annotation = column_ha,
#         column_split = factor(meta$Group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC")),
#         cluster_column_slices = FALSE,
#         column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
#         row_names_gp = gpar(fontsize = 10),
#         column_names_rot = 0, column_names_centered = TRUE,
#         col = c("#15099A", "white", "#B7151B") #colorRamp2(c(-1, 0, 1), c("#15099A", "white", "#B7151B")) #
# )
# p

# pdf(paste(target[index], "_whole_heatmap.pdf", sep = ""), width = 20, height = 10) # ,
# print(p)
# dev.off()

# choose feature
log2fc.threshold <- 0.5

# sig 1 dmCRC character
(idx1 <- which(sig_out[,3] & !sig_out[,2] & !sig_out[,1] & abs(heatmap_out[,3]) > log2fc.threshold))

# sig 2 dmCRC enhancer
(idx2.1 <- which(sig_out[,2] & sig_out[,3] & !sig_out[,1] & 
                   (heatmap_out[,3] > heatmap_out[,2] + log2fc.threshold & heatmap_out[,2] > 0)
))
(idx2.2 <- which(sig_out[,2] & sig_out[,3] & !sig_out[,1] & 
                   (heatmap_out[,3] < heatmap_out[,2] - log2fc.threshold & heatmap_out[,2] < 0)
))
(idx2 <- c(idx2.1, idx2.2))

# sig 3 dmCRC character
(idx3.1 <- which(sig_out[,1] & sig_out[,3] & !sig_out[,2] & 
                   (heatmap_out[,3] > heatmap_out[,1] + log2fc.threshold & heatmap_out[,1] > 0)
))
(idx3.2 <- which(sig_out[,1] & sig_out[,3] & !sig_out[,2] & 
                   (heatmap_out[,3] < heatmap_out[,1] - log2fc.threshold & heatmap_out[,1] < 0)
))
(idx3 <- c(idx3.1, idx3.2))

#
sig = 3
#

if(sig == 1) {
  idx = idx1
} else if (sig == 2) {
  idx = idx2
} else {
  idx = idx3
}
# idx <- c(idx1,idx2,idx3)

signif <- sig_out[idx,] # rownames(sig_out) %in% idx
heatmap <- heatmap_out[idx,] #1

dim(signif)
dim(heatmap)

write.table(rownames(heatmap), paste(target[index], "_signature", sig, "_q", p.threshold, "_log2fc_", log2fc.threshold,"_heatmap_sig_0716.txt", sep = ""))
# pdf(paste("result/", target[index], "_q", p.threshold, "_heatmap_sig_mean_0716_1.pdf", sep = ""), height = 60, width = 4)
# print(
#   Heatmap(heatmap, name = "log2FC",
#           # left_annotation = ha,
#           show_column_names = TRUE,# show_row_names = FALSE,
#           column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
#           row_names_gp = gpar(fontsize = 7),
#           column_names_gp = gpar(fontsize = 4),
#           column_names_rot = 0, column_names_centered = TRUE,
#           # col = colorRamp2(c(-max, -7, 0, 7, max), c("#063737", "#109393", "white", "#A8201A", "#A01F18")),
#           col = colorRamp2(c(-2,0,2), c("#0F8B8D", "white", "#A8201A")),
#           column_order  = c("T2DM", "CRC", "T2DM_CRC"), # "CTRL", 
#           row_order = rownames(heatmap), #rownames(heatmap_out[unique(idx),]),
#           # width = unit(10, "cm"), height = unit(40, "cm"),
#           cell_fun = function(j, i, x, y, width, height, fill) {
#             if(signif[i,j]) {
#               sig <- ifelse(heatmap[i,j] > 0, "+", "-")
#               grid.text(sig, x, y, gp = gpar(fontsize = 5))
#             }
#           }
#   )
# )
# 
# dev.off()

######################################################################################################
kotext <- read.table("ko_database.txt")
colnames(kotext) <- kotext[1,]
kotext <- kotext[-1,]
kotext <- kotext[,-1]

koname <- rownames(heatmap)
# koname <- read.table("ko_q5e-04_heatmap_sig_0716_1.txt", sep = "")[,1]
koname <- kotext[kotext$ko %in% koname,]
table(koname$type)
table(koname$subtypecode)
# 需要调整代码，一个 KO 基因匹配多个 pathway
ko.duplicated <- duplicated(koname$ko)
koname <- koname[!ko.duplicated,]
# write.table(koname[!duplicated(koname$ko),], paste(target[index], "_q", p.threshold, "_log2fc_", log2fc.threshold, "_ko_anno_0716.txt", sep = ""))

heatmap <- heatmap[koname$ko,]
# rownames(heatmap) <- paste(koname[,1], koname[,2], sep = " ")
koname$type <- factor(koname$type, levels = unique(koname$type))
table(koname$type)
colvec <- c("#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6")
colist <- rep(colvec, length.out = length(unique(koname$type)))

names(colist) <- unique(koname$type)

ha <- rowAnnotation(pathway = koname$type, 
                    col = list(pathway = colist))

split = data.frame(type = koname$type)
row_title_rot = 0

pdf(paste("result/", target[index], "_signature", sig,  "_q", p.threshold, "_log2fc_", log2fc.threshold,"_heatmap_anno_0718_type.pdf", sep = ""), height = nrow(heatmap)/9 + 1, width = 7)
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
          # col = colorRamp2(c(-max, -7, 0, 7, max), c("#063737", "#109393", "white", "#A8201A", "#A01F18")),
          col = colorRamp2(c(-2, 0, 2), c("#0F8B8D", "white", "#A8201A")),
          column_order  = c("T2DM", "CRC", "T2DM_CRC"),
          row_order = rownames(heatmap),
          # width = unit(10, "cm"), height = unit(40, "cm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(signif[i,j]) {
              sig <- ifelse(heatmap[i,j] > 0, "+", "-")
              grid.text(sig, x, y, gp = gpar(fontsize = 5))
            }
          }
  )
)
dev.off()

##################################################################################
## subtype
koname$subtype[is.na(koname$subtype)] <- "NA"
koname$subtype <- factor(koname$subtype, levels = unique(koname$subtype))
table(koname$subtype)
colvec <- c("#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6")
colist <- rep(colvec, length.out = length(unique(koname$subtype)))

names(colist) <- unique(koname$subtype)

ha <- rowAnnotation(pathway = koname$subtype, 
                    col = list(pathway = colist))

split = data.frame(type = koname$subtype)
row_title_rot = 0

koname$koname <- ifelse(is.na(koname$koname), koname$ko, koname$koname)

pdf(paste("result/", target[index], "_signature", sig, "_q", p.threshold, "_log2fc_", log2fc.threshold, "_heatmap_anno_subtype_0718.pdf", sep = ""), height =  nrow(heatmap)/9+1, width = 10)
print(
  Heatmap(heatmap, name = "log2FC",
          row_labels = koname$koname,
          left_annotation = ha,
          row_split = split, 
          row_title = unique(koname$subtype),
          row_title_rot = 0,
          row_title_gp = gpar(fontsize = 8),
          show_column_names = TRUE,# show_row_names = FALSE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 0, column_names_centered = TRUE,
          # col = colorRamp2(c(-max, -7, 0, 7, max), c("#063737", "#109393", "white", "#A8201A", "#A01F18")),
          col = colorRamp2(c(-2, 0, 2), c("#0F8B8D", "white", "#A8201A")),
          column_order  = c("T2DM", "CRC", "T2DM_CRC"),
          row_order = rownames(heatmap),
          # width = unit(10, "cm"), height = unit(40, "cm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(signif[i,j]) {
              sig <- ifelse(heatmap[i,j] > 0, "+", "-")
              grid.text(sig, x, y, gp = gpar(fontsize = 5))
            }
          }
  )
)
dev.off()
