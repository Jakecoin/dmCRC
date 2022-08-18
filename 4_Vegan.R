library(stringr)
library(devtools)
# install_github("microbiota/amplicon")
# install_github("vqv/ggbiplot")
# library(ggbiplot) # M1/R 版本用不了这个包
library(ade4)

# library(picante)	#用于计算 PD_whole_tree，若不计算它就无需加载。
library(ggplot2)	#用于 ggplot2 作图
library(ggsignif)
library(doBy) 	#用于分组统计
library(ggalt)	#用于绘制拟合曲线
library(ggpubr)
library(dplyr)
library(scales)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(tibble)
library(vegan)	#用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

alpha_index <- function(x,  tree = NULL, base = exp(1)) {
  result <- data.frame(richness = rowSums(x > 0), chao1 = estimateR(x)[3, ], ace = estimateR(x)[5, ],
                       shannon = diversity(x, index = 'shannon'), simpson = diversity(x, index = 'simpson'),
                       pielou = diversity(x, index = 'shannon') / log(estimateR(x)[1, ], base), 
                       goods_coverage = 1 - rowSums(x == 1) / rowSums(x))
  if(!is.null(tree)) {#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- cbind(result, pd[ ,1])
  }
  return(result)
}

# 读取输入文件

mat <- read.table("dm_adj_species.txt") # genus
# mat <- round(mat * 1000000)
meta <- read.table("dm_meta.txt")
# meta <- meta[colnames(mat), ]
table(rownames(meta) == colnames(mat) )

x <- round(t(mat))
alpha.index <- alpha_index(x)

alpha.index$group <- meta$Group

# alpha.index$group <- ifelse(is.na(str_extract(meta$Group, "CRC")), "CTRL", "CRC")
alpha.index$group <- meta$Group
alpha.index$metformin <- meta$metformin

table(rownames(alpha.index) == rownames(meta))

# manual_color_vector = c("#61C9A8", "#FFEEDB", "#A53860", "#7E627F")
manual_color_vector = c("#16CAB2", "#16697A", "#FF990A", "#FF0022")

gginput <- reshape2::melt(alpha.index, id.vars = c("group", "metformin"), measure.vars = c("richness", "chao1", "ace", "shannon", "simpson", "pielou", "goods_coverage"))
gginput$group <- factor(gginput$group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))

my_comparison <- combn(as.character(unique(meta$Group)), 2, simplify=FALSE)
unique(gginput$variable)
# my_comparison <- my_comparison[c(1,2,3,6,8,12,13,14,15)]
plotlist <- list()
wiltestpvalue <- NA

for(i in 1:length(unique(gginput$variable))){
  test <- gginput[gginput$variable==unique(gginput$variable)[i],]
  idj <- NA
  for(j in 1:length(my_comparison)) {
    wiltestpvalue[j] <- wilcox.test(test[test$group==my_comparison[[j]][1],]$value, test[test$group==my_comparison[[j]][2],]$value)$p.value
  }
  
  # wiltestpvalue <- p.adjust(wiltestpvalue, method = "bonferroni")
  idj <- ifelse(wiltestpvalue < 0.05, 1, 0)
  
  plotlist[[i]] <- ggplot(test,aes(x=group, y=value, fill=group))+
    geom_violin(trim=FALSE, aes(linetype=NA)) +
    geom_boxplot(width = 0.25, outlier.size = 0.25) +
    # geom_point(position = position_jitterdodge(),size=0.3)+
    stat_compare_means(comparisons =my_comparison[idj==1],
                       test = "wilcox.test",
                       label = "p.signif",
                       hide.ns = TRUE)+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))+
    labs(title = unique(gginput$variable)[i]) +
    scale_fill_manual(values = manual_color_vector)
}

length(plotlist)
filename <- "alpha_index_boxplot_species.pdf" # species genus
pdf(filename, width = 8, height = 8)
print(ggarrange(plotlist = plotlist[1:7], nrow = 2, ncol = 4))
dev.off()

##############
# CTRL & CRC

# mat <- read.table("data/dm_adj_species.txt")
# mat <- round(mat * 1000000)
# meta <- read.table("data/dm_meta.txt")
# meta <- meta[colnames(mat), ]
# table(rownames(meta) == colnames(mat) )
# 
# x <- round(t(mat))
# alpha.index <- alpha_index(x)
# alpha.index$group <- ifelse(is.na(str_extract(meta$Group, "CRC")), "CTRL", "CRC")
# 
# manual_color_vector = c("#16CAB2", "#FF990A")
# 
# gginput <- reshape2::melt(alpha.index, id.vars = c("group"), measure.vars = c("richness", "chao1", "ace", "shannon", "simpson", "pielou", "goods_coverage"))
# gginput$group <- factor(gginput$group, levels = c("CTRL", "CRC"))
# 
# my_comparison <- combn(as.character(unique(alpha.index$group)), 2, simplify=FALSE)
# unique(gginput$variable)
# plotlist <- list()
# wiltestpvalue <- NA
# 
# for(i in 1:length(unique(gginput$variable))){
#   test <- gginput[gginput$variable==unique(gginput$variable)[i],]
#   idj <- NA
#   for(j in 1:length(my_comparison)) {
#     wiltestpvalue[j] <- wilcox.test(test[test$group==my_comparison[[j]][1],]$value, test[test$group==my_comparison[[j]][2],]$value)$p.value
#   }
#   
#   # wiltestpvalue <- p.adjust(wiltestpvalue, method = "bonferroni")
#   idj <- ifelse(wiltestpvalue < 0.05, 1, 0)
#   
#   plotlist[[i]] <- ggplot(test,aes(x=group, y=value, fill=group))+
#     geom_violin(trim=FALSE, aes(linetype=NA)) +
#     geom_boxplot(width = 0.25, outlier.size = 0.25) +
#     geom_point(position = position_jitterdodge(),size=0.3)+
#     stat_compare_means(comparisons =my_comparison[idj==1],
#                        test = "wilcox.test",
#                        label = "p.signif",
#                        hide.ns = TRUE)+
#     theme_classic()+
#     theme(legend.position = "none",
#           axis.title = element_blank(),
#           axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))+
#     labs(title = unique(gginput$variable)[i]) +
#     scale_fill_manual(values = manual_color_vector)
# }
# 
# length(plotlist)
# filename <- "alpha_index_boxplot_species_2.pdf"
# pdf(filename, width = 8, height = 8)
# print(ggarrange(plotlist = plotlist[1:7], nrow = 2, ncol = 4))
# dev.off()
# 
# # metformin
# alpha.index$group <- ifelse(is.na(str_extract(meta$Group, "CRC")), "CTRL", "CRC")
# 
# {
#   alpha.index$met <- ifelse(meta$Metformin, "Met", "N_Met")
#   alpha.index <- alpha.index[meta$Batch == 2,]
# }
# # rownames(alpha.index) == rownames(meta)
# # table(rownames(alpha.index) == colnames(mat))
# 
# # manual_color_vector = c("#61C9A8", "#FFEEDB", "#A53860", "#7E627F")
# manual_color_vector = c("#16CAB2", "#16697A", "#FF990A", "#FF0022")
# 
# gginput <- reshape2::melt(alpha.index, id.vars = c("group", "met"), measure.vars = c("richness", "chao1", "ace", "shannon", "simpson", "pielou", "goods_coverage"))
# # gginput$group <- factor(gginput$group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))
# 
# my_comparison <- combn(as.character(unique(meta$Group)), 2, simplify=FALSE)
# unique(gginput$variable)
# # my_comparison <- my_comparison[c(1,2,3,6,8,12,13,14,15)]
# plotlist <- list()
# wiltestpvalue <- NA
# 
# for(i in 1:length(unique(gginput$variable))){
#   test <- gginput[gginput$variable==unique(gginput$variable)[i],]
#   idj <- NA
#   for(j in 1:length(my_comparison)) {
#     wiltestpvalue[j] <- wilcox.test(test[test$group==my_comparison[[j]][1],]$value, test[test$group==my_comparison[[j]][2],]$value)$p.value
#   }
#   
#   # wiltestpvalue <- p.adjust(wiltestpvalue, method = "bonferroni")
#   idj <- ifelse(wiltestpvalue < 0.05, 1, 0)
#   
#   plotlist[[i]] <- ggplot(test,aes(x=group, y=value, fill=met))+
#     geom_violin(trim=FALSE, aes(linetype=NA)) +
#     geom_boxplot(width = 0.25, outlier.size = 0.25) +
#     # geom_point(position = position_jitterdodge(),size=0.3)+
#     stat_compare_means(comparisons =my_comparison[idj==1],
#                        test = "wilcox.test",
#                        label = "p.signif",
#                        hide.ns = TRUE)+
#     theme_classic()+
#     theme(legend.position = "none",
#           axis.title = element_blank(),
#           axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))+
#     labs(title = unique(gginput$variable)[i]) +
#     scale_fill_manual(values = manual_color_vector)
# }
# 
# length(plotlist)
# filename <- "alpha_index_boxplot_genus_2.pdf"
# pdf(filename, width = 8, height = 8)
# print(ggarrange(plotlist = plotlist[1:7], nrow = 2, ncol = 4))
# dev.off()

