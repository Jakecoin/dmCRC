library(stringr)
library(devtools)
# install_github("microbiota/amplicon")
# install_github("vqv/ggbiplot")
library(ggbiplot)
library(vegan)	#用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样
library(picante)	#用于计算 PD_whole_tree，若不计算它就无需加载。
library(ggplot2)	#用于 ggplot2 作图
library(ggsignif)
library(doBy) 	#用于分组统计
library(ggalt)	#用于绘制拟合曲线
library(ggpubr)
library(dplyr)
library(scales)
library(ggsci)

library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(convertid)
library(ggrepel)

library(tibble)
library(dplyr)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

meta <- read.table("data/dm_meta.txt")

# mat <- read.table("data/dm_adj_species.txt")
# ### 从这里开始
# diff <- read.csv("result/1_MMUPHin/MMUPHin_species/CRC_Bacteria_Species_ra_unadj_MMUPHin_maaslin_result.csv")
# colnames(diff)
# s.diff <- diff[diff$Study2021_06.qval < 0.05 & diff$Study2022_01.qval > 0.05,]
# s.diff <- na.omit(s.diff)
# dim(s.diff)
# 
# mat <- mat[s.diff$Study2021_06.feature,]
# dim(mat)

process.mat <- function(target, testgroup) {
  testgroupname <- ifelse(length(testgroup)>1, "TOTAL", ifelse(testgroup == 1, "CTRL", "T2DM"))
  mat <- read.table(paste("data/dm_adj_", target, ".txt", sep=""))
  diff <- read.csv(paste("result/1_MMUPHin/MMUPHin_", target, "/MMUPHin_maaslin_result.csv", sep = ""))
  diff <- na.omit(diff)
  rownames(diff) <- diff[,2]
  diff <- diff[,c(8,18)]
  table((diff[,1] < 0.05 & diff[,2] > 0.05))
  table(diff[,1] > 0.05 & diff[,2] < 0.05)
  if(testgroup == 1){
    idx <- (diff[,1] < 0.05 & diff[,2] > 0.05)
  } else {
    idx <- (diff[,1] > 0.05 & diff[,2] < 0.05)
  }
  s.diff <- diff[idx,]
  s.diff <- rownames(s.diff)
  s.diff
  return(mat[s.diff,])
}

mat <- process.mat("species", 1)
# change
ggplot_input <- as.data.frame(reshape::melt(as.matrix(mat)) %>% dplyr::select(item=1,sample=2,expr=3) %>% 
  mutate(group=meta$Group[match(sample, rownames(meta))], expr=as.numeric(expr),
         status=meta$Group[match(sample, rownames(meta))]))

library(ggplot2)
manual_color_vector <- c("#16CAB2", "#16697A", "#FF990A", "#FF0022")

plotlist=list()
my_comparison <- combn(as.character(unique(meta$Group)), 2, simplify=FALSE)
ggplot_input$group <- factor(ggplot_input$group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))
wiltestpvalue <- NA

for(i in 1:nrow(mat)){
  test <- ggplot_input[ggplot_input$item == rownames(mat)[i],]
  idj <- NA
  wiltestpvalue <- NA
  
  for(j in 1:length(my_comparison)) {
    wiltestpvalue[j] <- wilcox.test(test[test$group==my_comparison[[j]][1],]$expr, test[test$group==my_comparison[[j]][2],]$expr)$p.value
  }
  
  wiltestpvalue <- p.adjust(wiltestpvalue, method = "bonferroni")
  idj <- ifelse(wiltestpvalue < 0.05, 1, 0)
  
  test$group <- factor(test$group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))
  
  plotlist[[i]] <- ggplot(test, aes(x=group, y=expr, fill=group))+
    geom_boxplot(width = 0.4, outlier.size = 0.3) +
    # geom_point(position = position_jitterdodge(), size=0.3)+
    stat_compare_means(comparisons =my_comparison[idj==1],
                       test = "wilcox.test",
                       label = "p.signif",
                       hide.ns = TRUE)+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))+
    labs(title = rownames(mat)[i]) +
    scale_fill_manual(values = manual_color_vector)
}

length(plotlist)

pdf("dm_boxplot_spe_ctrl.pdf", width = 12, height = 15, onefile = TRUE)

for(i in 1:7){
  x <- 20*i-19
  y <- 20*i
  print(ggarrange(plotlist = plotlist[x:y],nrow = 4,ncol = 5))
}

# print(ggarrange(plotlist = plotlist[361:length(plotlist)],nrow = 6,ncol = 3))

dev.off()

