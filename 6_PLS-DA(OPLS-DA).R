# PLS-DA
# https://www.blog4xiang.world/posts/3f4f7ea1.html
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

library(ade4)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)

# install.packages("BiocManager")
# BiocManager::install("ropls")
# install.packages("ropls_1.26.0.tgz", repos = NULL, type="source")
library(ropls)

rundate <- "result/4_PLSDA/"
# 读取输入文件
mat <- read.table("dm_adj_met.txt")
# metinfo <- read.csv("processed_data/metabolite_info_obesity.csv")
mat[1:5,1:5]
dim(mat)

meta <- read.table("dm_meta.txt")

s.mat <- as.data.frame(apply(mat, 1, function(x){scale(x)}))
rownames(s.mat) <- colnames(mat)
s.mat[1:5,1:5]

my_comparison <- combn(as.character(unique(meta$Group)), 2, simplify=FALSE)
my_comparison

mat <- as.data.frame(s.mat)

# 1:6
  i=6
  my_comparison[[i]]
  idx <- meta$Group %in% my_comparison[[i]]
  table(idx)
  met.plsda <- opls(s.mat[idx,], meta[idx,]$Group) # Group不能设定好level , permI = 1000, , predI = 1, orthoI = NA
  
  sample.score <- met.plsda@scoreMN %>% 
                      as.data.frame() %>%
                      mutate(group = meta[idx,]$Group)#, o1 = met.plsda@orthoScoreMN[,1])
    
  # sample.score$p2 <- 0 # i = 3
    p <- ggplot(sample.score, aes(p1, p2, color = group)) +
      geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
      geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
      geom_point() +
      # geom_point(aes(-10,-10), color = 'white') +
      ggtitle(paste(my_comparison[[i]][1], my_comparison[[i]][2], sep = "vs")) +
      labs(x = paste('t1(', round(met.plsda@modelDF$R2X[1]*100), '%)', sep = ''),
           y = paste('t2(', round(met.plsda@modelDF$R2X[2]*100), '%)', sep = '')) + ## 横纵坐标名字可能需要修改 round(met.plsda@modelDF$R2X[1]*100)
      stat_ellipse(level = 0.95, linetype = 'solid', 
                   size = 1, show.legend = FALSE) +
      scale_color_manual(values = c('#84DCC6','#EB8884')) +
      theme_bw() +
      theme(
            legend.position = "right", # none
            legend.title = element_blank(),
            legend.text = element_text(color = 'black',size = 15, face = 'plain'),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_text(color = 'black',size = 15, face = 'plain'),
            axis.title = element_text(color = 'black',size = 15, face = 'plain'),
            axis.ticks = element_line(color = 'black')
      )
    
    pdf(paste(rundate, my_comparison[[i]][1], "vs", my_comparison[[i]][2],".pdf", sep=""), width=6.5, height=5)
    p
    dev.off()
    
    out_met <- met.plsda@vipVn[order(met.plsda@vipVn, decreasing = T)]
    write.table(out_met, paste(rundate, my_comparison[[i]][1], "vs", my_comparison[[i]][2],"_VIP.txt", sep=""),
                quote = FALSE, col.names = FALSE)


# ggplot(sample.score, aes(p1, p2, color = group)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
#   geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
#   geom_point() +
#   # geom_point(aes(-10,-10), color = 'white') +
#   labs(x = 't1(6%)',y = 't2(3%)') + ## 横纵坐标名字可能需要修改
#   stat_ellipse(level = 0.95, linetype = 'solid', 
#                size = 1, show.legend = FALSE) +
#   scale_color_manual(values = c('#84DCC6','#EB8884')) +
#   theme_bw() +
#   theme(
#     #legend.position = c(0.1,0.85),
#     legend.title = element_blank(),
#     legend.text = element_text(color = 'black',size = 15, face = 'plain'),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     axis.text = element_text(color = 'black',size = 15, face = 'plain'),
#     axis.title = element_text(color = 'black',size = 15, face = 'plain'),
#     axis.ticks = element_line(color = 'black')
#   )
# 
# ggsave(paste(rundate, "template.pdf", sep=""))

# BiocManager::install('mixOmics')
library(mixOmics)
