rm(list = ls())
library(tidyverse)
library(dplyr)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

dm <- read.csv("result/7_vocano_fdr0.05_vip1_log2fc1/vocano_met_T2DM.csv", header=T)
crc <- read.csv("result/7_vocano_fdr0.05_vip1_log2fc1/vocano_met_CRC.csv", header=T)
dmcrc <- read.csv("result/7_vocano_fdr0.05_vip1_log2fc1/vocano_met_T2DM_CRC.csv", header=T)

df <- rbind(dm, crc) %>% rbind(dmcrc)

head(df,1)
df$cluster = ifelse(df$group == "T2DM", 1, ifelse(df$group == "CRC", 2, 3))

#相关包的载入
library(ggplot2)
library(ggrepel)
# 先画背景柱，根据数据log2FC的max值,min值来确定
#根据数据中log2FC区间确定背景柱长度：
col1<-data.frame(x=c(1,2,3),
                 y=c(max(dm$log2fc), max(crc$log2fc), max(dmcrc$log2fc)))
col2<-data.frame(x=c(1,2,3),
                 y=c(min(dm$log2fc), min(crc$log2fc), min(dmcrc$log2fc)))
# 绘制背景柱
p1 <- ggplot()+
  geom_col(data = col1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = col2,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

#把散点火山图叠加到背景柱上：
p2 <- ggplot()+
  # geom_col(data = col1,
  #          mapping = aes(x = x,y = y),
  #          fill = "#dcdcdc",alpha = 0.6)+
  # geom_col(data = col2,
  #          mapping = aes(x = x,y = y),
  #          fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = df,
              aes(x = cluster , y = log2fc, color = TPplotcolor),
              size = 1,
              width =0.4)+
  scale_color_manual(labels = c("Nosig0"="No Significance","UP1"="Significantly Up","UP0"="Up","DOWN1"="Significantly Down","DOWN0"="Down"),
                     values = c("UP0"='#FAC0AE',"DOWN0"='#9BCFF0',"UP1"='#FA2311',"DOWN1"='#6175DB',"Nosig0" ='gray'))+
  labs(x="",y="log2(FoldChange)")

p2

# 添加X轴的分组色块标签：
dfcol<-data.frame(x=c(1:3),
                  y=0,
                  label=c(1:3))
# 添加分组色块标签
dfcol$group <- c("T2DM", "CRC", "T2DM_CRC")
# 加载包
library(RColorBrewer)
library(MetBrewer)

# 自定义分组色块的颜色
tile_color <- met.brewer("Thomas",3)
# 在图中镶嵌色块
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=1,
                     color = "black",
                     fill = tile_color,
                     alpha = 0.6,
                     show.legend = F)+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=group),
            size =3.5,
            color ="white")
p3

p4 <- p3 + geom_text_repel(
  data = df,
  aes(x = cluster, y = log2fc, label = TPplotlabel),
  force = 1.2,
  arrow = arrow(length = unit(0.008, "npc"),
                type = "open", ends = "last"))
p4

# 去除背景，美化图片
p5 <- p4+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 10)
  )
p5



