library(stringr)
library(devtools)
library(ggpmisc)
library(Rmisc)
library(ggplot2)
library(ggpubr)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")
manual_color_vector <- c("#16CAB2", "#16697A", "#FF990A", "#FF0022")
col = manual_color_vector[2:4]

outlier_process <- function(data){
  temp <- apply(data, 1, function(x){
    q <- quantile(x, probs=c(.25, .75), na.rm = T)
    iqr <- IQR(x, na.rm = T)
    # caps <- quantile(x, probs=c(.05, .95), na.rm = T)
    # medianx <- median(x)
    x[(x > as.numeric(q[2] + 3 * iqr))] <- q[2] #as.numeric(caps[2]) # 1.5 or 3
    x[(x < as.numeric(q[1] - 3 * iqr))] <- q[1] #as.numeric(caps[1]) # 1.5 or 3
    return(x)
  })
  return(as.data.frame(t(temp)))
}

myggbox <- function(mat) {
  mat <- outlier_process(mat)
  ggplot_input <- as.data.frame(reshape::melt(as.matrix(mat)) %>% dplyr::select(item=1,sample=2,expr=3) %>% 
                                  mutate(group=meta$Group[match(sample, rownames(meta))], expr=as.numeric(expr)))
  ggplot_input$group <- factor(ggplot_input$group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))
  p <- ggplot(ggplot_input, aes(x=group, y=expr, fill=group))+
    geom_boxplot(width=0.8, outlier.size = 0.5) +
    scale_fill_manual(values = manual_color_vector) +
    # stat_compare_means(aes(group = group),
    #                    # method = "wilcox.test",
    #                    label = "p.signif",
    #                    hide.ns = TRUE) + 
    theme(
      legend.position="none",
      panel.background = element_blank(),
      panel.border = element_rect(color = 'black', fill = "transparent"),
      axis.text = element_text(color = 'black',size = 6, face = 'plain'),
      axis.title = element_text(color = 'black',size = 6, face = 'plain'),
    ) +
    labs(x="", y="") + #Relative Abundance
    ggtitle(rownames(mat)) # + #
    # scale_x_discrete(breaks = NULL)
  return(p)
}

myggbar <- function(mat) {
  heatmap_out <- data.frame(
    item = rownames(mat),
    T2DM = apply(mat,1,function(x){log2(mean(x[meta$Group=="T2DM"] + 1)/mean(x[meta$Group=="CTRL"] + 1))}),
    CRC = apply(mat,1,function(x){log2(mean(x[meta$Group=="CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1))}),
    T2DM_CRC = apply(mat,1,function(x){log2(mean(x[meta$Group=="T2DM_CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1))})
  )
  heatmap_out <- heatmap_out[,-1]
  rownames(heatmap_out) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(heatmap_out))
  
  ggbar <- as.data.frame(reshape::melt(as.matrix(heatmap_out)) %>% dplyr::select(name=1, group=2, value=3))
  
  p <- ggplot() + geom_bar(data = ggbar, 
                           aes(x = group, y = value, fill = group), stat = "identity",
                           position = position_dodge(0.8),
                           width=0.6) + 
    scale_fill_manual(values = col) +
    theme(
      legend.position = "none",
      panel.background = element_blank(),
      panel.border = element_rect(color = 'black', fill = "transparent"),
      axis.text = element_text(color = 'black',size = 6, face = 'plain'),
      axis.title = element_text(color = 'black',size = 6, face = 'plain'),
    ) +
    labs(x="", y="") + #log2FC
    ggtitle(rownames(mat)) # + #
    # scale_x_discrete(breaks = NULL)
  return(p)
}

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

meta <- read.table("dm_meta.txt")
target <- c("species", "met", "ko")
group.test <- c("CTRL", "T2DM", "CRC", "T2DM_CRC")

# set the index
index = 1
log2fc.threshold <- 1

mat <- read.table(paste("dm_adj_", target[index], ".txt", sep = ""))
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
  T2DM = apply(mat,1,function(x){log2(mean(x[meta$Group=="T2DM"] + 1)/mean(x[meta$Group=="CTRL"] + 1))}),
  CRC = apply(mat,1,function(x){log2(mean(x[meta$Group=="CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1))}),
  T2DM_CRC = apply(mat,1,function(x){log2(mean(x[meta$Group=="T2DM_CRC"] + 1)/mean(x[meta$Group=="CTRL"] + 1))})
)

# heatmap_out$Abundance = apply(mat,1,mean)

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

# sig_out$Abundance = apply(mat,1,mean)

rownames(sig_out) <- sig_out[,1]
sig_out <- sig_out[,-1]

# order(apply(mat,1,mean), decreasing = T)
# mat <- mat[order(apply(mat,1,mean), decreasing = T),]
# table(apply(mat,1,mean) > 100)

########################################
sig = 3

if(sig == 1) {
  (idx1 <- which(sig_out[,3] & !sig_out[,2] & !sig_out[,1] & abs(heatmap_out[,3]) > log2fc.threshold))
  otu <- mat[idx1,]
} else if(sig == 2) {
  (idx2.1 <- which(sig_out[,2] & sig_out[,3] & !sig_out[,1] & 
                     (heatmap_out[,3] > heatmap_out[,2] + log2fc.threshold & heatmap_out[,2] > 0)
  ))
  (idx2.2 <- which(sig_out[,2] & sig_out[,3] & !sig_out[,1] & 
                     (heatmap_out[,3] < heatmap_out[,2] - log2fc.threshold & heatmap_out[,2] < 0)
  ))
  (idx2 <- c(idx2.1, idx2.2))
  otu <- mat[idx2,]
} else if(sig == 3) {
  (idx3.1 <- which(sig_out[,1] & sig_out[,3] & !sig_out[,2] & 
                     (heatmap_out[,3] > heatmap_out[,1] + log2fc.threshold & heatmap_out[,1] > 0)
  ))
  (idx3.2 <- which(sig_out[,1] & sig_out[,3] & !sig_out[,2] & 
                     (heatmap_out[,3] < heatmap_out[,1] - log2fc.threshold & heatmap_out[,1] < 0)
  ))
  (idx3 <- c(idx3.1, idx3.2))
  otu <- mat[idx3,]
}

dim(otu)
otu <- otu[apply(otu,1,mean) > 50,]
dim(otu)

pdf(paste0(target[index], "_signature", sig, "_box_bar_plot.pdf"), onefile = T, width = 6.5, height = 3.5)
for(i in 1:nrow(otu)) {
  p1 <- myggbox(otu[i,])
  p2 <- myggbar(otu[i,])
  multiplot(p1, p2, cols = 2)
}
dev.off()



