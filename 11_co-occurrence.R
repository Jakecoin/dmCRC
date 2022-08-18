# BiocManager::install("preprocessCore")
# 加载包
# library(WGCNA)
library(igraph)
library(psych)
library(impute)
library(corrplot)

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

# start
otusig <- read.table("species_log2fc_1_0717.txt")
# otusig[grep(otusig[,1], pattern = "CAG|sp\\.|unclassified|bacterium"),] <- NA
otusig <- na.omit(otusig[,1])

# evsig <- read.table("met_sig_p0.05_0514.txt")
# kosig <- read.table("ko_sig_p0.05_full.txt")
table(rownames(mat) %in% otusig)
rownames(mat) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(mat))
mat <- mat[otusig, rownames(meta)]
mat <- t(mat)

group.test <- c("CTRL", "T2DM", "CRC", "T2DM_CRC")

# start
rvalue <- 0.4
pvalue <- 0.0005
point.del <- 0
withname <- FALSE
pdf(paste("result/8_co_occurence/", "r_", rvalue, "p_", pvalue, "del_", point.del, "_", withname,
          "_co-occurrence0729_big.pdf", sep = ""),
    width = 20, height = 20)
par(mfrow=c(2, 2)) #, mar=c(1,1,1,1)

lapply(1:4, function(i){
  OTU <- mat
  idy <- meta$Group %in% group.test[i]
  
  table(idy)
  OTU <- OTU[idy,]
  OTU <- na.omit(OTU)
  # 数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue, 但p值需要借助其他函数矫正
  occor = corr.test(OTU, use="pairwise", method="spearman", adjust="fdr", alpha = 0.05)
  occor.r = occor$r # 取相关性矩阵R值
  occor.p = occor$p # 取相关性矩阵p值
  # table(occor$p.adj < 0.05)
  write.table(occor.r, paste("result/8_co_occurence/", group.test[i],
          "_cooccurrence_matrix.txt", sep = ""))
  write.table(occor.p, paste("result/8_co_occurence/", group.test[i],
                             "_cooccurrence_matrix_pvalue.txt", sep = ""))
  
  # occor.r <- as.matrix(read.csv(paste("correlation/OTU_spearman_correlation_rmatrix_", group.test[i],".csv", sep=""), row.names = 1))
  # occor.p <- -as.matrix(read.csv(paste("correlation/OTU_spearman_correlation_pmatrix_", group.test[i],".csv", sep=""), row.names = 1))
  
  # occor.r <- occor.r[oture, oture]
  # occor.p <- occor.p[oture, oture]
  
  # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
  table(occor.p < pvalue & abs(occor.r) > rvalue)
  # table(occor$p.adj < pvalue & abs(occor.r) > rvalue)
  occor.r[occor.p > pvalue | abs(occor.r) < rvalue] = 0 
  
  igraph = graph_from_adjacency_matrix(occor.r, mode = "undirected", weighted=TRUE, diag=FALSE)
  # NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
  # 可以按下面命令转换数据
  # occor.r[occor.r!=0] = 1
  # igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)
  
  # # 是否去掉孤立顶点，根据自己实验而定
  # remove isolated nodes，即去掉和所有 otu 均无相关性的 otu 可省略，前期矩阵已处理过
  if (point.del > 0) {
    bad.vs = V(igraph)[degree(igraph) < point.del]
    igraph = delete.vertices(igraph, bad.vs)
  }
  # 将igraph weight属性赋值到igraph.weight
  igraph.weight = E(igraph)$weight
  # 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
  E(igraph)$weight = NA
  # 设定随机种子数，后续出图都从同一随机种子数出发，保证前后出图形状相对应
  set.seed(123)
  # 如果构建网络时，weighted=NULL,此步骤不能统计
  # sum(igraph.weight>0)# number of postive correlation
  # sum(igraph.weight<0)# number of negative correlation
  
  # set edge color，postive correlation 设定为red, negative correlation设定为blue
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#D04539",ifelse(E.color<0, "#3B94BA", "grey"))
  E(igraph)$color = as.character(E.color)
  
  # 可以设定edge的宽 度set edge width，例如将相关系数与edge width关联
  # summary((abs(igraph.weight) * 10 - rvalue * 10 + 1))
  # E(igraph)$width = (abs(igraph.weight) * 10 - rvalue * 10 + 1)
  
  # 添加OTU注释信息，如分类单元和丰度
  # 另外可以设置vertices size, vertices color来表征更多维度的数据
  # 注意otu_pro.txt文件为我随机产生的数据，因此网络图可能不会产生特定的模式或规律。
  
  # otu_pro = read.table("otu_pro.txt",head=T,row.names=1)
  
  # set vertices size 细菌丰度改变圆大小
  igraph.size = apply(OTU[,V(igraph)$name], 2, mean) # 筛选对应OTU属性
  igraph.size1 = log10(igraph.size + 1) + 1 # 原始数据是什么，为什么*100再取e对数
  summary(igraph.size)
  summary(igraph.size1)
  V(igraph)$size = igraph.size1
  
  # set vertices color
  # igraph.col = otusig[V(igraph)$name,]
  # ver.col <- factor(igraph.col$domain, labels = c("deepskyblue","deeppink"))
  # levels(igraph.col$phylum) = c("green","deeppink","deepskyblue","yellow","brown","pink","gray","cyan","peachpuff") # 直接修改levles可以连值全部对应替换
  # V(igraph)$color = as.character(ver.col)
  
  # V(igraph)$color <- "#EDD382"
  # 改变layout,layout有很多，具体查看igraph官方帮助文档。
  # set.seed(123)
  # plot(igraph,main="Co-occurrence network",layout=layout_with_kk,vertex.frame.color=NA,vertex.label=NA,
  #      edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
  
  # set.seed(123)
  # plot(igraph,main="Co-occurrence network",layout=layout.fruchterman.reingold,vertex.frame.color=NA,vertex.label=NA,
  #      edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
  
  # # 模块性 modularity
  # fc = cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  # modularity = modularity(igraph,membership(fc))
  # # 按照模块为节点配色
  # comps = membership(fc)
  # colbar = rainbow(max(comps))
  # V(igraph)$color = colbar[comps]
  # 最后添加删除color和label项可显示标签和点颜色边框
  set.seed(123)
  if(withname) {
    plot(igraph, main=paste(group.test[i], " Co-occurrence network", sep = ""),
       vertex.frame.color="black", # vertex.label = NA,
       edge.lty=1, edge.curved=FALSE, vertex.label.dist = -1)
  } else {
    plot(igraph, main=paste(group.test[i], " Co-occurrence network", sep = ""),
         vertex.frame.color="black", vertex.label = NA,
         edge.lty=1, edge.curved=FALSE, vertex.label.dist = -1)
  }
  # summary(E(igraph)$width)
  # summary(E(igraph)$size)
  # legend(x=-1.5,y=1.5,levels(factor(V(igraph)$location)),pch=21,col="#777777",pt.bg=vcolor)
})

dev.off()

# pdf(paste("cooccur/", group.test[i], "co-occurrence_legend.pdf", sep = ""), width = 6, height = 6)
# op <- par(mfrow = c(3, 1), mar = rep(1, 4))
# plot(1, type = "n", axes = FALSE, ann = FALSE)
# legend("center", legend=c("Positive", "Negative"), lwd = c(3,3),
#        col=c("#E0837B", "#73B6D3"), lty=1:1, cex=2.5) #, inset=c(0, 0.2)
# plot(1, type = "n", axes = FALSE, ann = FALSE)
# legend("center", legend=c("0.6", "0.8", "1"), lwd = c(1,3,5),
#        col=c("black"), lty=1:1, cex=2.5) #, inset=c(0,-0.2)
# plot(1, type = "n", axes = FALSE, ann = FALSE)
# legend("center", pt.cex = c(3,4,5,6),
#        pch = 19, c("3","4","5","6"), cex=2.5) #, inset=c(0,-0.2)
# dev.off()


# link
# pt.cex = c(min(link$value)/2, (min(link$value)+max(link$value))/4, max(link$value)/2),
# pch = 19, c(as.character(min(link$value)), as.character((min(link$value)+max(link$value))/2), 
#             as.character(max(link$value))), cex=2)

# 
# layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
# # Remove layouts that do not apply to our graph.
# layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
# par(mfrow=c(5,3), mar=c(1,1,1,1)) 
# for (layout in layouts) {
#   print(layout)
#   l <- do.call(layout, list(net_pc))
#   plot(net_pc, vertex.label="",edge.arrow.mode=0, layout=l, main=layout) }

#####################################################33333333333

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]
# start
mat <- read.table("species.txt")
otusig <- read.table("spe_sig.txt")
# evsig <- read.table("met_sig_p0.05_0514.txt")
# kosig <- read.table("ko_sig_p0.05_full.txt")
table(rownames(mat) %in% rownames(otusig))
rownames(mat) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(mat))
table(rownames(mat) %in% rownames(otusig))
mat <- mat[rownames(otusig), rownames(meta)]
mat <- t(mat)

taxinfo <- read.table("species_sig_taxonomy.txt")

occor = corr.test(mat, use="pairwise", method="spearman", adjust="fdr", alpha = 0.05)
occor.r = occor$r # 取相关性矩阵R值
occor.p = occor$p # 取相关性矩阵p值

write.table(occor.r, "sig_cooccurrence_matrix.txt")
write.table(occor.p, "sig_cooccurrence_matrix_pvalue.txt")

occor.r <- as.matrix(read.table("sig_cooccurrence_matrix.txt"))
occor.p <- as.matrix(read.table("sig_cooccurrence_matrix_pvalue.txt"))

# occor.r <- occor.r[spe.selected, spe.selected]
# occor.p <- occor.p[spe.selected, spe.selected]

(order.AOE = corrMatOrder(occor.r, order = 'AOE'))
occor.r = occor.r[order.AOE, order.AOE]
occor.p = occor.p[order.AOE, order.AOE]

pdf("species_correlation_whole.pdf", width = 30, height = 30)
corrplot(-occor.r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
         p.mat = occor.p, sig.level = 0.05, insig='label_sig', pch.cex = 0.9,
         col = COL2('RdYlBu'),
         # order = "hclust", hclust.method = "centroid"
         )
dev.off()

#####################################################33333333333

# group.test <- c("COb", "HOb", "COv", "HOv", "CN", "HN")
# idy <- meta$multigroup %in% group.test[1]

# cormat <- corr.test(mat[idy,], method="spearman", use="complete.obs")
pdf("species_correlation_274species_p0.05.pdf", width = 60, height = 90)
par(mfrow=c(3, 2)) #, mar=c(1,1,1,1)

lapply(1:6, function(i){
  group.test <- c("COb", "HOb", "COv", "HOv", "CN", "HN")
  idy <- meta$multigroup %in% group.test[i]
  # idy <- meta$Group %in% group.test[i]
  
  cormat <- corr.test(mat[idy,], method="spearman", use="complete.obs")
  r <- cormat$r
  p <- cormat$p
  
  r <- r[order.AOE, order.AOE]
  p <- p[order.AOE, order.AOE]
  
  corrplot(-r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = p, sig.level = 0.05, insig='label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu'),
           title = group.test[i],
  )
})

dev.off()

#######################################################################
bacte <- mat[,rownames(taxinfo)[taxinfo$domain == "Bacteria"]]
fungi <- mat[,rownames(taxinfo)[taxinfo$domain == "Eukaryota"]]

finalspe <- read.table("final_selected_species.txt") # species_heatmap_sig_0518.txt
finalspe <- finalspe[,1]
# 
# bacte <- bacte[,colnames(bacte) %in% finalspe]
table(colnames(fungi) %in% finalspe)
fungi <- fungi[,colnames(fungi) %in% finalspe]

# cormat <- corr.test(fungi,bacte, method="spearman", use="complete.obs")
# r <- cormat$r
# p <- cormat$p
# 
# coorder <- order(apply(r, 2, mean), decreasing = T)
# r <- r[,coorder]
# p <- p[,coorder]
# 
# pdf("bac_vs_fungi_correlation.pdf", width = 30, height = 10)
# corrplot(-r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
#                p.mat = p, sig.level = 0.05, insig='label_sig', pch.cex = 0.9,
#                col = COL2('RdYlBu')
#          )
# dev.off()

group.test <- c("COb", "HOb", "COv", "HOv", "CN", "HN")
idy <- meta$multigroup %in% group.test[1]

cormat <- corr.test(fungi[idy,], bacte[idy,], method="spearman", use="complete.obs")
r <- cormat$r
idx <- abs(apply(r, 2, mean)) > 0.4

bacte <- bacte[,idx]
r <- r[,idx]
co_order <- order(apply(r, 2, mean), decreasing = T)

pdf("bac_vs_fungi5_correlation_p0.05.pdf", width = 10, height = 10)
par(mfrow=c(3, 2)) #, mar=c(1,1,1,1)

lapply(1:6, function(i){
  group.test <- c("COb", "HOb", "COv", "HOv", "CN", "HN")
  idy <- meta$multigroup %in% group.test[i]
  # idy <- meta$Group %in% group.test[i]
  
  cormat <- corr.test(fungi[idy,], bacte[idy,], method="spearman", use="complete.obs")
  r <- cormat$r
  p <- cormat$p
  
  r <- r[,co_order]
  p <- p[,co_order]
  
  corrplot(-r, method="color", tl.cex = 0.5, tl.col = "black",# type = "lower",
           p.mat = p, sig.level = 0.05, insig='label_sig', pch.cex = 0.9,
           col = COL2('RdYlBu'),
           title = group.test[i],
  )
})

dev.off()
