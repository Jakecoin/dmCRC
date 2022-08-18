# BiocManager::install("sva")
suppressMessages(library(RColorBrewer))
suppressMessages(library(sva))
suppressMessages(library(limma))

suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
suppressMessages(library(amap))
suppressMessages(library(ggplot2))
suppressMessages(library(BiocParallel))
suppressMessages(library(YSX))
suppressMessages(library(ggfortify))
suppressMessages(library(patchwork))
suppressMessages(library(ggbeeswarm))
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

meta <- read.table("dm_meta.txt")
mat <- read.table("dm_species.txt")
mat[1:5,1:5]

biological_group = "Group"
batch = "Batch"

meta[[biological_group]] <- factor(meta[[biological_group]])
meta[[batch]] <- factor(meta[[batch]])

# 模型中引入关注的生物变量和其它非批次变量，保留生物差异和非批次差异
modcombat = model.matrix(as.formula(paste('~', biological_group, sep=" ")), data=meta)

# ComBat需要的是matrix
expr_mat_batch_correct <- ComBat(dat=as.matrix(mat), batch=meta[[batch]], mod=modcombat)
expr_mat_batch_correct <- as.data.frame(expr_mat_batch_correct)

expr_mat_batch_correct[1:3,1:4]


