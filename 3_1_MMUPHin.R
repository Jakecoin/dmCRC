library(Maaslin2)
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
library(stringr)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

target <- c("species", "met", "ko", "genus")
group.test <- c("CTRL", "CRC", "T2DM", "T2DM_CRC")
index = 1

meta <- read.table("dm_meta.txt")
dim(meta)
mat <- read.table(paste("dm_", target[index], ".txt", sep = ""))
dim(mat)
table(colnames(mat) == rownames(meta))

# filter.f <- function(dat, Num){
#   Num <- ncol(dat)
#   SD <- apply(dat,1,sd)
#   num_0 <- apply(dat, 1, function(x) length(which(x == 0)))
#   ave_abun <- apply(dat,1,mean)
#   tmp <- cbind(dat,data.frame(SD),data.frame(num_0),data.frame(ave_abun))
#   colnames(tmp)[(ncol(dat)+1):(ncol(dat)+3)] <- c("sd","count0",'aveabun')
#   #dat_filter <- tmp %>% filter(count0 <= as.numeric(Num*0.9) & sd >0) 
#   dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0) & (tmp$aveabun > 0.0001),]
#   #dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0),]
#   return(dat_filter[,1:ncol(dat)])
# }
# mat <- filter.f(mat, ncol(mat))
mat[1:5,1:5]
mat <- round(mat) #*10
meta$Batch <- factor(meta$Batch)
meta$Group <- factor(meta$Group)
levels(meta$Group) <- c("CTRL", "CRC", "T2DM", "T2DM_CRC")
# meta$Cancer <- "CRC"
# meta$Cancer[meta$Group == "CTRL" | meta$Group == "T2DM"] <- "CTRL"

fit_adjust_batch <- adjust_batch(feature_abd = mat,
                                 batch = "Batch",
                                 covariates = "Group",
                                 data = meta)

dim(fit_adjust_batch$feature_abd_adj)
write.table(fit_adjust_batch$feature_abd_adj, paste("dm_adj_", target[index], ".txt", sep = ""), quote = F)

# meta$Cancer <- ifelse(is.na(str_extract(meta$Group, "CRC")), "CTRL", "CRC")
# 
# i = 2
# idx <- meta$Group == "CTRL" | meta$Group == group.test[i]
# 
# lapply(2:4, function(i){
#   fit_meta <- lm_meta(feature_abd = mat[,idx],
#                       exposure = "Group", # 必须是二分类变量
#                       batch = "Batch",
#                       covariates = c("Sex", "Age", "BMI"),
#                       # control = list(rma_method="HS",transform="AST"),
#                       data = meta[idx,])
# 
#   meta_results <- fit_meta$meta_fits
#   maaslin_results <- fit_meta$maaslin_fits
# 
#   write.csv(meta_results, paste('result/CRC_Bacteria_', data,'_ra_unadj_MMUPHin_result.csv', sep = ""))
#   write.csv(maaslin_results, paste('result/CRC_Bacteria_', data, '_ra_unadj_MMUPHin_maaslin_result.csv', sep = ""))
# })

########### remove batch
# library(limma)
# 
# x <- removeBatchEffect(mat, batch=meta$Batch, covariates = meta[,2:4])
# write.table(x, "data/dm_reba_species.txt", quote = F)
# 
# # install_github("lichen-lab/GMPR")
# # BiocManager::install("Wrench")
# # BiocManager::install("HCBravoLab/Wrench")
# 
# library(Wrench)
# #extract count and group information for from the mouse microbiome data in the metagenomeSeq package
# data(mouseData)
# mouseData
