library(Maaslin2)
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

data <- "genus"

mat <- read.table(paste("data/dm_", data,".txt", sep = ""))
meta <- read.table("data/dm_meta.txt")

filter.f <- function(dat, Num){
  Num <- ncol(dat)
  SD <- apply(dat,1,sd)
  num_0 <- apply(dat, 1, function(x) length(which(x == 0)))
  ave_abun <- apply(dat,1,mean)
  tmp <- cbind(dat,data.frame(SD),data.frame(num_0),data.frame(ave_abun))
  colnames(tmp)[(ncol(dat)+1):(ncol(dat)+3)] <- c("sd","count0",'aveabun')
  #dat_filter <- tmp %>% filter(count0 <= as.numeric(Num*0.9) & sd >0) 
  dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0) & (tmp$aveabun > 0.0001),]
  #dat_filter <- tmp[(tmp$count0 <= as.numeric(Num*0.8)) & (tmp$sd >0),]
  return(dat_filter[,1:ncol(dat)])
}

mat <- filter.f(mat, ncol(mat))

##### ko
# mat <- read.table("data/dm_ko.txt")
# table(colnames(mat) == rownames(meta))
# mat <- mat[,rownames(meta)]

mat <- round(mat / colSums(mat), 6)

# mat <- round(mat)
meta$Batch <- ifelse(meta$Batch == 1, "B2021_06", "B2022_01")
meta$Batch <- factor(meta$Batch)
meta$Group <- factor(meta$Group)
levels(meta$Group) <- c("CTRL", "CRC", "T2DM", "T2DM_CRC")

fit_adjust_batch <- adjust_batch(feature_abd = mat,
                                 batch = "Batch",
                                 covariates = c("Gender", "AgeYear", "BMI"),
                                 data = meta)

write.table(fit_adjust_batch$feature_abd_adj, paste("data/dm_adj_", data,".txt", sep = ""), quote = F)

meta$Cancer <- ifelse(is.na(str_extract(meta$Group, "CRC")), "CTRL", "CRC")
fit_meta <- lm_meta(feature_abd = mat,
                      exposure = "Cancer", # 必须是二分类变量
                      batch = "Batch",
                      covariates = c("Gender", "AgeYear", "BMI"),
                      # control = list(rma_method="HS",transform="AST"),
                      data = meta)
  
meta_results <- fit_meta$meta_fits
maaslin_results <- fit_meta$maaslin_fits

write.csv(meta_results, paste('result/CRC_Bacteria_', data,'_ra_unadj_MMUPHin_result.csv', sep = ""))
write.csv(maaslin_results, paste('result/CRC_Bacteria_', data, '_ra_unadj_MMUPHin_maaslin_result.csv', sep = ""))

########### remove batch
library(limma)

x <- removeBatchEffect(mat, batch=meta$Batch, covariates = meta[,2:4])
write.table(x, "data/dm_reba_species.txt", quote = F)

# install_github("lichen-lab/GMPR")

# BiocManager::install("Wrench")
# BiocManager::install("HCBravoLab/Wrench")

library(Wrench)
#extract count and group information for from the mouse microbiome data in the metagenomeSeq package
data(mouseData)
mouseData


##########################

# Ob
group.test <- c("CTRL", "T2DM", "CRC", "T2DM_CRC")
idy <- df_input_metadata$Group %in% group.test[1:2]
table(idy)

fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "output",
  fixed_effects = c("Cancer", "BMI"),
  
  reference = c("Cancer,Health"),
  max_significance = 0.05,
  min_prevalence = 0.1,
  min_abundance = 1)

# Ov
idy <- df_input_metadata$multigroup %in% group.test[3:4]
table(idy)

fit_data = Maaslin2(
  input_data = df_input_data[,idy],
  input_metadata = df_input_metadata,
  output = "Ov_output",
  fixed_effects = c("Cancer", "BMI"),
  reference = c("Cancer,Health"),
  max_significance = 0.05,
  min_prevalence = 0.1,
  min_abundance = 1)

# N
idy <- df_input_metadata$multigroup %in% group.test[5:6]
table(idy)

fit_data = Maaslin2(
  input_data = df_input_data[,idy],
  input_metadata = df_input_metadata,
  output = "N_output",
  fixed_effects = c("Cancer", "BMI"),
  reference = c("Cancer,Health"),
  max_significance = 0.05,
  min_prevalence = 0.1,
  min_abundance = 1)

# Total
fit_data = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "Total_output_log",
  fixed_effects = c("Cancer", "Group", "BMI"),
  reference = c("Group,Normal"),
  max_significance = 0.05,
  min_prevalence = 0.1,
  min_abundance = 1,
  transform = "LOG")

df_input_metadata$CRC_BMI = (df_input_metadata$Cancer == "CRC") *
  df_input_metadata$BMI
df_input_metadata$H_BMI = (df_input_metadata$Cancer == "Health") *
  df_input_metadata$BMI

fit_data_interaction = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "Total_output_interaction",
  fixed_effects = c("Cancer", "Group", "BMI", "CRC_BMI", "H_BMI"),
  reference = c("Group,Normal"),
  max_significance = 0.2,
  min_prevalence = 0.1,
  min_abundance = 1,
  )

