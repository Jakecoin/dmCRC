library(Maaslin2)
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
library(tidyverse)

outlier_process <- function(data){
  temp <- apply(data, 1, function(x){
    q <- quantile(x)
    iqr <- q[4]-q[2]
    median <- median(x)
    x[x > as.numeric(q[4] + 3 * iqr) | x < as.numeric(q[2] - 3 * iqr)] <- NA # 1.5 or 3
    return(x)
  })
  return(t(temp))
}

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

df_input_metadata <- read.table("dm_meta.txt")
# df_input_metadata <- df_input_metadata[df_input_metadata$discover == 1, ]
# df_input_data <- df_input_data/colSums(df_input_data)

target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("T2DM", "CRC", "T2DM_CRC")

index <- 3

df_input_data <- read.table(paste("dm_adj_", target[index], ".txt", sep = ""))
df_input_data[1:5,1:5]
# df_input_data <- outlier_process(df_input_data)

for(i in 1:3) {
  idy <- df_input_metadata$Group %in% c(group.test[i], "CTRL")
  table(idy)
  
  fit_data <- Maaslin2(
    input_data = df_input_data[,idy],
    input_metadata = df_input_metadata,
    output = paste("result/", save_path[index], "/", group.test[i], "_output", sep = ""),
    fixed_effects = c("Group", "Sex", "Age", "BMI"), #纳入病程，空腹血糖
    # random_effects = "Batch", ##########################
    reference = c("Group, CTRL"),
    max_significance = 0.05,
    min_prevalence = 0.1,
    min_abundance = 1)
  
  maaslin2_all_results <- fit_data$results # Save results table
  maaslin2_results <- maaslin2_all_results %>% filter(metadata == 'Group') %>% arrange(pval) # Discard covariate associations
  maaslin2_results$qval <- p.adjust(maaslin2_results$pval, method = 'BH') # FDR correction using 'BH'
  write.table(maaslin2_results, paste("result/", save_path[index], "_", group.test[i], ".txt", sep = ""),
              quote = FALSE)
}

# # Total
# fit_data = Maaslin2(
#   input_data = df_input_data,
#   input_metadata = df_input_metadata,
#   output = "Total_output_log",
#   fixed_effects = c("Cancer", "Group", "BMI"),
#   reference = c("Group,Normal"),
#   max_significance = 0.05,
#   min_prevalence = 0.1,
#   min_abundance = 1,
#   transform = "LOG")
# 
# df_input_metadata$CRC_BMI = (df_input_metadata$Cancer == "CRC") *
#   df_input_metadata$BMI
# df_input_metadata$H_BMI = (df_input_metadata$Cancer == "Health") *
#   df_input_metadata$BMI
# 
# fit_data_interaction = Maaslin2(
#   input_data = df_input_data,
#   input_metadata = df_input_metadata,
#   output = "Total_output_interaction",
#   fixed_effects = c("Cancer", "Group", "BMI", "CRC_BMI", "H_BMI"),
#   reference = c("Group,Normal"),
#   max_significance = 0.2,
#   min_prevalence = 0.1,
#   min_abundance = 1,
#   )


# meta_analysis
# 样例
# meta.all <- group_input()
# rownames(meta.all) <- meta.all
# # meta.all$StudyID <- factor(meta.all$country)

# feat.abu <- metagenome_input()
# feat.abu <- feat.abu[,rownames(meta.all)]
# 
# feat.abu[is.na(feat.abu)] <- 0
# feat.abu <- round(feat.abu)
# 
# fit_meta <- lm_meta(feature_abd = feat.abu,
#                     exposure = "Cancer",
#                     batch = "Batch",
#                     covariates = c("Gender", "AgeYear"),
#                     control = list(rma_method="HS",transform="AST"),
#                     data = meta.all)
# 
# meta_results <- fit_meta$meta_fits
# maaslin_results <- fit_meta$maaslin_fits
# 
# write.csv(meta_results, 'CRC_Bacteria_Species_ra_unadj_MMUPHin_result.csv')
# write.csv(maaslin_results, 'CRC_Bacteria_Species_ra_unadj_MMUPHin_maaslin_result.csv')

#########################################################################################

# input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
# input_data
# 
# input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file
# input_metadata
# 
# #get the pathway (functional) data - place holder
# download.file("https://raw.githubusercontent.com/biobakery/biobakery_workflows/master/examples/tutorial/stats_vis/input/pathabundance_relab.tsv", "./pathabundance_relab.tsv")
# 
# df_input_data = read.table(file = input_data, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
# df_input_data[1:5, 1:5]
# 
# df_input_metadata = read.table(file = input_metadata, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
# df_input_metadata[1:5, ]
# 
# df_input_path = data.frame(read.csv("./pathabundance_relab.tsv", sep = "\t", stringsAsFactors = F, row.names = 1))
# df_input_path[1:5, 1:5]
# 
# #This can also be done with with the HUMAnN 3 untiliy `humann_split_stratified_table`
# unstrat_pathways <-function(dat_path){
#   temp = dat_path[!grepl("\\|",rownames(dat_path)),]
#   return(temp)
# }
# 
# df_input_path = unstrat_pathways(df_input_path)
# 
# fit_func = Maaslin2(
#   input_data = df_input_path, 
#   input_metadata = df_input_metadata, 
#   output = "demo_functional", 
#   fixed_effects = c("diagnosis", "dysbiosis"),
#   reference = c("diagnosis,nonIBD"),
#   min_abundance = 1,
#   min_prevalence = 0.1)
# 
# df_input_metadata$CD_dysbiosis = (df_input_metadata$diagnosis == "CD") *
#   df_input_metadata$dysbiosis
# df_input_metadata$UC_dysbiosis = (df_input_metadata$diagnosis == "UC") *
#   df_input_metadata$dysbiosis
# 
# fit_data_interaction = Maaslin2(
#   input_data = df_input_data, 
#   input_metadata = df_input_metadata, 
#   output = "demo_output_interaction", 
#   fixed_effects = c("diagnosis", "dysbiosis", "CD_dysbiosis", "UC_dysbiosis"), 
#   reference = c("diagnosis,nonIBD"),
#   max_significance = 0.05,
#   min_prevalence = 0.1,
#   min_abundance = 1)