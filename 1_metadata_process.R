library(lattice);library(Formula);library(ggplot2);
library(readxl)
library(tidyverse)
library(ggpubr)
library(grid)
library(vcd)
library(stringr)
library(jstable)
library(tableone)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

# metasample1 <- read_excel("DM_DMCRC.xlsx", 
#                          sheet = "DM", col_names = TRUE)
# metasample2 <- read_excel("DM_DMCRC.xlsx", 
#                           sheet = "DM_CRC", col_names = TRUE)
# 
# dm <- metasample1[,colnames(metasample1) %in% colnames(metasample2)]
# # colnames(dm) <- c("Sample", "Group", "Age", "Sex", "BMI", "Course", "Diagnosis", "Complication",
# #                   "FastingPlasmaGlucose", "HbA1c", "DrugsBefore", "Drugs", "LiverFunction", "KidneyFunction", 
# #                   "BloodLipid", "ComprehensiveMetabolicPanel")
# dmcrc <- metasample2[,colnames(metasample2) %in% colnames(metasample1)]
# dmcrc$入院前用药 <- dmcrc$控制情况
# exp <- rbind(dm, dmcrc)
# View(exp)
# 
# exp$Metformin <- ifelse(is.na(str_extract(exp$入院前用药, "二甲*")), 0, 1)
# exp$Batch <- 2
# exp$Sex <- ifelse(exp$Sex == "男", "Male", "Female")
# exp <- exp[,-c(2,4,9,11)]
# colnames(exp) <- c("SampleID", "Group", "AgeYear", "Gender", "BMI", "Course", "FPG", "Metformin", "Batch")
# 
# metasample3 <- read.table("data/Sample_obesity.txt", header = 1)
# metasample3 <- metasample3[, -c(3,4,8:19)]
# colnames(metasample3)[2] <- "Group"
# metasample3$Group <- ifelse(metasample3$Group == 1, "CRC", "CTRL")
# metasample3$Batch <- 1
# 
# crc <- sample(2, nrow(metasample3), replace = TRUE, 
#              prob = c(100, nrow(metasample3)-100))
# 
# summary(metasample3$BMI[crc == 1])
# table(metasample3$Group[crc == 1])
# 
# output <- dplyr::bind_rows(exp, metasample3[crc == 1,])
# table(output$Group)
# 
# summary(output$BMI[output$Batch==1])
# summary(output$BMI[output$Batch==2])
# 
# write.table(output, "data/dm_meta.txt", quote = F, row.names = F)

#######################
# groupinfo <- group_input()
# 
# myvars <- c("Gender", "AgeYear", "Location", "T", "N", "M", "CA199",
#             "CEA", "Differentiation", "Nerve_Invasion", "Lymphatic_Invasion",
#             "Vascular_Invasion", "Mutation", "Neoadjuvant", "MicroSatellite")
# 
# catvars <- c("Cancer","Group","Gender", "Location", "Differentiation","Size", "T", "N", "M", "Mutation", "Neoadjuvant", "MicroSatellite",
#              "CEA", "CA199", "Lymphatic_Invasion", "Vascular_Invasion", "Nerve_Invasion")
# 
# CRCtable <- CreateTableOne2(vars = myvars, strata = "Cancer", data = rawd$Group,
#                                factorVars = catvars, includeNA = TRUE, test = TRUE,
#                                showAllLevels = T, showpm = TRUE, Labels = F, printToggle = T,
#                                pDigits = 4)
# 
# out <- print(CRCtable)
# write.csv(out, file = "ObesityTable-total.csv")


