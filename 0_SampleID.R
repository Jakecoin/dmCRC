setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

crc <- read_excel("DM_DMCRC.xlsx", 
                sheet = "CRC", col_names = TRUE)
y <- read_excel("557宏基因组代谢组完善BMI.xlsx", 
                sheet = "Sheet2", col_names = TRUE)

crc$SampleID <- y$SampleID[match(crc$Index, y$Name)]
crc[,1] <- crc$SampleID
crc <- crc[,c(1,3,7:17)]
# table(is.na(crc$Sample_ID))

dm <- read_excel("DM_DMCRC.xlsx", 
                          sheet = "DM", col_names = TRUE)
dm <- dm[,c(1,3,5:17)]
dm_crc <- read_excel("DM_DMCRC.xlsx", 
                          sheet = "DM_CRC", col_names = TRUE)
dm_crc <- dm_crc[,c(1,3,7:22)]
colnames(dm_crc)
x <- merge(dm, dm_crc, by = c("Sample_ID","Group", "Age", "Sex", "BMI",
                              "Course", "Drugs", "FastingPlasmaGlucose"), all = TRUE)
# View(x)
colnames(crc)
x <- merge(x, crc, by = c("Sample_ID","Group", "Age", "Sex", "BMI",
                          "Diagnosis", "Location", "T", "N", "M",
                          "Pathogenesis", "RAS", "BRAF"), all = TRUE)
# View(x)

ctrl <- read.table("data/meta557_T2DM.txt", header = T)
ctrl <- ctrl[ctrl$Group=="CTRL",]
dim(ctrl)
ctrl <- ctrl[,c(1,2,5,6,7)]
colnames(ctrl) <- colnames(x)[1:5]

x <- merge(x, ctrl, by = c("Sample_ID","Group", "Age", "Sex", "BMI"), all = TRUE)
table(x$Group)

x$Batch <- 1
x$Batch[grep(pattern = "DM", x$Sample_ID)] <- 0

# Pre-process completed
x$metformin <- "no"
x$metformin[grep(pattern = "二甲", x$Drugs)] <- "metformin"

table(is.na(x$Sample_ID))
x <- x[!is.na(x$Sample_ID),]
write.table(x, "pre_meta.txt")

x <- x[,-c(50:53)] # 代谢物没有DM_54-57

###################################################
meta$Sex[meta$Sex == "男"] <- "Male"
meta$Sex[meta$Sex == "女"] <- "Female"
table(meta$Sex)
write.table(meta, "pre_meta.txt")
