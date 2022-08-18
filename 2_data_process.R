library(data.table)
library(readxl)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

# Species 数据
percentage_process <- function(data) { # 根据数据情况判定
  norm <- round(data)
  # norm <- t(data)/colSums(data,na=T)*100
  # # rowSums(norm)
  # idx <- colMeans(norm) >= 0.0001
  # HA <- norm[,idx]
  # tHA <- round(as.data.frame(t(HA)), 4)
  return(norm)
}

meta <- read.table("dm_meta.txt", header = T) # data/dm_meta.txt
meta$Sample_ID <- rownames(meta)

dmmat <- read_excel("data/rawdm_species.xlsx", sheet = "Sheet1", col_names = TRUE)
temp <- as.data.frame(dmmat)
rownames(temp) <- temp$Species
temp <- temp[,-1]
temp[1:5,1:5]

mat <- as.data.frame(fread("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome/processed_data/raw/Species_all.txt"))
  # read.table("data/Species_obesity.txt")
mat[1:5,1:5]
rownames(mat) <- mat$V1
mat <- mat[,-1]
table(colnames(mat) %in% meta$Sample_ID)
mat <- mat[, colnames(mat) %in% meta$Sample_ID]

table(rownames(temp) %in% rownames(mat))
table(rownames(mat) %in% rownames(temp))

temp <- temp[rownames(temp) %in% rownames(mat),]
mat <- mat[rownames(mat) %in% rownames(temp),]
mat <- mat[rownames(temp),]

mat <- cbind(mat, temp)
mat <- mat[,meta$Sample_ID]
mat[1:5,1:5]

mat <- as.data.frame(apply(mat, 2, function(x) x/sum(x))) * 1000000
mat[1:5,1:5]

fwrite(mat, file = paste("dm_species.txt", sep=""), sep = "\t", quote = FALSE, row.names = T)

############
# KO gene
dmmat <- read_excel("data/rawdm_ko.xlsx", sheet = "Sheet1", col_names = TRUE)
temp <- as.data.frame(dmmat)
rownames(temp) <- temp[,1]
temp <- temp[-1, -1]
temp[1:5,1:5]

mat <- as.data.frame(fread("/Users/lijinming/Documents/datamove/function_result/KEGG/4_KOEntry/KEGG_KOEntry_abund.txt")) #data/ko_obesity.txt
mat <- as.data.frame(mat)
mat[1:5,1:5]
mat <- mat[-1,]
rownames(mat) <- mat[,1]
mat <- mat[,-1]
table(colnames(mat) %in% meta$Sample_ID)
mat <- mat[, colnames(mat) %in% meta$Sample_ID]

table(rownames(temp) %in% rownames(mat))
table(rownames(mat) %in% rownames(temp))

temp <- temp[rownames(temp) %in% rownames(mat),]
mat <- mat[rownames(mat) %in% rownames(temp),]
mat <- mat[rownames(temp),]
mat <- cbind(mat, temp)

dim(mat)
mat <- mat[,meta$Sample_ID]
mat <- as.data.frame(apply(mat, 2, function(x) x/sum(x))) * 1000000

fwrite(mat, file = paste("dm_ko.txt", sep=""), sep = "\t",quote = FALSE, row.names = T)

##########################################################
# Genus
meta <- read.table("pre_meta.txt", header = T)
dmmat <- read_excel("data/rawdm_genus.xlsx", sheet = "Sheet1", col_names = TRUE)
temp <- as.data.frame(dmmat)
rownames(temp) <- temp[,1]
temp <- temp[,-1]
temp[1:5,1:5]

mat <- as.data.frame(fread("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome/taxonomy_data/taxonomy_Genus_abund.txt"))
# read.table("data/Species_obesity.txt")
mat[1:5,1:5]
rownames(mat) <- mat[,1]
mat <- mat[,-1]
table(colnames(mat) %in% meta$Sample_ID)
mat <- mat[, colnames(mat) %in% meta$Sample_ID]

table(rownames(temp) %in% rownames(mat))
table(rownames(mat) %in% rownames(temp))

temp <- temp[rownames(temp) %in% rownames(mat),]
mat <- mat[rownames(mat) %in% rownames(temp),]
mat <- mat[rownames(temp),]

mat <- cbind(mat, temp)
mat <- mat[,meta$Sample_ID]

mat <- as.data.frame(apply(mat, 2, function(x) x/sum(x))) * 1000000

fwrite(mat, file = paste("dm_genus.txt", sep=""), sep = "\t",quote = FALSE, row.names = T)

#################################################
# met
dmmat <- as.data.frame(read.csv("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC/data/summary 3/2.MetaboliteIdentification/identification.raw.intensity.csv"))
temp <- as.data.frame(dmmat)
dmmat[1:5,1:13]
write.table("dm_metinfo.txt", as.data.frame(dmmat[,1:13]))
# rownames(temp) <- temp[,1]
# temp <- temp[,-1]
temp[1:5,1:20]
table(temp$KEGG != "")
temp <- temp[temp$KEGG != "",]

# mat <- as.data.frame(read.csv("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome/metabolic_data/combine-norm-metaboAnalystInput（修正）.csv"))
mat <- as.data.frame(read_excel("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome/metabolic_data/rawintensity.xlsx",
                                sheet = "identification.raw.intensity"))
mat[1:5,1:5]
# rownames(mat) <- mat[,1]
# mat <- mat[,-1]
length(grep(mat$MS2kegg, pattern = "C"))
mat <- mat[grep(mat$MS2kegg, pattern = "C"),]
table(mat$MS2kegg %in% temp$KEGG)
table(temp$KEGG %in% mat$MS2kegg)
table(mat$MS2kegg)
table(temp$KEGG)

x <- temp[,c(2, 3, 9, 14:124)]
y <- mat[,c(2,3,10,17:635)]
colnames(y)[3] <- "KEGG"
y$RT <- round(y$RT, 3)
x <- x[,-c(1,2)]
y <- y[,-c(1,2)]

MAX <- by(x[,-1], x[,1],
          function(x) rownames(x)[which.max(rowMeans(x))]) 
MAX <- as.character(MAX) 
x <- x[rownames(x) %in% MAX,] 

MAX <- by(y[,-1], y[,1],
          function(x) rownames(x)[which.max(rowMeans(x))]) 
MAX <- as.character(MAX) 
y <- y[rownames(y) %in% MAX,] 

met <- merge(x, y, by = "KEGG", all = FALSE)
dim(met)

rownames(met) <- met[,1]
met <- met[,-1]

met <- met[,meta$Sample_ID]

met <- as.data.frame(apply(met, 2, function(x) x/sum(x))) * 10000000
met[1:5,1:5]

fwrite(met, file = "dm_met.txt", sep = "\t",quote = FALSE, row.names = T)

