library(stringr)
library(data.table)
library(dplyr)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  feature <- gsub(feature, pattern = "s__", replacement = "")
  return(feature)
}

unimatch <- read.table("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome/unigene_species_match.txt")
sigspe <- read.table("species_selected_fdr0.05_log2fc_0.txt") # species_heatmap_sig_0518.txt
sigspe <- sigspe[,1]
sigspe <- gsub(sigspe, pattern = "s__", replacement = "")

select <- unimatch[sigspe,]
names(select) <- sigspe

unigene_taxonomy <- fread("/Users/lijinming/Documents/Unigenes_taxonomy.txt")
unigene_taxonomy[1:5,]
dim(unigene_taxonomy)

unigene_taxonomy <- unigene_taxonomy[select,] #match[,1]

unigene_taxonomy$domain <- str_extract(unigene_taxonomy$Taxonomy, "(?<=d\\_\\_).*(?=\\;p\\_\\_)")
unigene_taxonomy$phylum <- str_extract(unigene_taxonomy$Taxonomy, "(?<=p\\_\\_).*(?=\\;c\\_\\_)")
unigene_taxonomy$class <- str_extract(unigene_taxonomy$Taxonomy, "(?<=c\\_\\_).*(?=\\;o\\_\\_)")
unigene_taxonomy$order <- str_extract(unigene_taxonomy$Taxonomy, "(?<=o\\_\\_).*(?=\\;f\\_\\_)")
unigene_taxonomy$family <- str_extract(unigene_taxonomy$Taxonomy, "(?<=f\\_\\_).*(?=\\;g\\_\\_)")
unigene_taxonomy$genus <- str_extract(unigene_taxonomy$Taxonomy, "(?<=g\\_\\_).*(?=\\;s\\_\\_)")
unigene_taxonomy$species <- str_extract(unigene_taxonomy$Taxonomy, "(?<=s\\_\\_).*")

unigene_taxonomy <- as.data.frame(unigene_taxonomy)
rownames(unigene_taxonomy) <- unigene_taxonomy$species
table(sigspe == unigene_taxonomy$species)
unigene_taxonomy[sigspe,]

write.table(unigene_taxonomy[sigspe,], "unigene_taxonomy.txt", quote = F)

# prepare for graphlan
unigene_taxonomy <- read.table("unigene_taxonomy.txt")
table(unigene_taxonomy$domain)

sigspe <- read.table("species_selected_fdr0.05_log2fc_0.txt") # species_heatmap_sig_0518.txt
sigspe$species <- gsub(sigspe$species, pattern = "s__", replacement = "")
table(sigspe$group)

unigene_taxonomy$group <- sigspe$group[match(unigene_taxonomy$species, sigspe$species)]

# unigene_taxonomy <- unigene_taxonomy[unigene_taxonomy$domain %in% c("Bacteria"),]
for(i in 4:9){  unigene_taxonomy[,i] <- gsub('\\.', replacement = '', unigene_taxonomy[,i])}

unigene_taxonomy$name <- apply(unigene_taxonomy[,c(10, 4:9)], 1, paste, collapse='.')
out <- apply(unigene_taxonomy[,4:9], 1, paste, collapse='.')

as.data.frame(out)[,1]
write.table(as.data.frame(out)[,1], 
            "/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC/result/_Graphlan/guide.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}
target <- c("species", "met", "ko")
index <- 1
meta <- read.table("dm_meta.txt")
meta$Group <- factor(meta$Group, levels = c("CTRL", "T2DM", "CRC", "T2DM_CRC"))

mat <- read.table("dm_adj_species.txt", sep = "")
mat <- mat[, rownames(meta)]
rownames(mat) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(mat))
rownames(mat) <- gsub(rownames(mat), pattern = "s__", replacement = "")
  
mat <- mat[rownames(unigene_taxonomy), rownames(meta)]

{
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
    Abundance = rep("#D00000", nrow(mat)),
    T2DM = apply(mat,1,function(x){log2(median(x[meta$Group=="T2DM"] + 1)/median(x[meta$Group=="CTRL"] + 1))}),
    CRC = apply(mat,1,function(x){log2(median(x[meta$Group=="CRC"] + 1)/median(x[meta$Group=="CTRL"] + 1))}),
    T2DM_CRC = apply(mat,1,function(x){log2(median(x[meta$Group=="T2DM_CRC"] + 1)/median(x[meta$Group=="CTRL"] + 1))})
  )
  heatmap_out$Signature <- sigspe$group[match(heatmap_out$item, sigspe$species)]
  heatmap_out <- heatmap_out[,-1]
  rownames(heatmap_out) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(heatmap_out))
  
  Abundance <- log10(apply(mat, 1, median)+1)
  Abundance <- Abundance/max(Abundance)
  summary(Abundance)
  
  sig_out <- data.frame(
    item = rownames(heatmap_out),
    Abundance = Abundance,
    T2DM = ifelse(rownames(heatmap_out) %in% t2dm.feat, TRUE, FALSE),
    CRC = ifelse(rownames(heatmap_out) %in% crc.feat, TRUE, FALSE),
    T2DM_CRC = ifelse(rownames(heatmap_out) %in% t2dmcrc.feat, TRUE, FALSE)
  )
  
  rownames(sig_out) <- sig_out[,1]
  sig_out <- sig_out[,-1]
  
  sig_out2 <- data.frame(
    item = rownames(heatmap_out),
    Abundance = Abundance,
    T2DM = apply(mat,1,function(x){log2(median(x[meta$Group=="T2DM"] + 1)/median(x[meta$Group=="CTRL"] + 1))}),
    CRC = apply(mat,1,function(x){log2(median(x[meta$Group=="CRC"] + 1)/median(x[meta$Group=="CTRL"] + 1))}),
    T2DM_CRC = apply(mat,1,function(x){log2(median(x[meta$Group=="T2DM_CRC"] + 1)/median(x[meta$Group=="CTRL"] + 1))})
  )
  
  rownames(sig_out2) <- sig_out2[,1]
  sig_out2 <- sig_out2[,-1]
  sig_out2 <- abs(sig_out2)
}

for(i in 2:4) {
  heatmap_out[,i] <- ifelse(sig_out[,i]==0, "#FFFFFF", ifelse(heatmap_out[,i]>0, "#D00000", "#009FFD"))
}

heatmap_out[,5] <-ifelse(heatmap_out[,5] == "dmCRC_Character", "#00BD9D",
                         ifelse(heatmap_out[,5] == "dmCRC_Enhancer", "#FFC09F", "#52489C"))

table(rownames(heatmap_out) == rownames(unigene_taxonomy))
rownames(heatmap_out) = unigene_taxonomy$name
rownames(sig_out) = unigene_taxonomy$name
rownames(sig_out2) = unigene_taxonomy$name
sig_out2 <- apply(sig_out2, 2, function(x){x <- x/max(x)})

out_color <- as.data.frame(reshape::melt(as.matrix(heatmap_out)) %>% dplyr::select(name=1, group=2, value=3))
out_color$ring <- "ring_color"
out_color$group <- as.numeric(out_color$group)
out_color <- out_color[,c(1,4,2,3)]

out_alpha <- as.data.frame(reshape::melt(as.matrix(sig_out2)) %>% dplyr::select(name=1, group=2, value=3))
out_alpha$ring <- "ring_alpha"
out_alpha$group <- as.numeric(out_alpha$group)
out_alpha <- out_alpha[,c(1,4,2,3)]
out_alpha <- out_alpha[out_alpha$group == 1,] #只取 abundance

out <- rbind(out_color, out_alpha)  %>% arrange(name, group) %>% as.data.frame()
write.table(out, "/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC/result/_Graphlan/annot_1.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# vspe <- read.table("final_selected_species.txt")
vspe <- unigene_taxonomy
vspe[,1]
vspe$species

vspe_out <- data.frame(
  name = rep(unigene_taxonomy[vspe$species,]$name, each = 3), # vspe[,1]
  ring = rep(c("ring_color", "ring_alpha", "ring_shape"), nrow(vspe)),# 
  layer = rep(6, nrow(vspe) * 3), # 
  arg = rep(c("#000000", "1", "v"), nrow(vspe)) # 
)

write.table(vspe_out, "/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC/result/_Graphlan/2_annot_vspe_total.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# cat(vspe_out, file="outfile.txt", append=TRUE)

anno_out <- data.frame(
  name = unigene_taxonomy[vspe$species,]$species, # species
  ring = rep("annotation", nrow(vspe)),
  arg = rep(c("*:*"), nrow(vspe))
)

anno_phylum_out <- data.frame(
  name = unique(unigene_taxonomy$phylum), # species
  ring = rep("annotation", length(unique(unigene_taxonomy$phylum))),
  arg = rep(c("*:*"), length(unique(unigene_taxonomy$phylum)))
)

write.table(anno_phylum_out, "/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC/result/_Graphlan/annot_2.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
