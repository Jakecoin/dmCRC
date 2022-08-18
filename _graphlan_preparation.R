
# prepare for graphlan
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202201DCRC")
unigene_taxonomy <- read.table("species_sig_taxonomy.txt")

unigene_taxonomy <- unigene_taxonomy[unigene_taxonomy$domain %in% c("Bacteria"),] #, "Eukaryota"
for(i in 4:10){  unigene_taxonomy[,i] <- gsub('\\.', replacement = '', unigene_taxonomy[,i])}

unigene_taxonomy$name <- apply(unigene_taxonomy[,4:10], 1,paste, collapse='.')
out <- apply(unigene_taxonomy[,4:10], 1,paste, collapse='.')

as.data.frame(out)[,1]
# write.table(as.data.frame(out)[,1], "/Users/lijinming/Downloads/guide.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("Normal", "Overweight", "Obesity")
index <- 1

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]

mat <- read.table(paste(target[index], ".txt", sep = ""))
rownames(mat) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(mat))

mat <- mat[rownames(unigene_taxonomy), rownames(meta)]

{
  ndiff <- read.table(paste("Result/", save_path[index], "/", group.test[1], "_maaslin2.txt", sep = ""), 
                      header = T, row.names = 1)
  ovdiff <- read.table(paste("Result/", save_path[index], "/", group.test[2], "_maaslin2.txt", sep = ""), 
                       header = T, row.names = 1)
  obdiff <- read.table(paste("Result/", save_path[index], "/", group.test[3], "_maaslin2.txt", sep = ""), 
                       header = T, row.names = 1)
  total <- read.table("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome/Result/species_maaslin2/whole_maaslin2.txt", 
                      header = T, row.names = 1)
  num <- 6 # 6 pvalue / 8 fdr
  threshold = 0.05
  table(ndiff[,num] < threshold)
  table(ovdiff[,num] < threshold)
  table(obdiff[,num] < threshold)
  
  ndiff <- ndiff[ndiff[,num] < threshold,]$feature
  ovdiff <- ovdiff[ovdiff[,num] < threshold,]$feature
  obdiff <- obdiff[obdiff[,num] < threshold,]$feature
  total <- total[total[,num] < threshold,]$feature
  heatmap_out <- data.frame(
    item = rownames(mat),
    Abundance = rep("#D00000", nrow(mat)),
    Normal = apply(mat,1,function(x){log2(median(x[meta$multigroup=="CN"])/median(x[meta$multigroup=="HN"] ))}), # +1
    Overweight = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COv"] )/median(x[meta$multigroup=="HOv"] ))}),
    Obesity = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COb"] )/median(x[meta$multigroup=="HOb"] ))}),
    Total = apply(mat,1,function(x){log2(median(x[meta$Cancer=="CRC"] )/median(x[meta$Cancer=="Health"] ))})
  )
  heatmap_out <- heatmap_out[,-1]
  rownames(heatmap_out) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(heatmap_out))
  
  ndiff <- convert_feature(ndiff)
  ovdiff <- convert_feature(ovdiff)
  obdiff <- convert_feature(obdiff)
  total <- convert_feature(total)
  heatmap_out <- na.omit(heatmap_out)
  
  Abundance <- log10(apply(mat, 1, median)+1)
  Abundance <- Abundance/max(Abundance)
  summary(Abundance)
  
  sig_out <- data.frame(
    item = rownames(heatmap_out),
    Abundance = Abundance,
    Normal = ifelse(rownames(heatmap_out) %in% ndiff, 1, 0),
    Overweight = ifelse(rownames(heatmap_out) %in% ovdiff, 1, 0),
    Obesity = ifelse(rownames(heatmap_out) %in% obdiff, 1, 0),
    Total = ifelse(rownames(heatmap_out) %in% total, 1, 0)
  )
  
  rownames(sig_out) <- sig_out[,1]
  sig_out <- sig_out[,-1]
  
  sig_out2 <- data.frame(
    item = rownames(heatmap_out),
    Abundance = Abundance,
    Normal = apply(mat,1,function(x){log2(median(x[meta$multigroup=="CN"] )/median(x[meta$multigroup=="HN"] ))}), # +1
    Overweight = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COv"] )/median(x[meta$multigroup=="HOv"] ))}),
    Obesity = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COb"] )/median(x[meta$multigroup=="HOb"] ))}),
    Total = apply(mat,1,function(x){log2(median(x[meta$Cancer=="CRC"] )/median(x[meta$Cancer=="Health"] ))})
  )
  
  rownames(sig_out2) <- sig_out2[,1]
  sig_out2 <- sig_out2[,-1]
  sig_out2 <- abs(sig_out2)
}

for(i in 2:5) {
  heatmap_out[,i] <- ifelse(sig_out[,i]==0, "#FFFFFF", ifelse(heatmap_out[,i]>0, "#D00000", "#009FFD"))
}

table(heatmap_out$Obesity)

rownames(heatmap_out) == rownames(unigene_taxonomy)
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

out <- rbind(out_color, out_alpha) # is.alpha
out <- out %>% arrange(name, group)
out <- as.data.frame(out)
write.table(out, "/Users/lijinming/Downloads/1_annot_heatmap_total_without_alpha.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

vspe <- read.table("final_selected_species.txt")
vspe[,1]
vspe_out <- data.frame(
  name = rep(unigene_taxonomy[vspe[,1],]$name, each = 3),
  ring = rep(c("ring_color", "ring_alpha", "ring_shape"), nrow(vspe)),
  layer = rep(6, nrow(vspe) * 3),
  arg = rep(c("#000000", "1", "v"), nrow(vspe))
)

write.table(vspe_out, "/Users/lijinming/Downloads/2_annot_vspe_total.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# cat(vspe_out, file="outfile.txt", append=TRUE)

anno_out <- data.frame(
  name = unigene_taxonomy[vspe[,1],]$spe,
  ring = rep("annotation", nrow(vspe)),
  arg = rep(c("*:*"), nrow(vspe))
)

write.table(anno_out, "/Users/lijinming/Downloads/3_annot_letter.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# cat(anno_out, file="outfile.txt", append=TRUE)
# mylist = list(out, vspe_out, anno_out)
# lapply(mylist, write, "test.txt", append=TRUE)
# sink("mylist.txt")
# print(out)
# sink()
