library(lattice);library(Formula)
library(readxl)
library(ggpubr)
library(grid)
library(vcd)
library(data.table)
library(RColorBrewer)
library(caret)
require(lattice)
require(Formula)
require(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
# library(devtools)
# devtools::install_github("jaspershen/MetNormalizer")
# library(MetNormalizer)
library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(convertid)

group.test <- c("CTRL", "T2DM", "CRC", "T2DM_CRC")

meta <- read.table("dm_meta.txt")
mat <- read.table("dm_adj_met.txt")
mat <- mat[,rownames(meta)]
metname <- read.csv("dm_metinfo.csv")
# mat$name <- rownames(mat)
mat[1:5,1:5]

# s.mat <- t(apply(mat, 1, function(x){scale(x)}))
# colnames(s.mat) <- colnames(mat)
s.mat <- mat

# i 2:4
i = 4
if(i==2) {
  VIP <- read.table("result/4_PLSDA/T2DMvsCTRL_VIP.txt")
} else if (i == 3) {
  VIP <- read.table("result/4_PLSDA/CRCvsCTRL_VIP.txt")
} else if (i == 4) {
  VIP <- read.table("result/4_PLSDA/T2DM_CRCvsCTRL_VIP.txt")
}

threshold = 1
diff_output <- data.frame(
  kegg = rownames(s.mat),
  name = metname$Metabolite[match(rownames(s.mat), metname$KEGG)],
  concentration = apply(s.mat, 1, mean),
  vip = VIP[,2][match(rownames(s.mat), VIP[,1])],
  pvalue = apply(s.mat,1,function(x){wilcox.test(unlist(x[meta$Group == group.test[i]]),unlist(x[meta$Group == "CTRL"]))$p.value}),
  log2fc = apply(s.mat,1,function(x){log2(mean(na.omit(x[meta$Group == group.test[i]]))/(mean(na.omit(x[meta$Group == "CTRL"]))))}),
  group = rep(group.test[i], nrow(s.mat))
  )  %>% arrange(pvalue) %>% mutate(fdr=p.adjust(pvalue,method = "BH")) %>%
  mutate(label=ifelse(log2fc>threshold, "UP", ifelse(log2fc<(-threshold),"DOWN","Nosig"))) %>%
  mutate(TPplotlabel = ifelse((fdr < 0.05 & abs(log2fc) > threshold & vip > 1), name, NA)) %>% # pvalue/fdr
  mutate(TPplotcolor = as.factor(paste(label,ifelse(is.na(TPplotlabel), 0, 1), sep = ""))) # 

diff_output <- diff_output[!duplicated(diff_output$TPplotlabel) | is.na(diff_output$TPplotlabel),]

# write.csv(diff_output,"Ob_met_diff_CRCvsH.csv")
x.axis.lim <- round(max(abs(diff_output[!is.na(diff_output$TPplotlabel),]$log2fc))+1)
x.axis.lim
min(abs(diff_output[!is.na(diff_output$TPplotlabel),]$concentration))

pdf(paste("result/7_vocano/vocano_met_", group.test[i], ".pdf", sep = ""), width = 8, height = 8)

ggplot(data=as.data.frame(diff_output), aes(x=log2fc, y=-log10(fdr))) + #
  geom_point(aes(color = TPplotcolor), alpha = 0.8) + # , size = log10(concentration) 
  scale_color_manual(labels = c("Nosig0"="No Significance","UP1"="Significantly Up","UP0"="Up","DOWN1"="Significantly Down","DOWN0"="Down"),
                     values = c("UP0"='#FAC0AE',"DOWN0"='#9BCFF0',"UP1"='#FA2311',"DOWN1"='#6175DB',"Nosig0" ='gray')) +
  labs(x="log2 Fold Change",  y="-log10 fdr") + #
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.position = c(0.12, 0.75)) +
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = 'transparent')) +
  geom_hline(yintercept = -log10(0.2),linetype=4, color = 'black', size = 0.5) +
  geom_vline(xintercept = c(-threshold, threshold), linetype = 4, color = 'black', size = 0.5) +
  geom_text_repel(data=diff_output, aes(x = log2fc, y = -log10(fdr), label = TPplotlabel), size=5, colour = "black") + # 
  xlim(-x.axis.lim, x.axis.lim)+ ##控制横坐标长度
  theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"))

dev.off()

write.csv(diff_output, paste("result/7_vocano/vocano_met_", group.test[i], ".csv", sep = ""))

# temp <- read.csv("_Result/4_Vocano/vocano_met_O_CRCvsO_CTRL.csv")
# temp <- temp[!is.na(temp$TPplotlabel),]
# selected <- temp[,1]
# 
# temp <- read.csv("_Result/4_Vocano/vocano_met_N_CRC.csv")
# temp <- temp[!is.na(temp$TPplotlabel),]
# selected <- c(selected, temp[,1])
# 
# temp <- read.csv("_Result/4_Vocano/vocano_met_O_CTRL.csv")
# temp <- temp[!is.na(temp$TPplotlabel),]
# selected <- c(selected, temp[,1])
# 
# selected
# write.table(unique(selected), "_Result/met_sig.txt")
