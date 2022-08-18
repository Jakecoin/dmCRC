# devtools::install_github("Tong-Chen/ImageGP")
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library("ImageGP")
 
filename = "species" #met

  data <- paste("dm_adj_", filename, ".txt", sep = "") #species
  metadata <- "dm_meta.txt"
  input_type = "normalized_OTUtable"
  dissimilarity_index = "bray"
  k = 3
  binary_dissimilarity_index = F
  data_transform = 'auto'
  group_variable = "Group"
  color_variable = "Group"
  color_variable_order = c("CTRL", "T2DM", "CRC", "T2DM_CRC")
  shape_variable = NULL
  shape_variable_order = NULL
  size_variable = NULL
  size_variable_order = NULL
  label_variable = NULL
  label_variable_order = NULL
  legend.position = 'right'
  draw_ellipse = 'auto'
  # manual_color_vector = c("#61C9A8", "#FFB35C", "#A53860", "#7E627F")
  manual_color_vector = c("#16CAB2", "#16697A", "#FF990A", "#FF0022")
  
  title = NULL
  label_font_size = NULL
  debug = FALSE
  type = 't'
  level = 0.95
  extra_ggplot2_cmd = NULL
  check_significance = T
  check_paired_significance = T
  facet_variable = NULL
  coord_fixed = T
  
  if (class(data) == "character") {
    file <- data
    data <- read.table(data)
    # data <- fit_adjust_batch$feature_abd_adj
    # data <- sp_readTable(data, row.names = 1)
  } else if ('data.frame' %in% class(data) |
             'dist' %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }
  
  if (class(metadata) == "character") {
    metadata <- read.table(metadata)
  } else if ('data.frame' %in% class(metadata) |
             'dist' %in% class(metadata)) {
    stop("Unknown input format for `metadata` parameter.")
  }
  
  # Keep same columns of data with rows of metadata
  matchedL <- match_two_df(data, metadata, way = "col-row")
  
  data <- matchedL$df1
  metadata <- matchedL$df2
  
  # metadata_colnames <- colnames(metadata)
  # if(!sp.is.null(group_variable) && (!group_variable %in% metadata_colnames) ){
  #   stop(paste(group_variable, 'must be column names of data!'))
  # }
  #
  # if(!sp.is.null(color_variable) && (!color_variable %in% metadata_colnames) ){
  #   stop(paste(color_variable, 'must be column names of data!'))
  # }
  if (input_type == "normalized_OTUtable" &&
      (!'dist' %in% class(data))) {
    data <- t(data)
    # copy and modified from vegan::metaMDSdist
    if (data_transform == "auto") {
      xam <- max(data)
      if (xam > 50) {
        data <- sqrt(data)
      }
      if (xam > 9) {
        data <- vegan::wisconsin(data)
      }
    } else if (data_transform != "None") {
      data <- vegan::decostand(data, method = data_transform)
    }
    
    dist_matrix <- vegan::vegdist(data, method = dissimilarity_index,
                                  binary = binary_dissimilarity_index)
  } else {
    dist_matrix <- as.dist(data)
  }
  
  # 设置维数
  ndimensions = k
  num_samples <- nrow(data)
  if (ndimensions >= num_samples) {
    ndimensions <- num_samples - 1
  }
  
  # 计算pcoa
  pcoa <- cmdscale(dist_matrix, k = ndimensions, eig = T)
  
  # 整理pcoa结果
  pcoa_points <- as.data.frame(pcoa$points)
  sum_eig <- sum(pcoa$eig)
  eig_percent <- round(pcoa$eig / sum_eig * 100, 1)
  
  colnames(pcoa_points) <- paste0("PCoA", 1:ndimensions)
  
  data <- cbind(pcoa_points, metadata)
  
  data$Row.names <- row.names(metadata)
  rownames(data) <- data$Row.names
  
  data_colnames <- colnames(data)
  #print(data_colnames)
  #print(data)
  
  if (!sp.is.null(color_variable)) {
    if (sp.is.null(group_variable)) {
      group_variable =  color_variable
    }
    if (!(color_variable %in% data_colnames)) {
      stop(paste(color_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, color_variable, color_variable_order)
  }
  
  
  if (!sp.is.null(shape_variable)) {
    if (sp.is.null(group_variable)) {
      group_variable =  shape_variable
    }
    if (!(shape_variable %in% data_colnames)) {
      stop(paste(shape_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, shape_variable, shape_variable_order)
    data[[shape_variable]] <- as.factor(data[[shape_variable]])
    shapes <- generate_shapes(data, shape_variable)
  }
  
  if (!sp.is.null(size_variable)) {
    if (!(size_variable %in% data_colnames)) {
      stop(paste(size_variable, 'must be column names of data!'))
    }
    data = sp_set_factor_order(data, size_variable, size_variable_order)
  }
  
  if (!sp.is.null(group_variable) &&
      group_variable != color_variable &&
      group_variable != shape_variable) {
    if (!(group_variable %in% data_colnames)) {
      stop(paste(group_variable, 'must be column names of data!'))
    }
    #data = sp_set_factor_order(data, group_variable, group_variable_order)
  }
  
  if (draw_ellipse == 'auto') {
    if (all(table(data[[group_variable]]) > 4)) {
      draw_ellipse = "confiden ellipse"
    } else {
      # library(ggalt)
      draw_ellipse = "encircle"
    }
  }
  
  library(ggplot2)
  
  group_variable_en = sym(group_variable)
  
  if (!sp.is.null(filename)) {
    sp_writeTable(
      data,
      file = paste0(filename, ".pcoas.txt"), 
      keep_rownames = F
    )
  }
  
  # table(data$PCoA2 > -0.2)
  # data <- data[data$PCoA2 > -0.2,]
  p <- ggplot(data, aes(x = PCoA1, y = PCoA2, group = !!group_variable_en)) +
    labs(
      x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
      y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
      title = title
    ) + geom_point() 
  
  if (!sp.is.null(color_variable)) {
    color_variable_en = sym(color_variable)
    p <- p + aes(color = !!color_variable_en)
    p <-
      sp_manual_color_ggplot2(p, data, color_variable, manual_color_vector)
  }
  
  
  if (!sp.is.null(shape_variable)) {
    shape_variable_en = sym(shape_variable)
    p <- p + aes(shape = !!shape_variable_en) +
      scale_shape_manual(values = c(19, 8))
  }
  
  if (!sp.is.null(size_variable)) {
    size_variable_en = sym(size_variable)
    p <- p + aes(size = !!size_variable_en)
  }
  
  if (!sp.is.null(label_variable)) {
    # For cloud platform usages.
    if (!(label_variable %in% data_colnames)) {
      label_variable = "Row.names"
    }
    label_variable_en = sym(label_variable)
    library(ggrepel)
    p <-
      p + geom_text_repel(
        aes(label = !!label_variable_en),
        show.legend = F,
        max.overlaps = 100
      )
  }
  
  if (draw_ellipse == "encircle") {
    p <-
      p + geom_encircle(alpha = 0.2,
                        show.legend = F,
                        aes(fill = !!group_variable_en))
    p <-
      sp_manual_fill_ggplot2(p, data, group_variable, manual_color_vector)
    
  } else if (draw_ellipse == "confiden ellipse") {
    p <-
      p + stat_ellipse(
        level = 0.95,
        type = type,
        na.rm = TRUE
      )
  }
  
  if (check_significance) {
    pcoa_adonis2 <-
      adonis2(as.formula(paste("dist_matrix", "~", group_variable)),
              data = metadata,
              permutations = 5999)
    
    dispersion <-
      betadisper(dist_matrix, group = metadata[[group_variable]])
    dispersion_test <- permutest(dispersion)
    dispersion_test_p <- dispersion_test$tab$`Pr(>F)`[1]
    
    title <- paste0(
      "adonis R2: ",
      round(pcoa_adonis2$R2, 2),
      "\nadonis P-value: ",
      round(pcoa_adonis2$`Pr(>F)`,6),
      "\ndispersion P-value: ",
      round(dispersion_test_p,6)
    )
    
    if (check_paired_significance) {
      # devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
      library(pairwiseAdonis)
      pairwise.adonis <-
        pairwise.adonis(
          x = dist_matrix,
          factors = metadata[[group_variable]],
          p.adjust.m = "BH",
          reduce = NULL,
          perm = 5999
        )
      if (!sp.is.null(filename)) {
        sp_writeTable(
          pairwise.adonis,
          file = paste0(filename, ".pairwiseAdonis.txt"),
          keep_rownames = F
        )
      }
      
      tukeyHSD <- TukeyHSD(dispersion)$group
      
      if (!sp.is.null(filename)) {
        sp_writeTable(
          tukeyHSD,
          file = paste0(filename, ".pairwiseDispersionCheck.txt")
        )
      }
      
    }
  }
  
  if (coord_fixed){
    p <- p +  coord_fixed(1)
  }
  
  p
  
  pdf(paste("dm_adj_", filename, "_PCoA.pdf", sep = ""), width = 7, height = 8) #species #_excluded_PCoA2<-0.2 
  p + theme_bw() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(color = 'black',size = 15, face = 'plain'),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(color = 'black',size = 15, face = 'plain'),
      axis.title = element_text(color = 'black',size = 15, face = 'plain'),
      axis.ticks = element_line(color = 'black')
    )
  dev.off()
  # }
  
  # pdf("dm_species_PCoA_text.pdf", width = 7, height = 8)
  # p + geom_text(data = data, aes(x = PCoA1, y = PCoA2, label = Row.names), size=2, colour = "black")
  # dev.off()
  