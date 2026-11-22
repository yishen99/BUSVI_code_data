##############################################################
################     Simulation Case III    ##################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(ggplot2)

font.size = 20


################# Visualization Functions #####################

visualize_by_subtype_no_batch <- function(Data, feature_ind,
                                           title_name,
                                           color_key_range = seq(0, 8, 0.8), 
                                           is.Y = NULL,
                                           subtype_list = NULL,
                                           batch_colors = NULL) {
  
  K <- length(Data)
  n_vec <- sapply(Data, nrow)
  total_n <- sum(n_vec)
  
  df <- do.call(cbind, lapply(Data, t))
  df <- df[feature_ind, ]
  
  if(is.Y==TRUE) subtype_labels <- rep(paste0("A", 1:K), times = n_vec)
  else subtype_labels <- rep(paste0("C", 1:K), times = n_vec)
  
  if (!is.null(subtype_list)) {
    batch_vec <- unlist(subtype_list)
    if (length(batch_vec) != total_n) stop("Subtypes don't match sample number.")
  } else {
    batch_vec <- rep("NA", total_n)
  }
  
  anno_df <- data.frame(
    Batch = batch_vec,
    Subtype = subtype_labels
  )
  
  if(is.Y==TRUE) subtype_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0")[1:K], paste0("A", 1:K))
  else subtype_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0")[1:K], paste0("C", 1:K))
  if (is.null(batch_colors)) {
    batch_levels <- unique(batch_vec)
    n_batch <- length(batch_levels)
    gray_shades <- gray.colors(n_batch, start = 0.75, end = 0.2)
    batch_colors <- setNames(gray_shades, batch_levels)
  }
  
  ha <- HeatmapAnnotation(
    Batch = anno_df$Batch,
    Subtype = anno_df$Subtype,
    col = list(
      Batch = batch_colors,
      Subtype = subtype_colors
    ),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      Subtype = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      ),
      Batch = list(
        legend = FALSE
      )
    )
  )
  
  ht <- Heatmap(df,
                name = "Value",
                col = colorRamp2(color_key_range, colorRampPalette(c("#3B4CC0", "white", "#B40426"))(length(color_key_range))),
                top_annotation = ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                column_title = title_name,
                heatmap_legend_param = list(
                  title = "Value",
                  title_gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = font.size),
                  labels_gp = gpar(fontsize = 16),
                  legend_height = unit(4, "cm")
                )
  )
  draw(ht, heatmap_legend_side = "left")
  decorate_heatmap_body("Value", {
    grid.rect(gp = gpar(col = "black", lwd = 1, fill = NA))
  })
}




visualize_by_batch <- function(Data, feature_ind,
                               title_name,
                               color_key_range = seq(0, 8, 0.8),
                               is.Y = NULL,
                               subtype_list = NULL,
                               subtype_colors = NULL) {
  
  B <- length(Data)
  K <- length(unique(unlist(subtype_list)))
  n_vec <- sapply(Data, nrow)
  total_n <- sum(n_vec)
  
  df <- do.call(cbind, lapply(Data, t))
  df <- df[feature_ind, ]
  
  batch_labels <- rep(1:B, times = n_vec)
  if(is.Y==TRUE) subtype_labels <- paste0("A", unlist(subtype_list))
  else subtype_labels <- paste0("C", unlist(subtype_list))
  
  anno_df <- data.frame(
    Batch = batch_labels,
    Subtype = subtype_labels
  )
  
  if (is.null(subtype_colors)) {
    subtype_levels <- unique(unlist(subtype_list))
    n_subtype <- length(subtype_levels)
    if(is.Y==TRUE) subtype_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0")[1:n_subtype], paste0("A", 1:n_subtype))
    else subtype_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0")[1:n_subtype], paste0("C", 1:n_subtype))
  }
  gray_shades <- gray.colors(B, start = 0.75, end = 0.2)
  batch_colors <- setNames(gray_shades, 1:B)
  
  ha <- HeatmapAnnotation(
    df = anno_df,
    col = list(Batch = batch_colors, Subtype = subtype_colors),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      Subtype = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      ),
      Batch = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      )
    )
  )
  
  ht <- Heatmap(df,
                name = "Value",
                col = colorRamp2(color_key_range, colorRampPalette(c("#3B4CC0", "white", "#B40426"))(length(color_key_range))),
                top_annotation = ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                column_title = title_name,
                heatmap_legend_param = list(
                  title = "Value",
                  title_gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = font.size),
                  labels_gp = gpar(fontsize = 16),
                  legend_height = unit(4, "cm")
                )
  )
  draw(ht, heatmap_legend_side = "left")
  decorate_heatmap_body("Value", {
    grid.rect(gp = gpar(col = "black", lwd = 1, fill = NA))
  })
}


visualize_by_subtype <- function(Data, feature_ind,
                                  title_name,
                                  color_key_range = seq(0, 8, 0.8), 
                                  is.Y = NULL,
                                  batch_list = NULL,
                                  batch_colors = NULL) {
  
  K <- length(Data)
  n_vec <- sapply(Data, nrow)
  total_n <- sum(n_vec)
  
  df <- do.call(cbind, lapply(Data, t))
  df <- df[feature_ind, ]
  
  if(is.Y==TRUE) subtype_labels <- rep(paste0("A", 1:K), times = n_vec)
  else subtype_labels <- rep(paste0("C", 1:K), times = n_vec)
  
  if (!is.null(batch_list)) {
    batch_vec <- unlist(batch_list)
    if (length(batch_vec) != total_n) stop("Subtypes don't match sample number.")
  } else {
    batch_vec <- rep("NA", total_n)
  }
  
  anno_df <- data.frame(
    Batch = batch_vec,
    Subtype = subtype_labels
  )
  
  if(is.Y==TRUE) subtype_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0")[1:K], paste0("A", 1:K))
  else subtype_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0")[1:K], paste0("C", 1:K))
  if (is.null(batch_colors)) {
    batch_levels <- sort(unique(batch_vec))
    n_batch <- length(batch_levels)
    gray_shades <- gray.colors(n_batch, start = 0.75, end = 0.2)
    batch_colors <- setNames(gray_shades, batch_levels)
  }
  
  ha <- HeatmapAnnotation(
    Batch = anno_df$Batch,
    Subtype = anno_df$Subtype,
    col = list(
      Batch = batch_colors,
      Subtype = subtype_colors
    ),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      Subtype = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      ),
      Batch = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      )
    )
  )
  
  ht <- Heatmap(df,
                name = "Value",
                col = colorRamp2(color_key_range, colorRampPalette(c("#3B4CC0", "white", "#B40426"))(length(color_key_range))),
                top_annotation = ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                column_title = title_name,
                heatmap_legend_param = list(
                  title = "Value",
                  title_gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = font.size),
                  labels_gp = gpar(fontsize = 16),
                  legend_height = unit(4, "cm")
                )
  )
  draw(ht, heatmap_legend_side = "left")
  decorate_heatmap_body("Value", {
    grid.rect(gp = gpar(col = "black", lwd = 1, fill = NA))
  })
}


batch_label_by_subtype  <- function(Z_list) {
  batch_label = list()
  for(k in 1:K){
    temp = c()
    for(b in 1:B){
      if(length(which(Z_list[[b]]==k))>0) temp = c(temp, rep(b, times = length(which(Z_list[[b]]==k))))
    }
    batch_label[[k]] = temp
  }
  return(batch_label)
}


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y
load("../input_data/Y.RData")
load("../input_data/Z_real.RData")
B <- length(Y)
G <- ncol(Y[[1]])
K <- 3 ### From BIC analysis
n_vec <- sapply(Y, nrow)

## Read true subtype mean
load("../input_data/truemean_by_subtype.RData")

## Read data Y_BUSVI
load("../result_data/Y_BUSVI.RData")
load("../result_data/Z_BUSVI.RData")

## Read data Y_BUS
load("../result_data/Y_BUS.RData")
load("../result_data/Z_BUS.RData")
# Reorder the clustering labels (without affecting the ARI result)
map <- c(2, 1, 3)
for(b in 1:B){
  Z_BUS[[b]] <- map[Z_BUS[[b]]]
}

## Read running times
load("../result_data/Running_Time_BUSVI.RData")
load("../result_data/Running_Time_BUS.RData")



######## Raw Data and Corrected Data Organized by Subtypes #########

Y_by_subtype <- list()
for(k in 1:K){
  temp = Y[[1]][which(Z_real[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y[[b]][which(Z_real[[b]]==k), ])
  }
  Y_by_subtype[[k]] = temp
}

YBUSVI_by_subtype <- list()
for(k in 1:K){
  temp = Y_BUSVI[[1]][which(Z_BUSVI[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y_BUSVI[[b]][which(Z_BUSVI[[b]]==k), ])
  }
  YBUSVI_by_subtype[[k]] = temp
}

YBUS_by_subtype <- list()
for(k in 1:K){
  temp = Y_BUS[[1]][which(Z_BUS[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y_BUS[[b]][which(Z_BUS[[b]]==k), ])
  }
  YBUS_by_subtype[[k]] = temp
}



################# Visualization #####################

## True subtype mean
png(filename="../figures/Figure4(a).png", width = 500, height = 300)
visualize_by_subtype_no_batch(Data = truemean_by_subtype, feature_ind = 1:G, title_name = "True subtype means", is.Y = TRUE)
dev.off()


### Y by subtype
png(filename="../figures/Figure4(b).png", width = 500, height = 300)
visualize_by_subtype(Data = Y_by_subtype, feature_ind = 1:G, title_name = "Raw simulated data", is.Y = TRUE, batch_list = batch_label_by_subtype(Z_real))
dev.off()


### Y_BUSVI by subtype
png(filename="../figures/Figure4(d).png", width = 500, height = 300)
visualize_by_subtype(Data = YBUSVI_by_subtype, feature_ind = 1:G, title_name = "BUSVI", is.Y = FALSE, batch_list = batch_label_by_subtype(Z_BUSVI))
dev.off()


### Y_BUS by subtype
png(filename="../figures/Figure4(e).png", width = 500, height = 300)
visualize_by_subtype(Data = YBUS_by_subtype, feature_ind = 1:G, title_name = "BUS", is.Y = FALSE, batch_list = batch_label_by_subtype(Z_BUS))
dev.off()


### Boxplot of running times
times_df <- data.frame(
  Method = c(
    rep("BUSVI", length(Running_Time_BUSVI)),
    rep("BUS", length(Running_Time_BUS))
  ),
  Time = c(
    Running_Time_BUSVI,
    Running_Time_BUS
  )
)

png(filename="../figures/Figure4(f).png", width = 500, height = 400)
ggplot(times_df,
       aes(x = Method,
           y = Time,
           fill = Method)) +
  geom_boxplot(width = 0.6,
               fill = "white",
               color = "black"
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    color = "black"
  ) +
  labs(
    x = "",
    y = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text=element_text(
      size=24,
      face="plain",
      family="Arial"
    ),
    axis.title=element_text(
      size=24,
      face="plain",
      family="Arial"
    ),
    panel.grid = element_blank() 
  )
# boxplot(times_df)

dev.off()






