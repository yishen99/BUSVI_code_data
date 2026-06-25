##############################################################
################    breast cancer data    ####################
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

visualize_by_batch_no_subgroup <- function(Data, feature_ind,
                                           title_name,
                                           color_key_range = seq(0, 8, 0.8)
) {
  
  B <- length(Data)
  n_vec <- sapply(Data, nrow)
  total_n <- sum(n_vec)
  
  df <- do.call(cbind, lapply(Data, t))
  df <- df[feature_ind, ]
  
  batch_labels <- rep(1:B, times = n_vec)
  
  anno_df <- data.frame(
    Batch = batch_labels
  )
  
  gray_shades <- gray.colors(B, start = 0.75, end = 0.2)
  batch_colors <- setNames(gray_shades, 1:B)
  
  ha <- HeatmapAnnotation(
    df = anno_df,
    col = list(Batch = batch_colors),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
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



visualize_by_subgroup <- function(Data, feature_ind,
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
  
  if(is.Y==TRUE) subgroup_labels <- rep(paste0("A", 1:K), times = n_vec)
  else subgroup_labels <- rep(paste0("C", 1:K), times = n_vec)
  
  if (!is.null(batch_list)) {
    batch_vec <- unlist(batch_list)
    if (length(batch_vec) != total_n) stop("Subgroups don't match sample number.")
  } else {
    batch_vec <- rep("NA", total_n)
  }
  
  anno_df <- data.frame(
    Batch = batch_vec,
    Subgroup = subgroup_labels
  )
  
  if(is.Y==TRUE) subgroup_colors <- setNames(c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")[1:K], paste0("A", 1:K))
  else subgroup_colors <- setNames(c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")[1:K], paste0("C", 1:K))
  if (is.null(batch_colors)) {
    batch_levels <- sort(unique(batch_vec))
    n_batch <- length(batch_levels)
    gray_shades <- gray.colors(n_batch, start = 0.75, end = 0.2)
    batch_colors <- setNames(gray_shades, batch_levels)
  }
  
  ha <- HeatmapAnnotation(
    Batch = anno_df$Batch,
    Subgroup = anno_df$Subgroup,
    col = list(
      Batch = batch_colors,
      Subgroup = subgroup_colors
    ),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      Subgroup = list(
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


batch_label_by_subgroup  <- function(Z_list) {
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
load("../input_data/integrativeGEdata_microarray.RData")
Y <- Data_list
Y <- lapply(Y, t)
B <- length(Y)
G <- ncol(Y[[1]])
K <- 4 ### From BIC analysis
n_vec <- sapply(Y, nrow)


## Read data Y_BUSVI
load("../result_data/Y_BUSVI.RData")
load("../result_data/Z_BUSVI.RData")

## Read data Y_BUS
load("../result_data/Y_BUS.RData")
load("../result_data/Z_BUS.RData")

## Read running times
load("../result_data/Running_Time_BUSVI.RData")
load("../result_data/Running_Time_BUS.RData")



################# Visualization #####################


### Y by batch
png(filename="../figures/Figure3(a).png", width = 500, height = 300)
visualize_by_batch_no_subgroup(Data = Y, feature_ind = 1:G, title_name = "Raw simulated data")
dev.off()


### Y_BUSVI by batch
png(filename="../figures/Figure3(c).png", width = 500, height = 300)
visualize_by_batch_no_subgroup(Data = Y_BUSVI, feature_ind = 1:G, title_name = "BUSVI")
dev.off()


### Y_BUS by batch
png(filename="../figures/Figure3(d).png", width = 500, height = 300)
visualize_by_batch_no_subgroup(Data = Y_BUS, feature_ind = 1:G, title_name = "BUS")
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

png(filename="../figures/Figure3(e).png", width = 500, height = 400)
ggplot(times_df,
       aes(x = Method,
           y = Time,
           fill = Method)) +
  geom_boxplot(width = 0.6) +
  geom_jitter(
    width = 0.15,
    size = 2
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



### Running time trajectory diagram

times_df <- data.frame(
  Method = c(
    rep("BUSVI", length(Running_Time_BUSVI)),
    rep("BUS", length(Running_Time_BUS))
  ),
  Time = c(
    Running_Time_BUSVI,
    Running_Time_BUS
  ),
  Rep = c(
    1:length(Running_Time_BUSVI),
    1:length(Running_Time_BUS))
)

png(filename="../figures/Figure3(f).png", width = 600, height = 400)

ggplot(times_df,
       aes(x = Rep,
           y = Time,
           color = Method,
           group = Method)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  scale_color_manual(
    values = c(
      "BUSVI" = "#B40426",
      "BUS"   = "#3B4CC0"
    )
  ) +
  
  scale_x_continuous(
    breaks = c(2, 4, 6, 8, 10)
  ) +
  
  labs(
    x = "Replications",
    y = "",
    color = ""
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "right",
    legend.text = element_text(
      size = 24,
      family = "Arial"
    ),
    legend.title = element_text(
      size = 24,
      family = "Arial"
    ),
    
    axis.text = element_text(
      size = 24,
      face = "plain",
      family = "Arial"
    ),
    axis.title = element_text(
      size = 24,
      face = "plain",
      family = "Arial"
    ),
    panel.grid = element_blank()
  )

dev.off()

