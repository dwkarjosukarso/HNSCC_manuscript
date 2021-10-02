#########################################
#         ADT&RNA data analysis         #
#           Markers & Plots             #
#          Plotting function            #
#########################################

#The name DFP stands for Density, Feature and Violin plot. 
#To use it, set the assay of interest from the Seurat object as default assay.

#The function produces four distinct visualizations, namely
# - A density plot using the package Nebulosa
# - A feature plot using the SCT scaled data
# - A feature plot using the TMM scaled data
# - A violin plot 

#The required arguments are:
# - seurat_object: the Seurat object containing the data of interest for the markers
# - slot: the data slot from the Seurat object to be considered.
# - features: the features from the Seurat object and assay of interest to be visualized
# - title: custom title
# - subtitle: custom subtitle
# - reduction_1 and reduction_2: the two available reductions with which the user wants to visualize the data.

DFV_plot <- function(seurat_object, slot, features, title, subtitle, reduction_1 = "umap_tmm", reduction_2 = "umap"){
  
  p1 <- FeaturePlot(seurat_object, 
                    slot = slot, 
                    label = TRUE, 
                    reduction = reduction_1, 
                    features = features, 
                    min.cutoff = "q9", 
                    ncol = 1,
                    cols = c("lightgrey", "red"))
  
  p2 <- FeaturePlot(seurat_object, 
                    slot = slot, 
                    label = TRUE, 
                    reduction = reduction_2, 
                    features = features, 
                    min.cutoff = "q9", 
                    ncol = 1,
                    cols = c("lightgrey", "red"))
  
  p3 <- VlnPlot(object = seurat_object, 
                features = features, 
                ncol = 1) 
  
  p4 <- plot_density(seurat_object,
                     reduction = reduction_2,
                     combine = FALSE,
                     size = 2.5,
                     features = features)
  
  p4 <- lapply(p4, function(x){x <- x + theme(plot.title = element_text(hjust = 0.5, face = "bold"))})
  p4 <- wrap_plots(p4, ncol = 1)
  
  plot <- (p4 | p2 | p1 | p3) + plot_annotation(
    title = title,
    subtitle = subtitle,
    theme = theme(
      plot.title = element_text(size = 25, hjust = 0.5),
      plot.subtitle = element_text(size = 20, hjust = 0.5)
    )
  )
  return(plot)
}
