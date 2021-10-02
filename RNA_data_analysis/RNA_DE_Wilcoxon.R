###################################################################################
#                               RNA data analysis                                 #
#                 - Differential expression analysis using                        #
#                             Wilcoxon Rank Sum Test                              #
###################################################################################

library(Seurat)
library(patchwork)
library(ggplot2)
library(Nebulosa)
library(viridis)

#This source contains the function needed for the feature+violin plots.
source("DFV_plot.R")

#This function is needed just for the appropriate selection of the heatmaps height.
height_sel <- function(markers){
  if (nrow(markers) > 100){
    return(25)
  }else{
    return(7)
  }
}

#This function is used to obtain the best size ratio for the heatmaps
  
#Cluster-based and condition-based markers analysis
# - type: can be "clusters" or "conditions"
# - seurat_object: Seurat object containing the assay for the RNA analysis
# - assay: Seurat assay for the RNA analysis
# - slot: slot for the RNA analysis
# - idents_1, idents_2: type of identities used (RNA clusters/conditions) set as default for the clusters and the conditions.
#   Can be changed upon the user's interests.
# - patient: patient name

scRNA_DE_analysis <- function(type, seurat_object, assay, slot, idents_1 = "rnaClusterID", idents_2 = "orig.ident", patient){
  
  #Setting the default assay
  DefaultAssay(seurat_object) <- assay
  
  if (type == "clusters"){
    
    #First etting the identities to be used: RNA clusters
    Idents(seurat_object) <- seurat_object[[idents_1]]
    
    #Computation of all the RNA markers
    #The object RNAmarkers is a data frame containing the markers - i.e. the entries enriched 
    #in each of the clusters. In particular, only genes that are detected in 25% of the cells,
    #that are positive markers, and which show, on average, at least 0.25 log-fold difference.
    #The test used is Wilcoxon Rank Sum test.
    RNAmarkers <- FindAllMarkers(seurat_object, assay = assay, 
                                 slot = slot, min.pct = 0.25, only.pos = TRUE, 
                                 logfc.threshold = 0.25, test.use = "wilcox")
    
    #Subdivision of the markers by clusters,
    #Ordering by significance and average difference in expression
    #And final filtering by adjusted p-value
    
    for (i in Idents(seurat_object)){
      markers <- RNAmarkers[RNAmarkers$cluster == i,]
      markers <- markers[order(markers$p_val_adj, decreasing = FALSE),]
      markers <- markers[order(markers$avg_diff, decreasing = TRUE),]
      markers <- markers[markers$p_val_adj < 0.001,]
      assign(paste0("markers_cluster_",i), markers, envir = .GlobalEnv)
    }
    
    #Top 10 markers per each cluster
    top10 <- c()
    for (markers_set in grep("^markers_cluster",names(.GlobalEnv),value=TRUE)){
      top10 <- rbind(top10, head(eval(as.symbol(markers_set)), n = 10))
    }
    
    #Plotting - Heatmap
    #Top 10 differentially expressed genes per each cluster - Wilcoxon Rank Sum test DE
    ggsave(DoHeatmap(seurat_object, 
                     assay = "SCT", 
                     slot = "scale.data", 
                     features = top10[order(top10$cluster),]$gene, 
                     group.bar = TRUE,
                     size = 3) + labs(title = patient) + theme(text = element_text(size = 8)),
           filename = paste0("RNA_top10_heatmap_markers_clusters_", patient, ".pdf"), height = 10, dpi = 700)
    
    #Plotting - Feature plots and VlnPlots
    RNA_plot_clusters <- DFV_plot(seurat_object, 
                                  "scale.data", 
                                  top10$gene, 
                                  paste0(patient," Wilcoxon RNA markers"), 
                                  "RNA cluster-based analysis")
    
    ggsave(paste0("RNA_top10RNAmarkers_per_cluster_Feature&Vln_", patient,".pdf"), RNA_plot_clusters, width = 25, height = length(top10$gene)*5.36, limitsize = FALSE, dpi = 700)
    
    #Plots for single markers objects - top 10 per each one
    for (i in as.vector(grep("^markers_cluster",names(.GlobalEnv),value=TRUE))){
      
      #Markers to be visualized
      markers_features <- eval(as.symbol(i))
      
      if (nrow(markers_features) > 0){
        
      #Download of the markers as .csv file
      write.csv(eval(as.symbol(i)), paste0("RNA_", as.character(i),"_",patient,".csv"))
      
      #Extraction of the cluster number for the title
      cluster <- substring(i, nchar(i))
      
      #Plotting heatmap
      ggsave(DoHeatmap(seurat_object, 
                       assay = "SCT", 
                       slot = "scale.data", 
                       features = markers_features$gene, 
                       group.bar = TRUE,
                       size = 5) + theme(axis.text.y = element_text(size = 7), text = element_text(size = 7)),
             filename = paste0("RNA_heatmap_", i,"_",patient,".pdf"), dpi = 700, height = height_sel(markers_features), width = 9, limitsize = FALSE)
      
      #Plotting density, feature plots and the violin plot  
      plot <- DFV_plot(seurat_object, 
                       "scale.data", 
                       markers_features$gene, 
                       paste0(patient, " Wilcoxon RNA markers"), 
                       paste0("RNA cluster-based analysis - Cluster ",cluster))
      
      #Saving the plot
      ggsave(paste0("RNA_Feat&Vln_", i,"_",patient,".pdf"), 
             plot, 
             width = 25,
             height = 5.36*nrow(eval(as.symbol(i))),
             limitsize = FALSE)
      
      }
    }
    return("Done")
  } else if (type == "conditions"){
    
    #Setting the identities to be used: conditions
    Idents(seurat_object) <- seurat_object[[idents_2]]
    
    #Computation of all the RNA markers
    RNAmarkers_conditions <- FindAllMarkers(seurat_object, assay = assay, 
                                            slot = slot, min.pct = 0.25, only.pos = TRUE, 
                                            logfc.threshold = 0.25, test.use = "wilcox")
    
    #Subdivision of the markers by condition,
    #Ordering by significance and average difference in expression
    #And final filtering by adjusted p-value
    for (i in Idents(seurat_object)){
      markers <- RNAmarkers_conditions[RNAmarkers_conditions$cluster == i,]
      markers <- markers[order(markers$p_val_adj, decreasing = FALSE),]
      markers <- markers[order(markers$avg_diff, decreasing = TRUE),]
      markers <- markers[markers$p_val_adj < 0.001,]
      assign(paste0("markers_condition_",i), markers, envir = .GlobalEnv)
    }  
    #Top 10 markers per each condition
    top10_conditions <- c()
    for (markers_set in grep("^markers_condition",names(.GlobalEnv),value=TRUE)){
      top10_conditions <- rbind(top10_conditions, head(eval(as.symbol(markers_set)), n = 10))
    }
      
    #Plotting - Heatmap
    #Top 10 differentially expressed genes per each cluster - Wilcoxon Rank Sum test DE
    
    #Setting order of the conditions to be represented
    #Specific to our case. 
    
    if ("NT" %in% levels(seurat_object)){
      levels(seurat_object) <- c("NT", "AG", "R", "RAG")
    }
    
    ggsave(DoHeatmap(seurat_object, 
                     assay = "SCT", 
                     slot = "scale.data", 
                     features = top10_conditions$gene, 
                     group.bar = TRUE,
                     size = 3) + labs(title = patient) + theme(text = element_text(size = 8)),
           filename = paste0("RNA_top10_heatmap_markers_conditions_",patient, ".pdf"), height = 10, dpi = 700)
      
    #Plotting - Feature plots and VlnPlots
    RNA_plot_conditions <- DFV_plot(seurat_object, 
                                    "scale.data", 
                                    top10_conditions$gene, 
                                    paste0(patient," Wilcoxon RNA markers"), 
                                    "RNA condition-based analysis")
    
    ggsave(paste0("RNA_top10RNAmarkers_per_condition_Feature&Vln_", patient,".pdf"), RNA_plot_conditions, width = 25, height = length(top10_conditions$gene)*5.36, limitsize = FALSE, dpi = 700)
    
    #Plots for single markers objects - top 10 per each one
    for (i in as.vector(grep("^markers_condition",names(.GlobalEnv),value=TRUE))){
      
      #Markers to be visualized
      markers_features <- eval(as.symbol(i))
      
      if (nrow(markers_features) > 0){
        
      #Download of the markers as .csv file
      write.csv(eval(as.symbol(i)), paste0("RNA_", as.character(i),"_",patient,".csv"))
      
      #Extraction of the cluster number for the title
      condition <- substring(i, 19, nchar(i))
      
      #Plotting heatmap
      ggsave(DoHeatmap(seurat_object, 
                       assay = "SCT", 
                       slot = "scale.data", 
                       features = markers_features$gene, 
                       group.bar = TRUE,
                       size = 5) + theme(axis.text.y = element_text(size = 7), text = element_text(size = 7)),
             filename = paste0("RNA_heatmap_", i,"_",patient,".pdf"), dpi = 700, height = height_sel(markers_features), width = 9, limitsize = FALSE)
      
      
      #Plotting density, feature plots and the violin plot  
      plot <- DFV_plot(seurat_object, 
                       "scale.data", 
                       markers_features$gene, 
                       paste0(patient, " Wilcoxon RNA markers"), 
                       paste0("RNA condition-based analysis - Condition ", condition))
      
      #Saving the plot
      ggsave(paste0("RNA_Feat&Vln_", i,"_",patient,".pdf"), 
             plot, 
             width = 25, 
             height = 5.36*nrow(eval(as.symbol(i))), 
             limitsize = FALSE)
      }
    } 
    return("Done")
  }
  
  else{return("Invalid input type")}
}
