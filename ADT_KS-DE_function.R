###################################################################################
#                               ADT data analysis                                 #
#        Kolmogorov-Smirnov differential expression analysis (function KS-DE)     #
###################################################################################

library(ggplot2)
library(Seurat)
library(patchwork)
library(Nebulosa)

#KS DE function ####
#Disclaimer: the following calculations are based on TMM normalization, 
#hence, the Seurat assay of choice is TMM.

# - object: Seurat object
# - cluster_1: cluster 1 for the comparison
# - cluster_2: cluster 2 for the comparison
# - slot: data slot from the Seurat object (assay TMM)
# - base: base of the logFC for the differential expression analysis

KS_DE <- function(object, cluster_1, cluster_2, slot, base) {
  DE_FC <- as.data.frame(FoldChange(
    object,
    ident.1 = cluster_1,
    ident.2 = cluster_2,
    slot = slot,
    base = base
  ))
  
  group <- as.factor(as.vector(Idents(object)))
  
  KS_p_val <- sapply(
    X = 1:nrow(object[["TMM"]]@scale.data),
    FUN = function(x) {
      return(stats::ks.test(object[["TMM"]]@scale.data[x, group == as.character(cluster_1)],
                            object[["TMM"]]@scale.data[x, group == as.character(cluster_2)],
                            exact = TRUE)$p.value)
    }
  )
  
  #Adjustment of the p-value
  KS_p_val <- p.adjust(KS_p_val, method = "BH")
  
  #Setting protein names
  names(KS_p_val) <- rownames(object[["TMM"]]@scale.data)
  
  #Creating the column for the KS p-values
  DE_FC$KS_adj_pval <- KS_p_val
  
  #Adding a column to keep track of the comparisons
  DE_FC$comparison <- rep(paste(as.character(cluster_1),
                                as.character(cluster_2),
                                sep = "vs"), nrow(DE_FC))
  
  #Filtering out downregulated genes (to find markers) and non significant results
  DE_FC <- DE_FC[DE_FC$avg_diff > 0.25 & DE_FC$KS_adj_pval < 0.005, ]
  
  return(DE_FC)
}

#Computation of the ADT markers
#The only argument is the Seurat (integrated) object from the databse of choice used.

ADT_markers_KS <- function(Seurat_split.integrated) {
  #P-value computed via Kolmogorov-Smirnov  - Cluster based
  DefaultAssay(Seurat_split.integrated) <- "TMM"
  Idents(Seurat_split.integrated) <- Seurat_split.integrated[["adtClusterID"]]
  
  #Performing all pairwise comparisons
  for (i in levels(Idents(Seurat_split.integrated))) {
    for (k in levels(Idents(Seurat_split.integrated))) {
      if (i != k) {
        #Computation of the markers with Kolmogorov-Smirnov
        markers <- KS_DE(Seurat_split.integrated, i, k, "scale.data", 2)
        
        #Saving the markers in the comparison in the environment
        name <- paste0("DE_", i, "vs", k)
        assign(name, markers, envir = .GlobalEnv)
      }
    }
  }
  
  #Generating a list of markers per each cluster: double strategy
  # - ADT_markers_cluster_X contains all the markers - i.e. features enriched in a specific cluster X
  #   with respect to all the others, as computed in the pairwise comparisons. When the same feature
  #   is recognized as enriched/positively DE in more than one comparison of the cluster X against
  #   another (e.g. KLF4 is found enriched in cluster 1 with respect to both clusters 3 and 4), the entry
  #   with the lowest adjusted p-value is retained in ADT_markers_cluster_X.
  #
  # - overlap_cluster_X contains all the features that have been found enriched in the cluster X
  #   in ALL the comparisons (Xvs{all the clusters except itself}). These will be the cluster
  #   or condition-specific markers.
  
  for (i in levels(Idents(Seurat_split.integrated))) {
    #Storing in a vector all the pairwise DE objects
    list_comparisons <-
      as.vector(grep(paste0("^DE_", i, "vs"), names(.GlobalEnv), value = TRUE))
    
    #Creating a temporary data frame variable storing DE data of the first comparison
    ADT <- eval(as.symbol(list_comparisons[1]))
    
    #Iterating through the remaining DE pairwise comparisons data frames for the cluster i under study
    for (k in list_comparisons[-1]) {
      k <- eval(as.symbol(k))
      
      #j is each ADT in the data frame k containing DE proteins with respect to the considered comparison
      for (j in rownames(k)) {
        #If the protein j is already present in the original data frame, then
        #If in the comparison under study the significance is higher, the entry is substituted
        if (j %in% rownames(ADT)) {
          if (k[j, "KS_adj_pval"] < ADT[j, "KS_adj_pval"]) {
            ADT[j, ] <- k[j, ]
          }
        }
        else{
          ADT[j, ] <- k[j, ]
        }
      }
      
      #Ordering of the markers by significance (adjusted p-value) and average difference in expression
      ADT <- ADT[order(ADT$KS_adj_pval,-ADT$avg_diff), ]
      assign(paste0("ADT_markers_cluster_", i), ADT, envir = .GlobalEnv)
      
      #Common markers
      assign(paste0("common_markers_", i),
             Reduce(intersect, lapply(list_comparisons, function(x) {
               rownames(eval(as.symbol(x)))
             })),
             envir = .GlobalEnv)
      
      #Cluster-specific DE proteins data frame creation - only if there are any overlapping markers
      common_markers <- eval(as.symbol(paste0("common_markers_", i)))
      if (length(common_markers) != 0) {
        #Columns of the future data frame per each cluster/condition: adjusted p-value and avg_diff per each comparison
        columns <- c()
        for (c in lapply(list_comparisons, function(x) {
          substring(x, 4, 7)
        })) {
          columns <- c(columns,
              paste0(c, "_KS_adj_pval"),
              paste0(c, "_avg_diff"))
        }
        
        #Creation of the preliminary matrix
        overlaps <- matrix(
            data = NA,
            nrow = length(common_markers),
            ncol = length(columns)
          )
        
        #Setting the column names
        colnames(overlaps) <- columns
        
        #Transforming the object as data frame
        overlaps <- data.frame(overlaps, check.names = FALSE)
        
        #Setting the rownames as the overlapping, cluster-specific markers
        rownames(overlaps) <- common_markers
        
        #Filling the data frame with the corresponding data in each of the comparison
        for (marker in common_markers) {
          for (comparison in list_comparisons) {
            column <- substring(comparison, 4, 7)
            comparison <- eval(as.symbol(comparison))
            overlaps[marker, paste0(column, "_KS_adj_pval")] <- comparison[marker, "KS_adj_pval"]
            overlaps[marker, paste0(column, "_avg_diff")] <- comparison[marker, "avg_diff"]
          }
        }
        
        #Assigning the object
        assign(paste0("ADT_cluster_", i,"_specific_markers"), overlaps, envir = .GlobalEnv)
      }
    }
  }
  
  #Condition-based (NT-AG-R-RAG) ####
  Idents(Seurat_split.integrated) <-
    Seurat_split.integrated[["orig.ident"]]
  
  #P-value computed via Kolmogorov-Smirnov  - Condition based ####
  for (i in levels(Idents(Seurat_split.integrated))) {
    for (k in levels(Idents(Seurat_split.integrated))) {
      if (i != k) {
        #Computation of the markers with Kolmogorov-Smirnov
        markers <- KS_DE(Seurat_split.integrated, i, k, "scale.data", 2)
        
        #Saving the markers in the comparison in the environment
        name <- paste0("DE_Conditions_", i, "vs", k)
        assign(name, markers, envir = .GlobalEnv)
      }
    }
  }
  
  #Generating a list of markers per each cluster: double strategy
  
  for (i in levels(Idents(Seurat_split.integrated))) {
    #Storing in a vector all the pairwise DE objects
    list_comparisons <- as.vector(grep(paste0("^DE_Conditions_", i, "vs"), names(.GlobalEnv), value = TRUE))
    
    #Creating a temporary data frame variable storing DE data of the first comparison
    ADT <- eval(as.symbol(list_comparisons[1]))
    
    #Iterating through the remaining DE pairwise comparisons data frames for the condition i under study
    for (k in list_comparisons[-1]) {
      k <- eval(as.symbol(k))
      
      #j is each ADT in the data frame k containing DE proteins with respect to the considered comparison
      for (j in rownames(k)) {
        #If the protein j is already present in the original data frame, then
        #If in the comparison under study the significance is higher, the entry is substituted
        if (j %in% rownames(ADT)) {
          if (k[j, "KS_adj_pval"] < ADT[j, "KS_adj_pval"]) {
            ADT[j, ] <- k[j, ]
          }
        }
        else{
          ADT[j, ] <- k[j, ]
        }
      }
      
      #Ordering of the markers by significance (adjusted p-value) and average difference in expression
      ADT <- ADT[order(ADT$KS_adj_pval,-ADT$avg_diff), ]
      assign(paste0("ADT_markers_condition_", i), ADT, envir = .GlobalEnv)
      
      #Common markers
      assign(paste0("common_markers_", i),
             Reduce(intersect, lapply(list_comparisons, function(x) {
               rownames(eval(as.symbol(x)))
             })))
      
      #Condition-specific DE proteins data frame creation - only if there are any overlapping markers
      common_markers <- eval(as.symbol(paste0("common_markers_", i)))
      if (length(common_markers) != 0) {
        #Columns of the future data frame per each cluster/condition: adjusted p-value and avg_diff per each comparison
        columns <- c()
        for (c in lapply(list_comparisons, function(x) {
          substring(x, 15, nchar(x))
        })) {
          columns <-c(columns,
              paste0(c, "_KS_adj_pval"),
              paste0(c, "_avg_diff"))
        }
        
        #Creation of the preliminary matrix
        overlaps <- matrix(
            data = NA,
            nrow = length(common_markers),
            ncol = length(columns)
          )
        
        #Setting the column names
        colnames(overlaps) <- columns
        
        #Transforming the object as data frame
        overlaps <- data.frame(overlaps, check.names = FALSE)
        
        #Setting the rownames as the overlapping markers
        rownames(overlaps) <- common_markers
        
        #Filling the data frame with the corresponding data in each of the comparison
        for (marker in common_markers) {
          for (comparison in list_comparisons) {
            column <- substring(comparison, 15, nchar(comparison))
            comparison <- eval(as.symbol(comparison))
            overlaps[marker, paste0(column, "_KS_adj_pval")] <- comparison[marker, "KS_adj_pval"]
            overlaps[marker, paste0(column, "_avg_diff")] <- comparison[marker, "avg_diff"]
          }
        }
        #Assigning the object
        assign(paste0("ADT_condition_", i,"_specific_markers"), overlaps, envir = .GlobalEnv)
      }
    }
  }
}
