#########################################
#          RNA data analysis            #
#        NMF robust expression          #
#         metagenes extraction          #
#########################################

#Load the required packages
library(Seurat)

#This function computes the metagenes of a cell line given 
#the NMF matrices (basis & coefficient) for a certain range of ranks.
# - path_nmf_rds: the path of the folder containing all the n NMF RDS files, one per rank in the range.
#   Each RDS file has to be a list with a basis matrix ($w_basis) for the gene scores
# - dataset_patient: a string name that will be given to all the objects generated.
#The output is a list having as elements the n ranks of the range, and for each rank, the metagenes identified.

NMFmetagenesSeurat <- function(path_nmf_rds, dataset_patient){
  
  for (file in list.files(path_nmf_rds, pattern = "^nmf", full.names = T)){
    nmf.rds <- readRDS(file)
    obj_name <- strsplit(basename(file), split = ".rds")[[1]]
    assign(paste0(dataset_patient, "_", obj_name), nmf.rds, envir = .GlobalEnv)
  }

  #Computing the metagenes:
  # - grepping all the files whose name has NMF - which in our case its the name we assigned to the RDS object
  # - assigning a specific name Rank_(no. of rank) for easier understanding of the data
  # - computing metagenes for the clusters of that rank
  
  metagenes <- list()
  for (nmf_file in grep(dataset_patient, names(.GlobalEnv), value=TRUE)){
    obj_name <- paste0("Rank_", strsplit(nmf_file, split = "_")[[1]][[4]])
    nmf <- eval(as.symbol(nmf_file))
    metagenes[[obj_name]] <- apply(nmf$w_basis, 2, function(x) names(sort(x, decreasing = T))[1:50])
    colnames(metagenes[[obj_name]]) <- seq(1, ncol(metagenes[[obj_name]]))
    for (c in colnames(metagenes[[obj_name]])){
      colnames(metagenes[[obj_name]])[colnames(metagenes[[obj_name]]) == c] <- paste0(obj_name, "_Cluster_", c)
    }
  }
  
  saveRDS(metagenes, file = paste0(dataset_patient, "_metagenes.rds"))
  return(metagenes)
}
