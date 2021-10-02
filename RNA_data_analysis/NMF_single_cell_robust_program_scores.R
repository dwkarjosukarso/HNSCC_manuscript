#################################
#       NMF data analysis       #   
#     single-cell NMF robust    #
#        programs scores        #
#           function            #
#################################

#This function computes the single cell scores for NMF (robust) programs. 
#It returns a table with the cell IDs as row names and the NMF (robust) programs to which the scores refer as
#column names. The scores are computing according to the scRNA-seq counts for each gene in the programs.
#In our case, we used SCT normalized and scaled counts.
# - seurat_object: Seurat data slot that has to be considered for the 
#   calculation of the scores with genes as row names and cell IDs as column names.
# - NMF_programs: list of programs where each element corresponds to a list of genes and
#   the name of the element corresponds to the name of the NMF program or gene set. 
#   The name will be reported as column name in the matrix that the function returns.

sc_NMF_robust_programs_scores <- function(seurat_object, NMF_programs){
  all_NMF_programs_scores <- c()
  for (program in names(NMF_programs)){
    print(paste0("Calculation scores for ", program))
    
    # #1. Creation of a matrix with row names the genes of the NMF program for which we want to compute the score.
    # #   Since not all the genes in the NMF programs are present with the SCT scale.data slot, we need to find the
    # #   Intersection between those genes and the NMF robust program under study.
    filtered_std_RNA <- seurat_object[intersect(rownames(seurat_object), NMF_programs[[program]]),]
    
    # #2. Calculation of the score for all the cells: the SCT counts of the available genes
    # #   in the NMF programs are summed together and divided by the number of genes considered in the score. 
    
    #In this list there are 50 dataframes containing each the SCT count of one gene per each cell
    #The 50 datafames correspond to the 50 genes in each program
    
    score_list <- lapply(rownames(filtered_std_RNA), function(x) filtered_std_RNA[x, ]) 
    if (length(score_list) != 0){
      scores <- as.matrix(score_list[[1]])
      # #All the SCT counts are summed, leading to a final score per each cell
      for (score_matrix in 2:length(score_list)){scores <- scores + as.matrix(score_list[[score_matrix]])}
      scores <- scores/length(score_list)
      
      print(paste0(length(score_list), " genes are used to compute the score of this program."))
      rownames(scores) <- program
      # 3. The scores are saved in a matrix and returned. Each program to which the score
      #  refers will be indicated in the column name.
      all_NMF_programs_scores <- cbind(all_NMF_programs_scores, t(scores))
      cat("\n")
    }else{
      print("No genes in this gene set are included in the input Seurat object assay & data slot.")
    }
  }
  print("Done.")
  return(all_NMF_programs_scores)
}