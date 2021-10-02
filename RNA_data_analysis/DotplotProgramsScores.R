#################################
#       NMF data analysis       #   
#     single-cell NMF robust    #
#        programs scores        #
#    only programs from same    #
#    experimental conditions    #
#       no negative values      #
#           transformed         #
#################################

#This function produces a dot plot representing the median or the mean score for a specific expression program within RNA clusters.
# - programs: a list object containing the list of genes to be considered as programs whose name corresponds to the name of the program
# - RNAclusters: vector of RNA clusters
# - dataset_programs_scores: dataset containing the program names as first n columns where n is the number of programs and
#                            the program scores for each cell
# - no.programs: number of programs considered
# - patient_dataset: string containing the name of the dataset and of the patient to be reported in the figures as title

dotplot_programs_scores_RNAclusters <- function(programs, RNAclusters, dataset_programs_scores, no.programs, patient_dataset){
  dotplot_df_clusters <- data.frame(RNA_Cluster = RNAclusters)
  for (program in names(programs)){
    ncells <- sapply(RNAclusters, function(cluster)
      nrow(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == cluster &
                                       dataset_programs_scores[,program] > (median(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == cluster,][ ,program]) - sd(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == cluster,][ ,program])) &
                                       dataset_programs_scores[,program] < (median(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == cluster,][ ,program]) + sd(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == cluster,][ ,program])),]))
    dotplot_df_clusters[,program] <- ncells
  }
  
  dotplot_df_clusters <- melt(dotplot_df_clusters)
  dotplot_df_clusters <- dotplot_df_clusters[-(1:length(RNAclusters)),]
  dotplot_df_clusters$RNA_cluster <- rep(RNAclusters, times = no.programs)
  colnames(dotplot_df_clusters) <- c("NMF Robust Program", "Cells", "RNA Clusters")

  dotplot_df_clusters$`% Cells` <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  dotplot_df_clusters$Median <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  dotplot_df_clusters$Mean <- c(rep(NA, times = nrow(dotplot_df_clusters)))

  #1. Finding the % of cells in each cluster that contribute to a specific score
  #2. Computing the mean of each RNA cluster per each score
  #3. Computing the median of each RNA cluster per each score
  for (row in 1:nrow(dotplot_df_clusters)){
    dotplot_df_clusters[row, "% Cells"] <- dotplot_df_clusters[row, "Cells"]/(nrow(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == dotplot_df_clusters[row, "RNA Clusters"],])) * 100
    dotplot_df_clusters[row, "Mean"] <- mean(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == dotplot_df_clusters[row, "RNA Clusters"], as.character(dotplot_df_clusters[row, "NMF Robust Program"])])
    dotplot_df_clusters[row, "Median"] <- median(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == dotplot_df_clusters[row, "RNA Clusters"], as.character(dotplot_df_clusters[row, "NMF Robust Program"])])
  }

  #Significance
  dotplot_df_clusters$pvalue <- c(rep(NA, times = nrow(dotplot_df_clusters)))

  for (p in colnames(dataset_programs_scores)[1:no.programs]) {
    control <- c()
    for (c in RNAclusters) {
      #A sample of the same size of the cluster under assessment is randomly draw from scores of cells from the other clusters
      control_cluster <- sample(dataset_programs_scores[!(dataset_programs_scores$RNA_Cluster == c), p], length(dataset_programs_scores[dataset_programs_scores$RNA_Cluster == c, p]), replace = TRUE)
      control <- rbind(control, cbind(control_cluster, rep(c, times = length(control_cluster))))
    }
    control <- cbind(control, c(rep("control", times = nrow(control))))
    score_type <- strsplit(p, "_")
    colnames(control) <- c(paste0(score_type[[1]][1], "_score"), "rnaClusterID", "Type")
    control <- as.data.frame(control)
    assign(paste0("control_", p), control)

    df_clusters <- subset(dataset_programs_scores, select = c(p, "RNA_Cluster"))
    pvalue_cluster <- c()
    for (c in RNAclusters) {
      cases <- df_clusters[df_clusters$RNA_Cluster == c,]
      cases$Type <- c(rep("case", times = nrow(cases)))
      control_df <- get(x = paste0("control_", p))
      colnames(cases) <- colnames(control_df)
      df <- rbind(cases, control_df[control_df$rnaClusterID == c,])

      #The p-value is calculated via pairwise Wilcoxon rank sum test
      pvalue_cluster <- c(pvalue_cluster, pairwise.wilcox.test(as.numeric(df[,1]), df$Type, p.adjust.method = "BH")$p.value)
    }
    dotplot_df_clusters[dotplot_df_clusters$`NMF Robust Program` == p, "pvalue"] <- pvalue_cluster
  }

  dotplot_df_clusters$`RNA Clusters` <- as.factor(dotplot_df_clusters$`RNA Clusters`)

  return(ggplot(dotplot_df_clusters,
                aes(x = `NMF Robust Program`,
                    y = `RNA Clusters`,
                    color = Mean,
                    size = -log10(pvalue))) +
           geom_point() +
           scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue", midpoint = mean(dotplot_df_clusters$Mean)) +
           scale_size(range = c(3, 8)) +
           labs(x = "Robust Expression Program", y = "RNA Cluster", title = paste0(patient_dataset, " NMF Robust Programs scores\nRNA Clusters"), color = "Mean", size = "Significance") +
           theme_bw() +
           theme(
             # Hide panel borders and remove grid lines
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             # Change axis line
             axis.line = element_line(colour = "black"),
             axis.text.x = element_text(angle = 90, vjust = 0.6)))
}

#This function produces a dot plot representing the median or the mean score for a specific expression program within ADT clusters.
# - programs: a list object containing the list of genes to be considered as programs whose name corresponds to the name of the program
# - ADTclusters: vector of ADT clusters
# - dataset_programs_scores: dataset containing the program names as first n columns where n is the number of programs and
#                            the program scores for each cell
# - no.programs: number of programs considered
# - patient_dataset: string containing the name of the dataset and of the patient to be reported in the figures as title

dotplot_programs_scores_ADTclusters <- function(programs, ADTclusters, dataset_programs_scores, no.programs, patient_dataset){
  #Preparing data frame for the dot plot
  dotplot_df_clusters <- data.frame(ADT_Cluster = ADTclusters)
  for (program in names(programs)){
    ncells <- sapply(ADTclusters, function(cluster)
      nrow(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == cluster &
                                       dataset_programs_scores[,program] > (median(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == cluster,][,program]) - sd(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == cluster,][,program])) &
                                       dataset_programs_scores[,program] < (median(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == cluster,][,program]) + sd(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == cluster,][,program])),]))
    dotplot_df_clusters[,program] <- ncells
  }
  
  dotplot_df_clusters <- melt(dotplot_df_clusters)
  dotplot_df_clusters <- dotplot_df_clusters[-(1:length(ADTclusters)),]
  dotplot_df_clusters$ADT_cluster <- rep(ADTclusters, times = no.programs)
  colnames(dotplot_df_clusters) <- c("NMF Robust Program", "Cells", "ADT Clusters")
  dotplot_df_clusters$`% Cells` <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  dotplot_df_clusters$Median <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  dotplot_df_clusters$Mean <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  
  #1. Finding the % of cells in each cluster that contribute to a specific score
  #2. Computing the mean of each ADT cluster per each score
  #3. Computing the median of each ADT cluster per each score
  for (row in 1:nrow(dotplot_df_clusters)){
    dotplot_df_clusters[row, "% Cells"] <- dotplot_df_clusters[row, "Cells"]/(nrow(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == dotplot_df_clusters[row, "ADT Clusters"],])) * 100
    dotplot_df_clusters[row, "Mean"] <- mean(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == dotplot_df_clusters[row, "ADT Clusters"], as.character(dotplot_df_clusters[row, "NMF Robust Program"])])
    dotplot_df_clusters[row, "Median"] <- median(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == dotplot_df_clusters[row, "ADT Clusters"], as.character(dotplot_df_clusters[row, "NMF Robust Program"])])
  }
  
  #Significance
  dotplot_df_clusters$pvalue <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  
  for (p in colnames(dataset_programs_scores)[1:no.programs]) {
    control <- c()
    for (c in ADTclusters) {
      control_cluster <- sample(dataset_programs_scores[!(dataset_programs_scores$ADT_Cluster == c), p], length(dataset_programs_scores[dataset_programs_scores$ADT_Cluster == c, p]), replace = TRUE)
      control <- rbind(control, cbind(control_cluster, rep(c, times = length(control_cluster))))
    }
    control <- cbind(control, c(rep("control", times = nrow(control))))
    score_type <- strsplit(p, "_")
    colnames(control) <- c(paste0(score_type[[1]][1], "_score"), "adtClusterID", "Type")
    control <- as.data.frame(control)
    assign(paste0("control_", p), control)
    
    df_clusters <- subset(dataset_programs_scores, select = c(p, "ADT_Cluster"))
    pvalue_cluster <- c()
    for (c in ADTclusters) {
      cases <- df_clusters[df_clusters$ADT_Cluster == c,]
      cases$Type <- c(rep("case", times = nrow(cases)))
      control_df <- get(x = paste0("control_", p))
      colnames(cases) <- colnames(control_df)
      df <- rbind(cases, control_df[control_df$adtClusterID == c,])
      pvalue_cluster <- c(pvalue_cluster, pairwise.wilcox.test(as.numeric(df[,1]), df$Type, p.adjust.method = "BH")$p.value)
    }
    dotplot_df_clusters[dotplot_df_clusters$`NMF Robust Program` == p, "pvalue"] <- pvalue_cluster
  }
  
  dotplot_df_clusters$`ADT Clusters` <- as.factor(dotplot_df_clusters$`ADT Clusters`)
  
  return(ggplot(dotplot_df_clusters,
                aes(x = `NMF Robust Program`,
                    y = `ADT Clusters`,
                    color = Mean,
                    size = -log10(pvalue))) +
           geom_point() +
           scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue", midpoint = mean(dotplot_df_clusters$Mean)) +
           scale_size(range = c(3, 8)) +
           labs(x = "Robust Expression Program", y = "ADT Cluster", title = paste0(patient_dataset, " NMF Robust Programs scores\nADT Clusters"), color = "Mean", size = "Significance") +
           theme_bw() +
           theme(
             # Hide panel borders and remove grid lines
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             # Change axis line
             axis.line = element_line(colour = "black"),
             axis.text.x = element_text(angle = 90, vjust = 0.6)))
}

#This function produces a dot plot representing the median or the mean score for a specific expression program within the conditions
# - programs: a list object containing the list of genes to be considered as programs whose name corresponds to the name of the program
# - conditions: vector of conditions indicated as strings/characters
# - dataset_programs_scores: dataset containing the program names as first n columns where n is the number of programs and
#                            the program scores for each cell
# - no.programs: number of programs considered
# - patient_dataset: string containing the name of the dataset and of the patient to be reported in the figures as title

dotplot_programs_scores_conditions <- function(programs, conditions, dataset_programs_scores, no.programs, patient_dataset){
  #Conditions
  #Preparing data frame for the dot plot
  dotplot_df_conditions <- data.frame(Condition = conditions)
  for (program in names(programs)){
    ncells <- sapply(conditions, function(condition)
      nrow(dataset_programs_scores[dataset_programs_scores$Condition == condition &
                                       dataset_programs_scores[,program] > (median(dataset_programs_scores[dataset_programs_scores$Condition == condition,][ ,program]) - sd(dataset_programs_scores[dataset_programs_scores$Condition == condition,][,program])) &
                                       dataset_programs_scores[,program] < (median(dataset_programs_scores[dataset_programs_scores$Condition == condition,][,program]) + sd(dataset_programs_scores[dataset_programs_scores$Condition == condition,][,program])),]))
    dotplot_df_conditions[,program] <- ncells
  }
  
  dotplot_df_conditions <- melt(dotplot_df_conditions)
  dotplot_df_conditions$Condition <- rep(conditions, times = no.programs)
  colnames(dotplot_df_conditions) <- c("Condition", "NMF Robust Program", "Cells")
  dotplot_df_conditions$`% Cells` <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  dotplot_df_conditions$Median <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  dotplot_df_conditions$Mean <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  
  #1. Finding the % of cells in each cluster that contribute to a specific score
  #2. Computing the mean of each RNA cluster per each score
  #3. Computing the median of each RNA cluster per each score
  for (row in 1:nrow(dotplot_df_conditions)){
    dotplot_df_conditions[row, "% Cells"] <- dotplot_df_conditions[row, "Cells"]/(nrow(dataset_programs_scores[dataset_programs_scores$Condition == dotplot_df_conditions[row, "Condition"],])) * 100
    dotplot_df_conditions[row, "Mean"] <- mean(dataset_programs_scores[dataset_programs_scores$Condition == dotplot_df_conditions[row, "Condition"], as.character(dotplot_df_conditions[row, "NMF Robust Program"])])
    dotplot_df_conditions[row, "Median"] <- median(dataset_programs_scores[dataset_programs_scores$Condition == dotplot_df_conditions[row, "Condition"], as.character(dotplot_df_conditions[row, "NMF Robust Program"])])
  }
  
  #Significance
  dotplot_df_conditions$pvalue <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  
  programs_conditions_pvalues <- list()
  
  if ("NT" %in% conditions){
    for (program in colnames(dataset_programs_scores)[1:no.programs]){
      program_score <- subset(dataset_programs_scores, select = c(program, "Condition"))
      pvalue_condition <- c()
      #We exclude NT because it will be the set of cells to be used as comparison for a specific program with respect to the other conditions
      for (c in c("AG", "R", "RAG")){
        cases <- program_score[program_score$Condition == c,]
        cases$Type <- c(rep("case", times = nrow(cases)))
        control_scores <- program_score[program_score$Condition == "NT",]
        control_scores$Type <- c(rep("control", times = nrow(control_scores)))
        colnames(cases) <- colnames(control_scores)
        df <- rbind(cases, control_scores)
        colnames(df) <- c("Score", "Conditions", "Type")
        pvalue_condition <- c(pvalue_condition, pairwise.wilcox.test(as.numeric(df$Score), df$Type, p.adjust.method = "BH")$p.value)
      }
      #Adding p-value = 1 for NT
      pvalue_condition <- c(1, pvalue_condition)
      names(pvalue_condition) <- conditions
      programs_conditions_pvalues[[program]] <- pvalue_condition
    }
  }else if ("NTDMSO" %in% conditions){
    for (program in colnames(dataset_programs_scores)[1:no.programs]){
      program_score <- subset(dataset_programs_scores, select = c(program, "Condition"))
      pvalue_condition <- c()
      #We exclude NT because it will be the set of cells to be used as comparison for a specific program with respect to the other conditions
      for (c in c("NTAG","AGDMSO","AGAG")){
        cases <- program_score[program_score$Condition == c,]
        cases$Type <- c(rep("case", times = nrow(cases)))
        control_scores <- program_score[program_score$Condition == "NTDMSO",]
        control_scores$Type <- c(rep("control", times = nrow(control_scores)))
        colnames(cases) <- colnames(control_scores)
        df <- rbind(cases, control_scores)
        colnames(df) <- c("Score", "Conditions", "Type")
        pvalue_condition <- c(pvalue_condition, pairwise.wilcox.test(as.numeric(df$Score), df$Type, p.adjust.method = "BH")$p.value)
      }
      #Adding p-value = 1 for NT
      pvalue_condition <- c(1, pvalue_condition)
      names(pvalue_condition) <- c("NTDMSO","NTAG","AGDMSO","AGAG")
      programs_conditions_pvalues[[program]] <- pvalue_condition
    }
  }
  
  dotplot_df_conditions$pvalue <- unlist(programs_conditions_pvalues)
  dotplot_df_conditions$Condition <- factor(dotplot_df_conditions$Condition, levels = conditions)
  
  return(ggplot(dotplot_df_conditions,
                aes(x = `NMF Robust Program`,
                    y = Condition,
                    color = Mean,
                    size = -log10(pvalue))) +
           geom_point(stat = "identity") +
           scale_color_gradient2(high = "red", mid = "lightgrey", low = "blue", midpoint = mean(dotplot_df_conditions$Mean)) +
           scale_size(range = c(3, 8)) +
           labs(x = "Robust Expression Program", y = "Condition", title = paste0(patient_dataset, " NMF Robust Programs scores\nConditions"), color = "Mean", size = "Significance") +
           theme_bw() +
           theme(
             # Hide panel borders and remove grid lines
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             # Change axis line
             axis.line = element_line(colour = "black"),
             axis.text.x = element_text(angle = 90, vjust = 0.6)))
}