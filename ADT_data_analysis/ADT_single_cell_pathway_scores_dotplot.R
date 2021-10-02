#################################
#       ADT data analysis       #
#   single-cell pathway scores  #
#           dot plots           #
#################################

# This function calculates the % of cells of each treatment group (condition) contributing to 
# The score of a certain signalling pathway, corresponding to the number of cells in a condition having 
# a score that is included within the median + or - the standard deviation of the score for that pathway
# in the condition of interest. Then, it assesses the mean and the median pathway scores for each pathway
# and treatment group. Based on these two metrics, it calculates the significance of the scores for each condition
# as compared to the NT group using the pairwise Wilcoxon rank sum test. 
# The results are delivered through a dot plot.

dotplot_pathway_scores_conditions <- function(pathway_scores, conditions = c("NT", "AG", "R", "RAG"), patient) {
  dotplot_df_conditions <- data.frame(Condition = conditions)
  for (pathway in colnames(pathway_scores)[1:9]){
    #The cells contributing of a specific score are computed as the cells whose score is included within median +- SD of
    #The scores for a specific pathway, in a specific condition.
    ncells <- sapply(conditions, function(condition)
      nrow(pathway_scores[pathway_scores$orig.ident == condition &
                            pathway_scores[,pathway] > (median(pathway_scores[pathway_scores$orig.ident == condition,][ ,pathway]) - sd(pathway_scores[pathway_scores$orig.ident == condition,][ ,pathway])) &
                            pathway_scores[,pathway] < (median(pathway_scores[pathway_scores$orig.ident == condition,][ ,pathway]) + sd(pathway_scores[pathway_scores$orig.ident == condition,][ ,pathway])),]))
    dotplot_df_conditions[,pathway] <- ncells
  }
  
  dotplot_df_conditions <- melt(dotplot_df_conditions)
  colnames(dotplot_df_conditions) <- c("Condition", "Pathway", "Cells")
  
  dotplot_df_conditions$`% Cells` <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  dotplot_df_conditions$Median <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  dotplot_df_conditions$Mean <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  
  #1. Finding the % of cells in each condition that contribute to a specific score
  #2. Computing the mean of each condition per each score
  #3. Computing the median of each condition per each score
  for (row in 1:nrow(dotplot_df_conditions)){
    dotplot_df_conditions[row, "% Cells"] <- dotplot_df_conditions[row, "Cells"]/(nrow(pathway_scores[pathway_scores$orig.ident == dotplot_df_conditions[row, "Condition"],])) * 100
    dotplot_df_conditions[row, "Mean"] <- mean(pathway_scores[pathway_scores$orig.ident == dotplot_df_conditions[row, "Condition"], as.character(dotplot_df_conditions[row, "Pathway"])])
    dotplot_df_conditions[row, "Median"] <- median(pathway_scores[pathway_scores$orig.ident == dotplot_df_conditions[row, "Condition"], as.character(dotplot_df_conditions[row, "Pathway"])])
  }
  
  #Significance
  dotplot_df_conditions$pvalue <- c(rep(NA, times = nrow(dotplot_df_conditions)))
  pathways_conditions_pvalues <- list()
  
  for (pathway in colnames(pathway_scores)[1:9]){
    pathway_score <- subset(pathway_scores, select = c(pathway, "orig.ident"))
    pvalue_condition <- c()
    #We exclude NT (the control condition has always to be the first element in the conditions' list) because it will be the set of cells to be used as comparison for a specific pathway with respect to the other conditions
    for (c in conditions[-1]){
      cases <- pathway_score[pathway_score$orig.ident == c,]
      cases$Type <- c(rep("case", times = nrow(cases)))
      control_scores <- pathway_score[pathway_score$orig.ident == conditions[1],]
      control_scores$Type <- c(rep("control", times = nrow(control_scores)))
      colnames(cases) <- colnames(control_scores)
      df <- rbind(cases, control_scores)
      colnames(df) <- c("Score", "Conditions", "Type")
      pvalue_condition <- c(pvalue_condition, pairwise.wilcox.test(as.numeric(df$Score), df$Type, p.adjust.method = "BH")$p.value)
    }
    #Adding p-value = 1 for NT
    pvalue_condition <- c(1, pvalue_condition)
    names(pvalue_condition) <- conditions
    pathways_conditions_pvalues[[pathway]] <- pvalue_condition
  }
  
  dotplot_df_conditions$pvalue <- unlist(pathways_conditions_pvalues)
  dotplot_df_conditions$Condition <- factor(dotplot_df_conditions$Condition, levels = conditions)
  
  #Nice names for the pathways
  #Function for first letter of a string capitalized
  capFirst <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
  }
  dotplot_df_conditions$Pathway <- sapply(dotplot_df_conditions$Pathway, function(x) capFirst(strsplit(as.character(x), "_")[[1]][1]))
  
  plot <- ggplot(dotplot_df_conditions,
                 aes(
                   x = Pathway,
                   y = Condition,
                   color = Median,
                   size = -log10(pvalue)
                 )) +
    geom_point() +
    scale_colour_gradient2(high = "red", mid = "lightgrey", low = "blue", limits = c(min(dotplot_df_conditions$Median), max(dotplot_df_conditions$Median))) +
    scale_size(range = c(3, 8)) +
    labs(x = "Pathway", y = "Condition", title = patient, color = "Median", size = "Significance") +
    theme_bw() +
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.6)
    )
  return(plot)
}

# Cluster-based pathway scores dotplot with custom list of numeric clusters
# This function calculates the % of cells of each (ADT) cluster contributing to 
# The score of a certain signalling pathway, corresponding to the number of cells in a condition having 
# a score that is included within the median + or - the standard deviation of the score for that pathway
# in the cluster of interest. Then, it assesses the mean and the median pathway scores for each pathway
# and cluster. Based on these two metrics, it calculates the significance of the scores for each cluster
# as compared to a set of randomly sampled scores from the other clusters but for the same pathway. 
# The results are delivered through a dot plot.

dotplot_pathway_scores_clusters <- function(pathway_scores, clusters, patient) {
  dotplot_df_clusters <- data.frame(Cluster = clusters)
  for (pathway in colnames(pathway_scores)[1:9]){
    #The cells contributing of a specific score are computed as the cells whose score is included within median +- SD of
    #The scores for a specific pathway, in a specific cluster.
    ncells <- sapply(clusters, function(cluster)
      nrow(pathway_scores[pathway_scores$adtClusterID == cluster &
                            pathway_scores[,pathway] > (median(pathway_scores[pathway_scores$adtClusterID == cluster,][ ,pathway]) - sd(pathway_scores[pathway_scores$adtClusterID == cluster,][ ,pathway])) &
                            pathway_scores[,pathway] < (median(pathway_scores[pathway_scores$adtClusterID == cluster,][ ,pathway]) + sd(pathway_scores[pathway_scores$adtClusterID == cluster,][ ,pathway])),]))
    dotplot_df_clusters[,pathway] <- ncells
  }
  dotplot_df_clusters <- melt(dotplot_df_clusters)
  dotplot_df_clusters <- na.omit(dotplot_df_clusters[length(clusters)+1:nrow(dotplot_df_clusters),])
  dotplot_df_clusters$adtClusterID <- rep(clusters, times = 9)
  colnames(dotplot_df_clusters) <- c("Pathway", "Cells")
  dotplot_df_clusters$adtClusterID <- rep(clusters, times = 9)
  dotplot_df_clusters$`% Cells` <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  dotplot_df_clusters$Median <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  dotplot_df_clusters$Mean <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  
  #1. Finding the % of cells in each cluster that contribute to a specific score
  #2. Computing the mean of each adtClusterID per each score
  #3. Computing the median of each adtClusterID per each score
  for (row in 1:nrow(dotplot_df_clusters)){
    dotplot_df_clusters[row, "% Cells"] <- dotplot_df_clusters[row, "Cells"]/(nrow(pathway_scores[pathway_scores$adtClusterID == dotplot_df_clusters[row, "adtClusterID"],])) * 100
    dotplot_df_clusters[row, "Mean"] <- mean(pathway_scores[pathway_scores$adtClusterID == dotplot_df_clusters[row, "adtClusterID"], as.character(dotplot_df_clusters[row, "Pathway"])])
    dotplot_df_clusters[row, "Median"] <- median(pathway_scores[pathway_scores$adtClusterID == dotplot_df_clusters[row, "adtClusterID"], as.character(dotplot_df_clusters[row, "Pathway"])])
  }
  
  #Significance
  dotplot_df_clusters$pvalue <- c(rep(NA, times = nrow(dotplot_df_clusters)))
  
  #For each cluster and for each pathway score we generate a control group of scores from cells not included in the cluster of interest
  for (p in colnames(pathway_scores)[1:9]) {
    control <- c()
    for (c in clusters) {
      control_cluster <- sample(pathway_scores[!(pathway_scores$adtClusterID == c), p], length(pathway_scores[pathway_scores$adtClusterID == c, p]), replace = TRUE)
      control <- rbind(control, cbind(control_cluster, rep(c, times = length(control_cluster))))
    }
    control <- cbind(control, c(rep("control", times = nrow(control))))
    score_type <- strsplit(p, "_")
    colnames(control) <- c(paste0(score_type[[1]][1], "_score"), "adtClusterID", "Type")
    control <- as.data.frame(control)
    assign(paste0("control_", p), control)
    
    df_clusters <- subset(pathway_scores, select = c(p, "adtClusterID"))
    pvalue_cluster <- c()
    for (c in clusters) {
      cases <- df_clusters[df_clusters$adtClusterID == c,]
      cases$Type <- c(rep("case", times = nrow(cases)))
      control_df <- get(x = paste0("control_", p))
      colnames(cases) <- colnames(control_df)
      df <- rbind(cases, control_df[control_df$adtClusterID == c,])
      
      #The p-value is calculated via pairwise Wilcoxon rank sum test
      pvalue_cluster <- c(pvalue_cluster, pairwise.wilcox.test(as.numeric(df[,1]), df$Type, p.adjust.method = "BH")$p.value)
    }
    dotplot_df_clusters[dotplot_df_clusters$Pathway == p, "pvalue"] <- pvalue_cluster
  }
  
  #Nice names for the pathways
  #Function for fist letter of a string capitalized
  capFirst <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
  }
  dotplot_df_clusters$Pathway <- sapply(dotplot_df_clusters$Pathway, function(x) capFirst(strsplit(as.character(x), "_")[[1]][1]))
  
  plot <- ggplot(dotplot_df_clusters,
                 aes(
                   x = Pathway,
                   y = adtClusterID,
                   color = Median,
                   size = -log10(pvalue)
                 )) +
    geom_point() +
    scale_colour_gradient2(high = "red", mid = "lightgrey", low = "blue", limits = c(min(dotplot_df_clusters$Median), max(dotplot_df_clusters$Median)))+
    scale_size(range = c(3, 8)) +
    scale_y_continuous(breaks = as.numeric(levels(pathway_scores$adtClusterID))) +
    labs(x = "Pathway", y = "ADT Cluster", title = patient, color = "Median", size = "Significance") +
    theme_bw() +
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.6)
    )
  return(plot)
}
