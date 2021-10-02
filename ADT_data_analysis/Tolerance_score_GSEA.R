#########################################
#          ADT data analysis            #
#           Tolerance score             #
#########################################

library(hypeR)
library(ggplot2)
library(ggridges)
library(viridis)
library(patchwork)
library(ggpubr)

#Calculation tolerance score for DK114
#Loading ADT results
load(file = "SJG026_DK114_ADT.RData")
DK114 <- readRDS("DK114_Seurat_split.integrated.RDS")

#The tolerance score is defined according to the overlapping markers of ADT Cluster 2 in SJG026-DK114
ADT_tolerance_score <- rownames(ADT_cluster_2_specific_markers)
std_DK114 <- DK114[["TMM"]]@scale.data

#Creating directory for output files
dir.create(file.path("SJG026_DK114_GSEA_ToleranceScore"), showWarnings = FALSE)
setwd(file.path("SJG026_DK114_GSEA_ToleranceScore"))

#The tolerance score is computed according to the aforementioned definition
tolerance_score_DK114 <- std_DK114[ADT_tolerance_score[1],]
for (adt in ADT_tolerance_score[2:length(ADT_tolerance_score)]){
  tolerance_score_DK114 <- tolerance_score_DK114 + std_DK114[adt,]
}

tolerance_score_DK114 <- tolerance_score_DK114/length(ADT_tolerance_score)
tolerance_score_DK114 <- as.data.frame(tolerance_score_DK114)
tolerance_score_DK114$Cell <- rownames(tolerance_score_DK114)

#as.numeric is needed to avoid conflicts when saving the tolerance score data frame to csv
tolerance_score_DK114$ADT_cluster <- as.numeric(lapply(tolerance_score_DK114$Cell, function(x) DK114[["adtClusterID"]][x,])) - 1 
tolerance_score_DK114$RNA_cluster <- as.numeric(lapply(tolerance_score_DK114$Cell, function(x) DK114[["rnaClusterID"]][x,])) - 1
tolerance_score_DK114$Condition <- sapply(tolerance_score_DK114$Cell, function(x) DK114[["orig.ident"]][x,])
tolerance_score_DK114 <- tolerance_score_DK114[order(tolerance_score_DK114$tolerance_score_DK114, decreasing = TRUE),]

#Save tolerance score to xlsx
write.csv(tolerance_score_DK114, file = "SJG026_DK114_tolerance_score.csv")

#Ranked character vector of symbols for the KS test
DK114_ranked.signature <- tolerance_score_DK114$Cell

#ADT Clusters

#List of vectors corresponding to the cells sets
DK114_ADTClusters_GSEA_cellsets <- list("ADT Cluster 0" = tolerance_score_DK114[tolerance_score_DK114$ADT_cluster == 0,]$Cell,
                 "ADT Cluster 1" = tolerance_score_DK114[tolerance_score_DK114$ADT_cluster == 1,]$Cell,
                 "ADT Cluster 2" = tolerance_score_DK114[tolerance_score_DK114$ADT_cluster == 2,]$Cell,
                 "ADT Cluster 3" = tolerance_score_DK114[tolerance_score_DK114$ADT_cluster == 3,]$Cell,
                 "ADT Cluster 4" = tolerance_score_DK114[tolerance_score_DK114$ADT_cluster == 4,]$Cell)

#Hyper enrichment
DK114_ADTClusters_GSEA <- hypeR(DK114_ranked.signature, DK114_ADTClusters_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK114_ADTClusters_GSEA$plots[[1]], filename = "DK114_GSEA_ADTClusters_Plot1.pdf", dpi = 700)
ggsave(DK114_ADTClusters_GSEA$plots[[2]], filename = "DK114_GSEA_ADTClusters_Plot2.pdf", dpi = 700)
ggsave(DK114_ADTClusters_GSEA$plots[[3]], filename = "DK114_GSEA_ADTClusters_Plot3.pdf", dpi = 700)
ggsave(DK114_ADTClusters_GSEA$plots[[4]], filename = "DK114_GSEA_ADTClusters_Plot4.pdf", dpi = 700)
ggsave(DK114_ADTClusters_GSEA$plots[[5]], filename = "DK114_GSEA_ADTClusters_Plot5.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK114_ADTClusters_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK114_ADTClusters_GSEA), filename = "DK114_GSEA_ADTClusters_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK114_ADTClusters_GSEA)

write.csv(DK114_ADTClusters_GSEA$data, file = "DK114_GSEA_ADTClusters.csv")

#List of vectors corresponding to the cells sets
DK114_conditions_GSEA_cellsets <- list("Condition NT" = tolerance_score_DK114[tolerance_score_DK114$Condition == "NT",]$Cell,
                                  "Condition AG" = tolerance_score_DK114[tolerance_score_DK114$Condition == "AG",]$Cell,
                                  "Condition R" = tolerance_score_DK114[tolerance_score_DK114$Condition == "R",]$Cell,
                                  "Condition RAG" = tolerance_score_DK114[tolerance_score_DK114$Condition == "RAG",]$Cell)

#Hyper enrichment
DK114_conditions_GSEA <- hypeR(DK114_ranked.signature, DK114_conditions_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK114_conditions_GSEA$plots[[1]], filename = "DK114_GSEA_conditions_Plot1.pdf", dpi = 700)
ggsave(DK114_conditions_GSEA$plots[[2]], filename = "DK114_GSEA_conditions_Plot2.pdf", dpi = 700)
ggsave(DK114_conditions_GSEA$plots[[3]], filename = "DK114_GSEA_conditions_Plot3.pdf", dpi = 700)
ggsave(DK114_conditions_GSEA$plots[[4]], filename = "DK114_GSEA_conditions_Plot4.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK114_conditions_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK114_conditions_GSEA), filename = "DK114_GSEA_conditions_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK114_conditions_GSEA)

write.csv(DK114_conditions_GSEA$data, file = "DK114_GSEA_Conditions.csv")

#Calculation tolerance score for DK118
#Setting as directory the source file directory
setwd('..')

#Load the files that are needed
load(file = "SJG017_DK118_ADT.RData")
DK118 <- readRDS("DK118_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG017_DK118_GSEA_ToleranceScore"), showWarnings = FALSE)
setwd(file.path("SJG017_DK118_GSEA_ToleranceScore"))

std_DK118 <- DK118[["TMM"]]@scale.data
tolerance_score_DK118 <- std_DK118[ADT_tolerance_score[1],]
for (adt in ADT_tolerance_score[2:length(ADT_tolerance_score)]){
  tolerance_score_DK118 <- tolerance_score_DK118 + std_DK118[adt,]
}

tolerance_score_DK118 <- tolerance_score_DK118/length(ADT_tolerance_score)
tolerance_score_DK118 <- as.data.frame(tolerance_score_DK118)
tolerance_score_DK118$Cell <- rownames(tolerance_score_DK118)
tolerance_score_DK118$ADT_cluster <- as.numeric(lapply(tolerance_score_DK118$Cell, function(x) DK118[["adtClusterID"]][x,])) -1
tolerance_score_DK118$RNA_cluster <- as.numeric(lapply(tolerance_score_DK118$Cell, function(x) DK118[["rnaClusterID"]][x,])) -1
tolerance_score_DK118$Condition <- sapply(tolerance_score_DK118$Cell, function(x) DK118[["orig.ident"]][x,])
tolerance_score_DK118 <- tolerance_score_DK118[order(tolerance_score_DK118$tolerance_score_DK118, decreasing = TRUE),]

#Save tolerance score to xlsx
write.csv(tolerance_score_DK118, file = "SJG017_DK118_tolerance_score.csv")

#Ranked character vector of symbols for the KS test
DK118_ranked.signature <- tolerance_score_DK118$Cell

#ADT Clusters

#List of vectors corresponding to the cells sets
DK118_ADTClusters_GSEA_cellsets <- list("ADT Cluster 0" = tolerance_score_DK118[tolerance_score_DK118$ADT_cluster == 0,]$Cell,
                                  "ADT Cluster 1" = tolerance_score_DK118[tolerance_score_DK118$ADT_cluster == 1,]$Cell,
                                  "ADT Cluster 2" = tolerance_score_DK118[tolerance_score_DK118$ADT_cluster == 2,]$Cell,
                                  "ADT Cluster 3" = tolerance_score_DK118[tolerance_score_DK118$ADT_cluster == 3,]$Cell)

#Hyper enrichment
DK118_ADTClusters_GSEA <- hypeR(DK118_ranked.signature, DK118_ADTClusters_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK118_ADTClusters_GSEA$plots[[1]], filename = "DK118_GSEA_ADTClusters_Plot1.pdf", dpi = 700)
ggsave(DK118_ADTClusters_GSEA$plots[[2]], filename = "DK118_GSEA_ADTClusters_Plot2.pdf", dpi = 700)
ggsave(DK118_ADTClusters_GSEA$plots[[3]], filename = "DK118_GSEA_ADTClusters_Plot3.pdf", dpi = 700)
ggsave(DK118_ADTClusters_GSEA$plots[[4]], filename = "DK118_GSEA_ADTClusters_Plot4.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK118_ADTClusters_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK118_ADTClusters_GSEA), filename = "DK118_GSEA_ADTClusters_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK118_ADTClusters_GSEA)

write.csv(DK118_ADTClusters_GSEA$data, file = "DK118_GSEA_ADTClusters.csv")

#List of vectors corresponding to the cells sets
DK118_conditions_GSEA_cellsets <- list("Condition NT" = tolerance_score_DK118[tolerance_score_DK118$Condition == "NT",]$Cell,
                                 "Condition AG" = tolerance_score_DK118[tolerance_score_DK118$Condition == "AG",]$Cell,
                                 "Condition R" = tolerance_score_DK118[tolerance_score_DK118$Condition == "R",]$Cell,
                                 "Condition RAG" = tolerance_score_DK118[tolerance_score_DK118$Condition == "RAG",]$Cell)

#Hyper enrichment
DK118_conditions_GSEA <- hypeR(DK118_ranked.signature, DK118_conditions_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK118_conditions_GSEA$plots[[1]], filename = "DK118_GSEA_conditions_Plot1.pdf", dpi = 700)
ggsave(DK118_conditions_GSEA$plots[[2]], filename = "DK118_GSEA_conditions_Plot2.pdf", dpi = 700)
ggsave(DK118_conditions_GSEA$plots[[3]], filename = "DK118_GSEA_conditions_Plot3.pdf", dpi = 700)
ggsave(DK118_conditions_GSEA$plots[[4]], filename = "DK118_GSEA_conditions_Plot4.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK118_conditions_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK118_conditions_GSEA), filename = "DK118_GSEA_conditions_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK118_conditions_GSEA)

write.csv(DK118_conditions_GSEA$data, file = "DK118_GSEA_Conditions.csv")

#Preparation ridgeplot
#Setting as directory the source file directory
setwd('..')

tolerance_score_DK114$Dataset <- rep(x = "DK114", times = nrow(tolerance_score_DK114))
tolerance_score_DK118$Dataset <- rep(x = "DK118", times = nrow(tolerance_score_DK118))
                                          
tolerance_score_DK114$ADT_cluster <- as.character(tolerance_score_DK114$ADT_cluster)
tolerance_score_DK114$RNA_cluster <- as.character(tolerance_score_DK114$RNA_cluster)
tolerance_score_DK114$Condition <- as.factor(tolerance_score_DK114$Condition)
tolerance_score_DK114$Condition <- factor(tolerance_score_DK114$Condition, levels = c("NT", "AG", "R", "RAG"), ordered = TRUE)

tolerance_score_DK118$ADT_cluster <- as.character(tolerance_score_DK118$ADT_cluster)
tolerance_score_DK118$RNA_cluster <- as.character(tolerance_score_DK118$RNA_cluster)
tolerance_score_DK118$Condition <- as.factor(tolerance_score_DK118$Condition)
tolerance_score_DK118$Condition <- factor(tolerance_score_DK118$Condition, levels = c("NT", "AG", "R", "RAG"), ordered = TRUE)

#DK114
p1<- ggplot(tolerance_score_DK114, aes(x = tolerance_score_DK114, y = Condition, fill = Condition))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(tolerance_score_DK114, aes(x = tolerance_score_DK114, y = ADT_cluster, fill = ADT_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "ADT Cluster") +
  ylim(NA, 1) +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(tolerance_score_DK114, aes(x = tolerance_score_DK114, y = RNA_cluster, fill = RNA_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "RNA Cluster") +
  ylim(NA, 1) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = wrap_plots(p1, p2, p3, ncol = 3) + plot_annotation(title = "DK114"), 
       filename = "DK114_tolerance_score.pdf", width = 15, heigh = 4, dpi = 700)

#DK118
p1<- ggplot(tolerance_score_DK118, aes(x = tolerance_score_DK118, y = Condition, fill = Condition))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(tolerance_score_DK118, aes(x = tolerance_score_DK118, y = ADT_cluster, fill = ADT_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "ADT Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(tolerance_score_DK118, aes(x = tolerance_score_DK118, y = RNA_cluster, fill = RNA_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "RNA Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = wrap_plots(p1, p2, p3, ncol = 3) + plot_annotation(title = "DK118"), 
       filename = "DK118_tolerance_score.pdf", width = 15, heigh = 4, dpi = 700)

#Change directory
setwd("..")

#Dot plot conditions tolerance scores
#DK114
dotplot_df_conditions <- data.frame(Condition = c("NT", "AG", "R", "RAG"))
dotplot_df_conditions$Median <- rep(NA, times = 4)
dotplot_df_conditions$Mean <- rep(NA, times = 4)
dotplot_df_conditions$pvalue <- rep(NA, times = 4)
for (row in 1:nrow(dotplot_df_conditions)){
  dotplot_df_conditions[row, "Median"] <- median(tolerance_score_DK114[tolerance_score_DK114$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  dotplot_df_conditions[row, "Mean"] <- mean(tolerance_score_DK114[tolerance_score_DK114$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  #Kolmogorov-Smirnov test
  dotplot_df_conditions[row, "pvalue"] <- format.pval(stats::ks.test(tolerance_score_DK114[tolerance_score_DK114$Condition == "NT",]$tolerance_score,
                                                                     tolerance_score_DK114[tolerance_score_DK114$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score,
                                                                     alternative = "two.sided", exact = TRUE)$p.value)
  
  if (dotplot_df_conditions[row, "pvalue"] == "< 2.22e-16"){
    dotplot_df_conditions[row, "pvalue"] <- 2.22e-16
  }
}

#Adjustment p-value with BH
dotplot_df_conditions$adj_pval <- sapply(dotplot_df_conditions$pvalue, function(x) p.adjust(x, method = "BH"))
dotplot_df_conditions$log10_adj_pval <- sapply(dotplot_df_conditions$adj_pval, function(x) -log10(x))

#Dataset indication
dotplot_df_conditions$dataset <- rep("DK114", times = 4)
dotplot_df_conditions$patient <- rep("SJG-026", times = 4)
dotplot_df_conditions$Selection <- rep("Long treatment", times = 4)
dotplot_df_conditions_DK114 <- dotplot_df_conditions

#DK118
#Dotplot conditions tolerance scores
dotplot_df_conditions <- data.frame(Condition = c("NT", "AG", "R", "RAG"))
dotplot_df_conditions$Median <- rep(NA, times = 4)
dotplot_df_conditions$Mean <- rep(NA, times = 4)
dotplot_df_conditions$pvalue <- rep(NA, times = 4)
for (row in 1:nrow(dotplot_df_conditions)){
  dotplot_df_conditions[row, "Median"] <- median(tolerance_score_DK118[tolerance_score_DK118$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  dotplot_df_conditions[row, "Mean"] <- mean(tolerance_score_DK118[tolerance_score_DK118$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  #Kolmogorov-Smirnov test
  dotplot_df_conditions[row, "pvalue"] <- format.pval(stats::ks.test(tolerance_score_DK118[tolerance_score_DK118$Condition == "NT",]$tolerance_score,
                                                                     tolerance_score_DK118[tolerance_score_DK118$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score,
                                                                     alternative = "two.sided", exact = TRUE)$p.value)
  if (dotplot_df_conditions[row, "pvalue"] == "< 2.22e-16"){
    dotplot_df_conditions[row, "pvalue"] <- 2.22e-16
  }
}

#Adjustment p-value with BH
dotplot_df_conditions$adj_pval <- sapply(dotplot_df_conditions$pvalue, function(x) p.adjust(x, method = "BH"))
dotplot_df_conditions$log10_adj_pval <- sapply(dotplot_df_conditions$adj_pval, function(x) -log10(x))

#Dataset indication
dotplot_df_conditions$dataset <- rep("DK118", times = 4)
dotplot_df_conditions$patient <- rep("SJG-017", times = 4)
dotplot_df_conditions$Selection <- rep("Long treatment", times = 4)
dotplot_df_conditions_DK118 <- dotplot_df_conditions

dotplot_df_conditions_longtreatment <- rbind(dotplot_df_conditions_DK114,
                                             dotplot_df_conditions_DK118)

levels_order <- c("NT", "AG", "R", "RAG")
dotplot_df_conditions_longtreatment$patient <- factor(dotplot_df_conditions_longtreatment$patient, levels = c("SJG-026", "SJG-017"))
  
p <- ggplot(dotplot_df_conditions_longtreatment,
             aes(
               x = factor(Condition, level = levels_order),
               y = Selection,
               color = Mean,
               size = log10_adj_pval
             )) +
  geom_point() +
  scale_size(breaks = waiver(), limits = c(0,17), range = c(3,10), labels = c(0, 4, 8, 12, "> 16")) + 
  scale_colour_gradient2(high = "red", mid = "lightgrey", low = "blue", midpoint = 0) +
  labs(x = "", y = NULL, color = "Mean", size = "-log10(adj. p-value)") +
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "white"),
    axis.text.x = element_text(angle = 90, vjust = 0.6),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
  ) + facet_grid(.~patient)

#Saving dot plot
svg(filename = "DK114_DK118_Tolerance_score_dotplot_KS_conditions.svg", height = 2.5, width = 7)
p
dev.off()
