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

#Calculation tolerance score for DK126
#Setting as directory the source file directory
setwd('..')

#Load the files that are needed
load(file = "SJG017_DK126_ADT.RData")
DK126 <- readRDS("DK126_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG017_DK126_GSEA_ToleranceScore"), showWarnings = FALSE)
setwd(file.path("SJG017_DK126_GSEA_ToleranceScore"))

std_DK126 <- DK126[["TMM"]]@scale.data
tolerance_score_DK126 <- std_DK126[ADT_tolerance_score[1],]
for (adt in ADT_tolerance_score[2:length(ADT_tolerance_score)]){
  tolerance_score_DK126 <- tolerance_score_DK126 + std_DK126[adt,]
}

tolerance_score_DK126 <- tolerance_score_DK126/length(ADT_tolerance_score)
tolerance_score_DK126 <- as.data.frame(tolerance_score_DK126)
tolerance_score_DK126$Cell <- rownames(tolerance_score_DK126)
tolerance_score_DK126$ADT_cluster <- as.numeric(lapply(tolerance_score_DK126$Cell, function(x) DK126[["adtClusterID"]][x,])) -1
tolerance_score_DK126$RNA_cluster <- as.numeric(lapply(tolerance_score_DK126$Cell, function(x) DK126[["rnaClusterID"]][x,])) -1
tolerance_score_DK126$Condition <- sapply(tolerance_score_DK126$Cell, function(x) DK126[["orig.ident"]][x,])
tolerance_score_DK126 <- tolerance_score_DK126[order(tolerance_score_DK126$tolerance_score_DK126, decreasing = TRUE),]

#Save tolerance score 
write.csv(tolerance_score_DK126, file = "SJG017_DK126_tolerance_score.csv")

#Ranked character vector of symbols for the KS test
DK126_ranked.signature <- tolerance_score_DK126$Cell

#ADT Clusters
#List of vectors corresponding to the cells sets
DK126_ADTClusters_GSEA_cellsets <- list("ADT Cluster 0" = tolerance_score_DK126[tolerance_score_DK126$ADT_cluster == 0,]$Cell,
                                        "ADT Cluster 1" = tolerance_score_DK126[tolerance_score_DK126$ADT_cluster == 1,]$Cell,
                                        "ADT Cluster 2" = tolerance_score_DK126[tolerance_score_DK126$ADT_cluster == 2,]$Cell,
                                        "ADT Cluster 3" = tolerance_score_DK126[tolerance_score_DK126$ADT_cluster == 3,]$Cell,
                                        "ADT Cluster 4" = tolerance_score_DK126[tolerance_score_DK126$ADT_cluster == 4,]$Cell)

#Hyper enrichment
DK126_ADTClusters_GSEA <- hypeR(DK126_ranked.signature, DK126_ADTClusters_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK126_ADTClusters_GSEA$plots[[1]], filename = "DK126_GSEA_ADTClusters_Plot1.pdf", dpi = 700)
ggsave(DK126_ADTClusters_GSEA$plots[[2]], filename = "DK126_GSEA_ADTClusters_Plot2.pdf", dpi = 700)
ggsave(DK126_ADTClusters_GSEA$plots[[3]], filename = "DK126_GSEA_ADTClusters_Plot3.pdf", dpi = 700)
ggsave(DK126_ADTClusters_GSEA$plots[[4]], filename = "DK126_GSEA_ADTClusters_Plot4.pdf", dpi = 700)
ggsave(DK126_ADTClusters_GSEA$plots[[5]], filename = "DK126_GSEA_ADTClusters_Plot5.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK126_ADTClusters_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK126_ADTClusters_GSEA), filename = "DK126_GSEA_ADTClusters_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK126_ADTClusters_GSEA)

write.csv(DK126_ADTClusters_GSEA$data, file = "DK126_GSEA_ADTClusters.csv")

#List of vectors corresponding to the cells sets
DK126_conditions_GSEA_cellsets <- list("Condition NTDMSO" = tolerance_score_DK126[tolerance_score_DK126$Condition == "NTDMSO",]$Cell,
                                       "Condition NTAG" = tolerance_score_DK126[tolerance_score_DK126$Condition == "NTAG",]$Cell,
                                       "Condition AGDMSO" = tolerance_score_DK126[tolerance_score_DK126$Condition == "AGDMSO",]$Cell,
                                       "Condition AGAG" = tolerance_score_DK126[tolerance_score_DK126$Condition == "AGAG",]$Cell)

#Hyper enrichment
DK126_conditions_GSEA <- hypeR(DK126_ranked.signature, DK126_conditions_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK126_conditions_GSEA$plots[[1]], filename = "DK126_GSEA_conditions_Plot1.pdf", dpi = 700)
ggsave(DK126_conditions_GSEA$plots[[2]], filename = "DK126_GSEA_conditions_Plot2.pdf", dpi = 700)
ggsave(DK126_conditions_GSEA$plots[[3]], filename = "DK126_GSEA_conditions_Plot3.pdf", dpi = 700)
ggsave(DK126_conditions_GSEA$plots[[4]], filename = "DK126_GSEA_conditions_Plot4.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK126_conditions_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK126_conditions_GSEA), filename = "DK126_GSEA_conditions_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK126_conditions_GSEA)

write.csv(DK126_conditions_GSEA$data, file = "DK126_GSEA_Conditions.csv")

#Calculation tolerance score for DK104
#Setting as directory the source file directory
setwd('..')

#Load the files that are needed
load(file = "SJG026_DK104_ADT.RData")
DK104 <- readRDS("DK104_Seurat_split.integrated.rds")

#Creating directory for output files
dir.create(file.path("SJG026_DK104_GSEA_ToleranceScore"), showWarnings = FALSE)
setwd(file.path("SJG026_DK104_GSEA_ToleranceScore"))

std_DK104 <- DK104[["TMM"]]@scale.data
tolerance_score_DK104 <- std_DK104[ADT_tolerance_score[1],]
for (adt in ADT_tolerance_score[2:length(ADT_tolerance_score)]){
  tolerance_score_DK104 <- tolerance_score_DK104 + std_DK104[adt,]
}

tolerance_score_DK104 <- tolerance_score_DK104/length(ADT_tolerance_score)
tolerance_score_DK104 <- as.data.frame(tolerance_score_DK104)
tolerance_score_DK104$Cell <- rownames(tolerance_score_DK104)
tolerance_score_DK104$ADT_cluster <- as.numeric(lapply(tolerance_score_DK104$Cell, function(x) DK104[["adtClusterID"]][x,])) -1
tolerance_score_DK104$RNA_cluster <- as.numeric(lapply(tolerance_score_DK104$Cell, function(x) DK104[["rnaClusterID"]][x,])) -1
tolerance_score_DK104$Condition <- sapply(tolerance_score_DK104$Cell, function(x) DK104[["orig.ident"]][x,])
tolerance_score_DK104 <- tolerance_score_DK104[order(tolerance_score_DK104$tolerance_score_DK104, decreasing = TRUE),]

#Save tolerance score to xlsx
write.csv(tolerance_score_DK104, file = "SJG026_DK104_tolerance_score.csv")

#Ranked character vector of symbols for the KS test
DK104_ranked.signature <- tolerance_score_DK104$Cell

#ADT Clusters
#List of vectors corresponding to the cells sets
DK104_ADTClusters_GSEA_cellsets <- list("ADT Cluster 0" = tolerance_score_DK104[tolerance_score_DK104$ADT_cluster == 0,]$Cell,
                                        "ADT Cluster 1" = tolerance_score_DK104[tolerance_score_DK104$ADT_cluster == 1,]$Cell,
                                        "ADT Cluster 2" = tolerance_score_DK104[tolerance_score_DK104$ADT_cluster == 2,]$Cell)

#Hyper enrichment
DK104_ADTClusters_GSEA <- hypeR(DK104_ranked.signature, DK104_ADTClusters_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK104_ADTClusters_GSEA$plots[[1]], filename = "DK104_GSEA_ADTClusters_Plot1.pdf", dpi = 700)
ggsave(DK104_ADTClusters_GSEA$plots[[2]], filename = "DK104_GSEA_ADTClusters_Plot2.pdf", dpi = 700)
ggsave(DK104_ADTClusters_GSEA$plots[[3]], filename = "DK104_GSEA_ADTClusters_Plot3.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK104_ADTClusters_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK104_ADTClusters_GSEA), filename = "DK104_GSEA_ADTClusters_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK104_ADTClusters_GSEA)

write.csv(DK104_ADTClusters_GSEA$data, file = "DK104_GSEA_ADTClusters.csv")

#List of vectors corresponding to the cells sets
DK104_conditions_GSEA_cellsets <- list("Condition NTDMSO" = tolerance_score_DK104[tolerance_score_DK104$Condition == "NTDMSO",]$Cell,
                                       "Condition NTAG" = tolerance_score_DK104[tolerance_score_DK104$Condition == "NTAG",]$Cell,
                                       "Condition AGDMSO" = tolerance_score_DK104[tolerance_score_DK104$Condition == "AGDMSO",]$Cell,
                                       "Condition AGAG" = tolerance_score_DK104[tolerance_score_DK104$Condition == "AGAG",]$Cell)

#Hyper enrichment
DK104_conditions_GSEA <- hypeR(DK104_ranked.signature, DK104_conditions_GSEA_cellsets, test = "kstest", plotting = TRUE)
ggsave(DK104_conditions_GSEA$plots[[1]], filename = "DK104_GSEA_conditions_Plot1.pdf", dpi = 700)
ggsave(DK104_conditions_GSEA$plots[[2]], filename = "DK104_GSEA_conditions_Plot2.pdf", dpi = 700)
ggsave(DK104_conditions_GSEA$plots[[3]], filename = "DK104_GSEA_conditions_Plot3.pdf", dpi = 700)
ggsave(DK104_conditions_GSEA$plots[[4]], filename = "DK104_GSEA_conditions_Plot4.pdf", dpi = 700)

#Show interactive table 
hyp_show(DK104_conditions_GSEA)

#Plots dot plots
ggsave(hyp_dots(DK104_conditions_GSEA), filename = "DK104_GSEA_conditions_dotplot.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(DK104_conditions_GSEA)

write.csv(DK104_conditions_GSEA$data, file = "DK104_GSEA_Conditions.csv")

#Preparation ridgeplot
#Setting as directory the source file directory
setwd('..')

tolerance_score_DK114$Dataset <- rep(x = "DK114", times = nrow(tolerance_score_DK114))
tolerance_score_DK104$Dataset <- rep(x = "DK104", times = nrow(tolerance_score_DK104))
tolerance_score_DK118$Dataset <- rep(x = "DK118", times = nrow(tolerance_score_DK118))
tolerance_score_DK126$Dataset <- rep(x = "DK126", times = nrow(tolerance_score_DK126))

tolerance_score_DK114$ADT_cluster <- as.character(tolerance_score_DK114$ADT_cluster)
tolerance_score_DK114$RNA_cluster <- as.character(tolerance_score_DK114$RNA_cluster)
tolerance_score_DK114$Condition <- as.factor(tolerance_score_DK114$Condition)
tolerance_score_DK114$Condition <- factor(tolerance_score_DK114$Condition, levels = c("NT", "AG", "R", "RAG"), ordered = TRUE)

tolerance_score_DK104$ADT_cluster <- as.character(tolerance_score_DK104$ADT_cluster)
tolerance_score_DK104$RNA_cluster <- as.character(tolerance_score_DK104$RNA_cluster)
tolerance_score_DK104$Condition <- as.factor(tolerance_score_DK104$Condition)
tolerance_score_DK104$Condition <- factor(tolerance_score_DK104$Condition, levels = c("NTDMSO", "NTAG", "AGDMSO", "AGAG"), ordered = TRUE)

tolerance_score_DK118$ADT_cluster <- as.character(tolerance_score_DK118$ADT_cluster)
tolerance_score_DK118$RNA_cluster <- as.character(tolerance_score_DK118$RNA_cluster)
tolerance_score_DK118$Condition <- as.factor(tolerance_score_DK118$Condition)
tolerance_score_DK118$Condition <- factor(tolerance_score_DK118$Condition, levels = c("NT", "AG", "R", "RAG"), ordered = TRUE)

tolerance_score_DK126$ADT_cluster <- as.character(tolerance_score_DK126$ADT_cluster)
tolerance_score_DK126$RNA_cluster <- as.character(tolerance_score_DK126$RNA_cluster)
tolerance_score_DK126$Condition <- as.factor(tolerance_score_DK126$Condition)
tolerance_score_DK126$Condition <- factor(tolerance_score_DK126$Condition, levels = c("NTDMSO", "NTAG", "AGDMSO", "AGAG"), ordered = TRUE)

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

#DK104
p1<- ggplot(tolerance_score_DK104, aes(x = tolerance_score_DK104, y = Condition, fill = Condition))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(tolerance_score_DK104, aes(x = tolerance_score_DK104, y = ADT_cluster, fill = ADT_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "ADT Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(tolerance_score_DK104, aes(x = tolerance_score_DK104, y = RNA_cluster, fill = RNA_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "RNA Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = wrap_plots(p1, p2, p3, ncol = 3) + plot_annotation(title = "DK104"), 
       filename = "DK104_tolerance_score.pdf", width = 15, heigh = 4, dpi = 700)

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

#DK126
p1<- ggplot(tolerance_score_DK126, aes(x = tolerance_score_DK126, y = Condition, fill = Condition))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(tolerance_score_DK126, aes(x = tolerance_score_DK126, y = ADT_cluster, fill = ADT_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "ADT Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(tolerance_score_DK126, aes(x = tolerance_score_DK126, y = RNA_cluster, fill = RNA_cluster))+
  geom_density_ridges_gradient(scale = 4, show.legend = FALSE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0.3, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "RNA Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = wrap_plots(p1, p2, p3, ncol = 3) + plot_annotation(title = "DK126"), 
       filename = "DK126_tolerance_score.pdf", width = 15, heigh = 4, dpi = 700)

#DK114+DK104

#Creating directory for output files
dir.create(file.path("SJG026_DK114+DK104_GSEA_ToleranceScore"), showWarnings = FALSE)
setwd(file.path("SJG026_DK114+DK104_GSEA_ToleranceScore"))

#Clusters
#Changing the elements name to distinguish the different clusters from the different datasets
names(DK104_ADTClusters_GSEA_cellsets) <- lapply(names(DK104_ADTClusters_GSEA_cellsets), function(x) paste0("DK104 ", x))
names(DK114_ADTClusters_GSEA_cellsets) <- lapply(names(DK114_ADTClusters_GSEA_cellsets), function(x) paste0("DK114 ", x))

#Merging the individual cellsets 
all_data_ADTclusters_GSEA_cellsets_SJG026 <- c(DK104_ADTClusters_GSEA_cellsets, DK114_ADTClusters_GSEA_cellsets)
reduced_DK104_tol_score <- subset(tolerance_score_DK104, select = c(1,2))
reduced_DK114_tol_score <- subset(tolerance_score_DK114, select = c(1,2))
names(reduced_DK104_tol_score) <- names(reduced_DK114_tol_score) <- c("tolerance_score", "Cell")
all_data_SJG026_ranked.signature <- rbind(reduced_DK114_tol_score, reduced_DK104_tol_score)

#Ordering the dataset according to the score, obtaining an intergrative ranked signature
all_data_SJG026_ranked.signature <- all_data_SJG026_ranked.signature[order(all_data_SJG026_ranked.signature$tolerance_score, decreasing = TRUE), ]
all_data_SJG026_ranked.signature <- all_data_SJG026_ranked.signature$Cell

#GSEA
#Hyper enrichment
all_data_ADTClusters_GSEA_SJG026 <- hypeR(all_data_SJG026_ranked.signature, all_data_ADTclusters_GSEA_cellsets_SJG026, test = "kstest", plotting = TRUE)
all_data_ADTClusters_GSEA_SJG026$plots

#Show interactive table 
hyp_show(all_data_ADTClusters_GSEA_SJG026)

#Plots dot plots
ggsave(hyp_dots(all_data_ADTClusters_GSEA_SJG026), filename = "all_data_GSEA_ADTClusters_dotplot_SJG026.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(all_data_ADTClusters_GSEA_SJG026)

write.csv(all_data_ADTClusters_GSEA_SJG026$data, file = "all_data_GSEA_ADTClusters_SJG026.csv")

#Conditions
#Changing the elements name to distinguish the different clusters from the different datasets
names(DK104_conditions_GSEA_cellsets) <- lapply(names(DK104_conditions_GSEA_cellsets), function(x) paste0("DK104 ", x))
names(DK114_conditions_GSEA_cellsets) <- lapply(names(DK114_conditions_GSEA_cellsets), function(x) paste0("DK114 ", x))

#Merging the individual cellsets 
all_data_conditions_GSEA_cellsets_SJG026 <- c(DK104_conditions_GSEA_cellsets, DK114_conditions_GSEA_cellsets)

#GSEA
#Hyper enrichment
all_data_conditions_GSEA_SJG026 <- hypeR(all_data_SJG026_ranked.signature, all_data_conditions_GSEA_cellsets_SJG026, test = "kstest", plotting = TRUE)
all_data_conditions_GSEA_SJG026$plots

#Show interactive table 
hyp_show(all_data_conditions_GSEA_SJG026)

#Plots dot plots
ggsave(hyp_dots(all_data_conditions_GSEA_SJG026), filename = "all_data_GSEA_conditions_dotplot_SJG026.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(all_data_conditions_GSEA_SJG026)

write.csv(all_data_conditions_GSEA_SJG026$data, file = "all_data_GSEA_conditions_SJG026.csv")

setwd("..")

#DK118+DK126
#Creating directory for output files
dir.create(file.path("SJG017_DK118+DK126_GSEA_ToleranceScore"), showWarnings = FALSE)
setwd(file.path("SJG017_DK118+DK126_GSEA_ToleranceScore"))

#Clusters
#Changing the elements name to distinguish the different clusters from the different datasets
names(DK118_ADTClusters_GSEA_cellsets) <- lapply(names(DK118_ADTClusters_GSEA_cellsets), function(x) paste0("DK118 ", x))
names(DK126_ADTClusters_GSEA_cellsets) <- lapply(names(DK126_ADTClusters_GSEA_cellsets), function(x) paste0("DK126 ", x))

#Merging the individual cellsets 
all_data_ADTclusters_GSEA_cellsets_SJG017 <- c(DK126_ADTClusters_GSEA_cellsets, DK118_ADTClusters_GSEA_cellsets)
reduced_DK126_tol_score <- subset(tolerance_score_DK126, select = c(1,2))
reduced_DK118_tol_score <- subset(tolerance_score_DK118, select = c(1,2))
names(reduced_DK126_tol_score) <- names(reduced_DK118_tol_score) <- c("tolerance_score", "Cell")
all_data_SJG017_ranked.signature <- rbind(reduced_DK118_tol_score, reduced_DK126_tol_score)

#Ordering the dataset according to the score, obtaining an intergrative ranked signature
all_data_SJG017_ranked.signature <- all_data_SJG017_ranked.signature[order(all_data_SJG017_ranked.signature$tolerance_score, decreasing = TRUE), ]
all_data_SJG017_ranked.signature <- all_data_SJG017_ranked.signature$Cell

#GSEA
#Hyper enrichment
all_data_ADTClusters_GSEA_SJG017 <- hypeR(all_data_SJG017_ranked.signature, all_data_ADTclusters_GSEA_cellsets_SJG017, test = "kstest", plotting = TRUE)
all_data_ADTClusters_GSEA_SJG017$plots

#Show interactive table 
hyp_show(all_data_ADTClusters_GSEA_SJG017)

#Plots dot plots
ggsave(hyp_dots(all_data_ADTClusters_GSEA_SJG017), filename = "all_data_GSEA_ADTClusters_dotplot_SJG017.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(all_data_ADTClusters_GSEA_SJG017)

write.csv(all_data_ADTClusters_GSEA_SJG017$data, file = "all_data_GSEA_ADTClusters_SJG017.csv")

#Conditions
#Changing the elements name to distinguish the different clusters from the different datasets
names(DK118_conditions_GSEA_cellsets) <- lapply(names(DK118_conditions_GSEA_cellsets), function(x) paste0("DK118 ", x))
names(DK126_conditions_GSEA_cellsets) <- lapply(names(DK126_conditions_GSEA_cellsets), function(x) paste0("DK126 ", x))

#Merging the individual cellsets 
all_data_conditions_GSEA_cellsets_SJG017 <- c(DK126_conditions_GSEA_cellsets, DK118_conditions_GSEA_cellsets)

#GSEA
#Hyper enrichment
all_data_conditions_GSEA_SJG017 <- hypeR(all_data_SJG017_ranked.signature, all_data_conditions_GSEA_cellsets_SJG017, test = "kstest", plotting = TRUE)
all_data_conditions_GSEA_SJG017$plots

#Show interactive table 
hyp_show(all_data_conditions_GSEA_SJG017)

#Plots dot plots
ggsave(hyp_dots(all_data_conditions_GSEA_SJG017), filename = "all_data_GSEA_conditions_dotplot_SJG017.pdf", dpi = 700, width = 5, height = 5)

#Plot enrichment map
hyp_emap(all_data_conditions_GSEA_SJG017)

write.csv(all_data_conditions_GSEA_SJG017$data, file = "all_data_GSEA_conditions_SJG017.csv")

#Change directory
setwd("..")

#Ridge plots for the difference in the conditions' tolerance score
colnames(tolerance_score_DK104) <- c("tolerance_score", colnames(tolerance_score_DK104)[2:ncol(tolerance_score_DK104)])
colnames(tolerance_score_DK126) <- c("tolerance_score", colnames(tolerance_score_DK126)[2:ncol(tolerance_score_DK126)])
colnames(tolerance_score_DK114) <- c("tolerance_score", colnames(tolerance_score_DK114)[2:ncol(tolerance_score_DK114)])
colnames(tolerance_score_DK118) <- c("tolerance_score", colnames(tolerance_score_DK118)[2:ncol(tolerance_score_DK118)])

tolerance_score_DK104$Patient <- rep("SJG-026", times = nrow(tolerance_score_DK104))
tolerance_score_DK114$Patient <- rep("SJG-026", times = nrow(tolerance_score_DK114))

tolerance_score_DK118$Patient <- rep("SJG-017", times = nrow(tolerance_score_DK118))
tolerance_score_DK126$Patient <- rep("SJG-017", times = nrow(tolerance_score_DK126))

short_treatment_tolerance <- rbind(tolerance_score_DK104, tolerance_score_DK126)
long_treatment_tolerance <- rbind(tolerance_score_DK114, tolerance_score_DK118)
all_datasets_tolerance_score <- rbind(long_treatment_tolerance, short_treatment_tolerance)

#Saving ridge plot

svg("RidgePlot_ToleranceScore_Conditions.svg")
ggplot(all_datasets_tolerance_score, aes(x = tolerance_score, y = Dataset, fill = Condition))+
  geom_density_ridges_gradient(scale = 1.5, show.legend = TRUE) + theme_ridges() +
  scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.7) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  labs(x = "Tolerance Score", y = "Condition") +
  theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(Patient~.)
dev.off()

#Dotplot conditions tolerance scores
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

#DK126
#Dotplot conditions tolerance scores
dotplot_df_conditions <- data.frame(Condition = c("NTDMSO", "NTAG", "AGDMSO", "AGAG"))
dotplot_df_conditions$Median <- rep(NA, times = 4)
dotplot_df_conditions$Mean <- rep(NA, times = 4)
dotplot_df_conditions$pvalue <- rep(NA, times = 4)
for (row in 1:nrow(dotplot_df_conditions)){
  dotplot_df_conditions[row, "Median"] <- median(tolerance_score_DK126[tolerance_score_DK126$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  dotplot_df_conditions[row, "Mean"] <- mean(tolerance_score_DK126[tolerance_score_DK126$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  #Kolmogorov-Smirnov test
  dotplot_df_conditions[row, "pvalue"] <- format.pval(stats::ks.test(tolerance_score_DK126[tolerance_score_DK126$Condition == "NTDMSO",]$tolerance_score,
                                                                     tolerance_score_DK126[tolerance_score_DK126$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score,
                                                                     alternative = "two.sided", exact = TRUE)$p.value)
  if (dotplot_df_conditions[row, "pvalue"] == "< 2.22e-16"){
    dotplot_df_conditions[row, "pvalue"] <- 2.22e-16
  }
}

#Adjustment p-value with BH
dotplot_df_conditions$adj_pval <- sapply(dotplot_df_conditions$pvalue, function(x) p.adjust(x, method = "BH"))
dotplot_df_conditions$log10_adj_pval <- sapply(dotplot_df_conditions$adj_pval, function(x) -log10(x))

#Dataset indication
dotplot_df_conditions$dataset <- rep("DK126", times = 4)
dotplot_df_conditions$patient <- rep("SJG-017", times = 4)
dotplot_df_conditions$Selection <- rep("Short treatment", times = 4)
dotplot_df_conditions_DK126 <- dotplot_df_conditions

#DK104
#Dotplot conditions tolerance scores
dotplot_df_conditions <- data.frame(Condition = c("NTDMSO", "NTAG", "AGDMSO", "AGAG"))
dotplot_df_conditions$Median <- rep(NA, times = 4)
dotplot_df_conditions$Mean <- rep(NA, times = 4)
dotplot_df_conditions$pvalue <- rep(NA, times = 4)
for (row in 1:nrow(dotplot_df_conditions)){
  dotplot_df_conditions[row, "Median"] <- median(tolerance_score_DK104[tolerance_score_DK104$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  dotplot_df_conditions[row, "Mean"] <- mean(tolerance_score_DK104[tolerance_score_DK104$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score)
  #Kolmogorov-Smirnov test
  dotplot_df_conditions[row, "pvalue"] <- format.pval(stats::ks.test(tolerance_score_DK104[tolerance_score_DK104$Condition == "NTDMSO",]$tolerance_score,
                                                                     tolerance_score_DK104[tolerance_score_DK104$Condition == dotplot_df_conditions[row, "Condition"],]$tolerance_score,
                                                                     alternative = "two.sided", exact = TRUE)$p.value)
  if (dotplot_df_conditions[row, "pvalue"] == "< 2.22e-16"){
    dotplot_df_conditions[row, "pvalue"] <- 2.22e-16
  }
}

#Adjustment p-value with BH
dotplot_df_conditions$adj_pval <- sapply(dotplot_df_conditions$pvalue, function(x) p.adjust(x, method = "BH"))
dotplot_df_conditions$log10_adj_pval <- sapply(dotplot_df_conditions$adj_pval, function(x) -log10(x))

#Dataset indication
dotplot_df_conditions$dataset <- rep("DK104", times = 4)
dotplot_df_conditions$patient <- rep("SJG-026", times = 4)
dotplot_df_conditions$Selection <- rep("Short treatment", times = 4)
dotplot_df_conditions_DK104 <- dotplot_df_conditions

dotplot_df_conditions_longtreatment <- rbind(dotplot_df_conditions_DK114,
                                             dotplot_df_conditions_DK118)

dotplot_df_conditions_shorttreatment <- rbind(dotplot_df_conditions_DK104,
                                              dotplot_df_conditions_DK126)

p1 <- ggplot(dotplot_df_conditions_shorttreatment,
             aes(
               x = Condition,
               y = Selection,
               color = Mean,
               size = log10_adj_pval
             )) +
  geom_point() +
  scale_size(breaks = waiver(), limits = c(0,17), range = c(3,10), labels = c(0, 4, 8, 12, "> 16")) + 
  scale_colour_gradient2(high = "red", mid = "lightgrey", low = "blue", midpoint = 0) +
  labs(x = "", y = "", color = "Mean", size = "-log10(adj. p-value)") +
  theme_bw() +
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.6),
  ) + facet_grid(.~patient)

levels_order <- c("NT", "AG", "R", "RAG")
dotplot_df_conditions_longtreatment$patient <- factor(dotplot_df_conditions_longtreatment$patient, levels = c("SJG-026", "SJG-017"))
  
p2 <- ggplot(dotplot_df_conditions_longtreatment,
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

final_plot <- ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right", align = "hv")

#Saving KS dot plot
svg(filename = "Tolerance_score_dotplot_KS_conditions.svg", height = 4)
annotate_figure(final_plot)
dev.off()

svg(filename = "DK114_DK118_only_Tolerance_score_dotplot_KS_conditions.svg", height = 2.5, width = 7)
p2
dev.off()
