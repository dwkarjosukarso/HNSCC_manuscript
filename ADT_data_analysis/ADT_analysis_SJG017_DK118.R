#########################################
#          ADT data analysis            #
#           Markers & Plots             #
#            SJG017 (DK118)             #
#               dataset                 #
#########################################

#This source contains the functions used for the Kolmogorov-Smirnov-based 
#differential expression analysis of ADT data.
source("ADT_KS-DE_function.R")

#This source contains the function needed for the feature+violin plots.
source("DFV_plot.R")

#This source contains the function needed for the computation of the 
#single cell pathway scores. 
source("ADT_single_cell_pathway_score.R")

#This source contains the function needed for the dotplots of the 
#single cell pathway scores. 
source("ADT_single_cell_pathway_scores_dotplot.R")

#ADT analysis of DK118 (SJG017)
DK118 <- readRDS("DK118_Seurat_split.integrated.RDS")
ADT_markers_KS(DK118)
DefaultAssay(DK118) <- "TMM"

#Clusters - Visualizations
Idents(DK118) <- "adtClusterID"

#Creating directory for output files
dir.create(file.path("SJG017_DK118"), showWarnings = FALSE)
setwd(file.path("SJG017_DK118"))

#Heatmaps & Feature plot + Violin plot - Clusters - top 10
for(object in as.vector(grep("ADT_markers_cluster_",names(.GlobalEnv),value=TRUE))){
  markers <- eval(as.symbol(object))
  ggsave(DoHeatmap(DK118, slot = "scale.data", assay = "TMM", features = rownames(head(markers[order(-markers$avg_diff),], n = 10))) + labs(title = paste0("SJG017_DK118 top 10 ",object)),
         filename = paste0(object, "_heatmap_top10_SJG017_DK118.pdf"), dpi = 700)
  ggsave(DFV_plot(seurat_object = DK118, 
                  slot = "scale.data", 
                  features = rownames(head(markers[order(-markers$avg_diff),], n = 10)), 
                  title = paste0("SJG017_DK118 top 10 ",object), subtitle = "ADT cluster-based analysis",
                  reduction_1 = "umap", reduction_2 = "umap_tmm"),
         limitsize = FALSE, width = 25, height = 53.6, dpi = 700, #5.36*10
         filename = paste0(object, "_Feat&Vln_top10_SJG017_DK118.pdf"))
}

ggsave(DoHeatmap(DK118, assay = "TMM", slot = "scale.data", 
                 features = c(rownames(ADT_markers_cluster_0[order(ADT_markers_cluster_0$comparison),]), 
                              rownames(ADT_markers_cluster_1[order(ADT_markers_cluster_1$comparison),]),
                              rownames(ADT_markers_cluster_2[order(ADT_markers_cluster_2$comparison),]),
                              rownames(ADT_markers_cluster_3[order(ADT_markers_cluster_3$comparison),])), 
                 group.bar = TRUE, 
                 group.by = "adtClusterID"), filename = "ADT_heatmap_markers_clusters_SJG017_DK118.pdf", dpi = 700, height = 15)


#Cluster-specific markers analysis

ggsave(DoHeatmap(DK118, assay = "TMM", slot = "scale.data", 
                 features = c(rownames(ADT_cluster_0_specific_markers[order(ADT_cluster_0_specific_markers$`0vs1_avg_diff`),]),
                              rownames(ADT_cluster_1_specific_markers[order(ADT_cluster_1_specific_markers$`1vs0_avg_diff`),]), 
                              rownames(ADT_cluster_2_specific_markers[order(ADT_cluster_2_specific_markers$`2vs0_avg_diff`),]), 
                              rownames(ADT_cluster_3_specific_markers[order(ADT_cluster_3_specific_markers$`3vs0_avg_diff`),])), 
                 group.bar = TRUE, 
                 group.by = "adtClusterID"), filename = "ADT_heatmap_cluster-specific_markers_SJG017_DK118.pdf", dpi = 700)

ggsave(DFV_plot(seurat_object = DK118, 
                slot = "scale.data", 
                features = c(rownames(ADT_cluster_0_specific_markers), rownames(ADT_cluster_1_specific_markers), rownames(ADT_cluster_2_specific_markers), rownames(ADT_cluster_3_specific_markers)), 
                title = "SJG017_DK118 Kolmogorov-Smirnov ADT cluster-specific markers among pairwise comparisons",
                reduction_1 = "umap", reduction_2 = "umap_tmm"),
       limitsize = FALSE, width = 25, height = 5.36*length(c(rownames(ADT_cluster_0_specific_markers), rownames(ADT_cluster_1_specific_markers), rownames(ADT_cluster_2_specific_markers), rownames(ADT_cluster_3_specific_markers))), 
       filename = "ADT_cluster-specific_markers_Feature&Vln_SJG017_DK118.pdf")


#Conditions - Visualizations
Idents(DK118) <- "orig.ident"
#Conditions re-assignment to have the desired order of the groups
levels(DK118) <- c("NT", "AG", "R", "RAG")

#Heatmaps & Feature plot + Violin plot- Conditions - top 10
for(object in as.vector(grep("ADT_markers_condition_",names(.GlobalEnv),value=TRUE))){
  markers <- eval(as.symbol(object))
  ggsave(DoHeatmap(DK118, features = rownames(head(markers[order(-markers$avg_diff),], n = 10))) + labs(title = paste0("SJG017_DK118 ",object)),
         filename = paste0(object, "_heatmap_top10_SJG017_DK118.pdf"), dpi = 700)
  ggsave(DFV_plot(seurat_object = DK118, 
                  slot = "scale.data", 
                  features = rownames(head(markers[order(-markers$avg_diff),], n = 10)), 
                  title = paste0("SJG017_DK118 top 10 ",object), subtitle = "ADT condition-based analysis",
                  reduction_1 = "umap", reduction_2 = "umap_tmm"),
         limitsize = FALSE, width = 25, height = 53.6, dpi = 700, #5.36*10
         filename = paste0(object, "_Feat&Vln_top10_SJG017_DK118.pdf"))
}

ggsave(DoHeatmap(DK118, assay = "TMM", slot = "scale.data", 
                 features = c(rownames(ADT_markers_condition_NT[order(ADT_markers_condition_NT$comparison),]), 
                              rownames(ADT_markers_condition_AG[order(ADT_markers_condition_AG$comparison),]),
                              rownames(ADT_markers_condition_R[order(ADT_markers_condition_R$comparison),]),
                              rownames(ADT_markers_condition_RAG[order(ADT_markers_condition_RAG$comparison),])), 
                 group.bar = TRUE) + labs(title = "SJG017_DK118 differentially expressed\nADT per each condition", subtitle = "Kolmogorov-Smirnov DE, markers ordered by comparison"),
       filename = "ADT_heatmap_markers_conditions_SJG017_DK118.pdf", dpi = 700, height = 15)


#Condition-specific markers analysis

ggsave(DoHeatmap(DK118, assay = "TMM", slot = "scale.data", 
                 features = c(rownames(ADT_condition_NT_specific_markers), 
                              rownames(ADT_condition_AG_specific_markers),
                              rownames(ADT_condition_R_specific_markers),
                              rownames(ADT_condition_RAG_specific_markers)), 
                 group.bar = TRUE), filename = "ADT_heatmap_condition-specific_markers_SJG017_DK118.pdf", dpi = 700)

ggsave(DFV_plot(seurat_object = DK118, 
                slot = "scale.data", 
                features = c(rownames(ADT_condition_NT_specific_markers), 
                             rownames(ADT_condition_AG_specific_markers),
                             rownames(ADT_condition_R_specific_markers),
                             rownames(ADT_condition_RAG_specific_markers)), 
                title = "SJG017_DK118 Kolmogorov-Smirnov ADT condition-specific markers among pairwise comparisons",
                reduction_1 = "umap", reduction_2 = "umap_tmm"),
       limitsize = FALSE, width = 25, height = 5.36*length(c(rownames(ADT_condition_NT_specific_markers), 
                                                             rownames(ADT_condition_AG_specific_markers),
                                                             rownames(ADT_condition_R_specific_markers),
                                                             rownames(ADT_condition_RAG_specific_markers))), 
       filename = "ADT_condition-specific_markers_Feature&Vln_SJG017_DK118.pdf")


#Saving csv data

for (file in as.vector(grep("ADT_markers_cluster",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG017_DK118.csv"))
}

for (file in as.vector(grep("ADT_markers_condition",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG017_DK118.csv"))
}

for (file in as.vector(grep("ADT_cluster_",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG017_DK118.csv"))
}

for (file in as.vector(grep("ADT_condition_",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG017_DK118.csv"))
}

#Single-cell pathway scores
SJG017_DK118_pathway_scores <- sc_pathway_scores(DK118)

#Saving csv file
write.csv(SJG017_DK118_pathway_scores, file = "ADT_pathway_scores_SJG017_DK118.csv")

#Packages for plots
library(reshape2)
library(viridis)
library(patchwork)
library(ggplot2)

plot_SJG017_DK118_pathway_scores <- melt(SJG017_DK118_pathway_scores)
colnames(plot_SJG017_DK118_pathway_scores) <- c("Conditions", "ADT Cluster", "Sample", "Pathway", "Score")

#Violin plot
ggsave(ggplot(plot_SJG017_DK118_pathway_scores, aes(Pathway, Score, fill=Conditions)) + 
         geom_violin(position = position_dodge(width = 0.9))+
         labs(title="SJG017_DK118 ADT pathway scores", y = "Score")+
         geom_violin(trim=FALSE)+
         scale_fill_brewer(palette="RdBu") + theme_bw() +
         theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank()), width = 15, filename = "ADT_pathway_scores_violinplot_SJG017_DK118.pdf", limitsize = FALSE)

ggsave(ggplot(plot_SJG017_DK118_pathway_scores, aes(Pathway, Score, fill = `ADT Cluster`)) + 
         geom_violin(position = position_dodge(width = 0.9))+
         labs(title="SJG017_DK118 ADT pathway scores", y = "Score")+
         geom_violin(trim=FALSE)+
         scale_fill_brewer(palette="RdBu") + theme_bw() +
         theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank()), width = 15, filename = "ADT_pathway_scores_violinplot_clusters_SJG017_DK118.pdf", limitsize = FALSE)

#Dot plots
dotplots_clusters <- dotplot_pathway_scores_clusters(pathway_scores = SJG017_DK118_pathway_scores, clusters = c(0, 1, 2, 3), patient = "SJG017_DK118")
ggsave(filename = "ADT_Dotplot_pathwayscores_clusters_SJG017_DK118.pdf", plot = dotplots_clusters, dpi = 700, width = 4.4)

dotplots_conditions <- dotplot_pathway_scores_conditions(pathway_scores = SJG017_DK118_pathway_scores, patient = "SJG017_DK118")
ggsave(filename = "ADT_Dotplot_pathwayscores_conditions_SJG017_DK118.pdf", plot = dotplots_conditions, dpi = 700, width = 4.4)

save.image("SJG017_DK118_ADT.RData")

#Adding metadata referring to the pathway scores
for (pathway in colnames(SJG017_DK118_pathway_scores)[1:9]){
  DK118 <- AddMetaData(DK118, metadata = SJG017_DK118_pathway_scores[, pathway], col.name = pathway)
}

#Feature plot
plot_data = function (data) {
  FeaturePlot(DK118, 
              features = data, 
              pt.size = 0.5,
              combine = TRUE, cols = c("lightgrey", "red"), min.cutoff = 0)
}
myplots <- lapply(colnames(DK118@meta.data)[36:ncol(DK118@meta.data)], plot_data)
ggsave(filename = "DK118_FeaturePlot_pathway_scores.pdf", dpi = 700, plot = wrap_plots(myplots), width = 17, height = 15)