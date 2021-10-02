#########################################
#          ADT data analysis            #
#           Markers & Plots             #
#            SJG026 (DK114)             #
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

#ADT analysis of DK114 (SJG026)
DK114 <- readRDS("DK114_Seurat_split.integrated.RDS")
ADT_markers_KS(DK114)
DefaultAssay(DK114) <- "TMM"

#Clusters - Visualizations
Idents(DK114) <- "adtClusterID"

#Creating directory for output files
dir.create(file.path("SJG026_DK114"), showWarnings = FALSE)
setwd(file.path("SJG026_DK114"))

#Heatmaps & Feature plot + Violin plot - Clusters - top 10
for(object in as.vector(grep("ADT_markers_cluster_",names(.GlobalEnv),value=TRUE))){
  markers <- eval(as.symbol(object))
  ggsave(DoHeatmap(DK114, slot = "scale.data", assay = "TMM", features = rownames(head(markers[order(-markers$avg_diff),], n = 10))) + labs(title = paste0("SJG026_DK114 top 10 ",object)),
         filename = paste0(object, "_heatmap_top10_SJG026_DK114.pdf"), dpi = 700)
  ggsave(DFV_plot(seurat_object = DK114, 
                  slot = "scale.data", 
                  features = rownames(head(markers[order(-markers$avg_diff),], n = 10)), 
                  title = paste0("SJG026_DK114 top 10 ",object), subtitle = "ADT cluster-based analysis",
                  reduction_1 = "umap", reduction_2 = "umap_tmm"),
         limitsize = FALSE, width = 25, height = 53.6, dpi = 700, #5.36*10
         filename = paste0(object, "_Feat&Vln_top10_SJG026_DK114.pdf"))
}

ggsave(DoHeatmap(DK114, assay = "TMM", slot = "scale.data", 
                 features = c(rownames(ADT_markers_cluster_0[order(ADT_markers_cluster_0$comparison),]), 
                              rownames(ADT_markers_cluster_1[order(ADT_markers_cluster_1$comparison),]),
                              rownames(ADT_markers_cluster_2[order(ADT_markers_cluster_2$comparison),]),
                              rownames(ADT_markers_cluster_3[order(ADT_markers_cluster_3$comparison),]),
                              rownames(ADT_markers_cluster_4[order(ADT_markers_cluster_4$comparison),])), 
                 group.bar = TRUE, 
                 group.by = "adtClusterID"), filename = "ADT_heatmap_markers_clusters_SJG026_DK114.pdf", dpi = 700, height = 15)

#Cluster-specific markers analysis 

ggsave(DoHeatmap(DK114, assay = "TMM", slot = "scale.data", 
          features = c(rownames(ADT_cluster_1_specific_markers[order(ADT_cluster_1_specific_markers$`1vs0_avg_diff`),]), 
                       rownames(ADT_cluster_2_specific_markers[order(ADT_cluster_2_specific_markers$`2vs0_avg_diff`),]), 
                       rownames(ADT_cluster_3_specific_markers[order(ADT_cluster_3_specific_markers$`3vs0_avg_diff`),])), 
          group.bar = TRUE, 
          group.by = "adtClusterID"), filename = "ADT_heatmap_cluster-specific_markers_SJG026_DK114.pdf", dpi = 700)

ggsave(DFV_plot(seurat_object = DK114, 
                slot = "scale.data", 
                features = c(rownames(ADT_cluster_1_specific_markers), rownames(ADT_cluster_2_specific_markers), rownames(ADT_cluster_3_specific_markers)), 
                title = "SJG026_DK114 Kolmogorov-Smirnov ADT cluster-specific markers among pairwise comparisons", subtitle = "ADT cluster-based analysis",
                reduction_1 = "umap", reduction_2 = "umap_tmm"),
       limitsize = FALSE, width = 25, height = 5.36*length(c(rownames(ADT_cluster_1_specific_markers), rownames(ADT_cluster_2_specific_markers), rownames(ADT_cluster_3_specific_markers))), 
       filename = "ADT_cluster-specific_markers_Feature&Vln_SJG026_DK114.pdf")

#Conditions - Visualizations
Idents(DK114) <- "orig.ident"
#Conditions re-assignment to have the desired order of the groups
levels(DK114) <- c("NT", "AG", "R", "RAG")

#Heatmaps & Feature plot + Violin plot- Conditions - top 10
for(object in as.vector(grep("ADT_markers_condition_",names(.GlobalEnv),value=TRUE))){
  markers <- eval(as.symbol(object))
  ggsave(DoHeatmap(DK114, features = rownames(head(markers[order(-markers$avg_diff),], n = 10))) + labs(title = paste0("SJG026_DK114 ",object)),
         filename = paste0(object, "_heatmap_top10_SJG026_DK114.pdf"), dpi = 700)
  ggsave(DFV_plot(seurat_object = DK114, 
                  slot = "scale.data", 
                  features = rownames(head(markers[order(-markers$avg_diff),], n = 10)), 
                  title = paste0("SJG026_DK114 top 10 ",object), subtitle = "ADT condition-based analysis",
                  reduction_1 = "umap", reduction_2 = "umap_tmm"),
         limitsize = FALSE, width = 25, height = 53.6, dpi = 700, #5.36*10
         filename = paste0(object, "_Feat&Vln_top10_SJG026_DK114.pdf"))
}

ggsave(DoHeatmap(DK114, assay = "TMM", slot = "scale.data", 
                 features = c(rownames(ADT_markers_condition_NT[order(ADT_markers_condition_NT$comparison),]), 
                              rownames(ADT_markers_condition_AG[order(ADT_markers_condition_AG$comparison),]),
                              rownames(ADT_markers_condition_R[order(ADT_markers_condition_R$comparison),]),
                              rownames(ADT_markers_condition_RAG[order(ADT_markers_condition_RAG$comparison),])), 
                 group.bar = TRUE), filename = "ADT_heatmap_markers_conditions_SJG026_DK114.pdf", dpi = 700, height = 15)


#Condition-specific markers analysis

ggsave(DoHeatmap(DK114, assay = "TMM", slot = "scale.data", 
          features = c(rownames(ADT_condition_NT_specific_markers), 
                       rownames(ADT_condition_AG_specific_markers),
                       rownames(ADT_condition_R_specific_markers),
                       rownames(ADT_condition_RAG_specific_markers)), 
          group.bar = TRUE), filename = "ADT_heatmap_condition-specific_markers_SJG026_DK114.pdf", dpi = 700)

ggsave(DFV_plot(seurat_object = DK114, 
                slot = "scale.data", 
                features = c(rownames(ADT_condition_NT_specific_markers), 
                             rownames(ADT_condition_AG_specific_markers),
                             rownames(ADT_condition_R_specific_markers),
                             rownames(ADT_condition_RAG_specific_markers)), 
                title = "SJG026_DK114 Kolmogorov-Smirnov ADT condition-specific markers among pairwise comparisons", subtitle = "ADT condition-based analysis",
                reduction_1 = "umap", reduction_2 = "umap_tmm"),
       limitsize = FALSE, width = 25, height = 5.36*length(c(rownames(ADT_condition_NT_specific_markers), 
                                                             rownames(ADT_condition_AG_specific_markers),
                                                             rownames(ADT_condition_R_specific_markers),
                                                             rownames(ADT_condition_RAG_specific_markers))), 
       filename = "ADT_condition-specific_markers_Feature&Vln_SJG026_DK114.pdf")

#Saving csv data
for (file in as.vector(grep("ADT_markers_cluster",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG026_DK114.csv"))
}

for (file in as.vector(grep("ADT_markers_condition",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG026_DK114.csv"))
}

for (file in as.vector(grep("ADT_cluster",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG026_DK114.csv"))
}

for (file in as.vector(grep("ADT_condition_",names(.GlobalEnv),value=TRUE))){
  write.csv(eval(as.symbol(file)), file = paste0(file, "_SJG026_DK114.csv"))
}

#Single-cell pathway scores
SJG026_DK114_pathway_scores <- sc_pathway_scores(DK114)

#Saving csv file
write.csv(SJG026_DK114_pathway_scores, file = "ADT_pathway_scores_SJG026_DK114.csv")

#Packages for plots
library(reshape2)
library(viridis)
library(patchwork)
library(ggplot2)

plot_SJG026_DK114_pathway_scores <- melt(SJG026_DK114_pathway_scores)
colnames(plot_SJG026_DK114_pathway_scores) <- c("Conditions", "ADT Cluster", "Sample", "Pathway", "Score")

#Violin plot
ggsave(ggplot(plot_SJG026_DK114_pathway_scores, aes(Pathway, Score, fill=Conditions)) + 
         geom_violin(position = position_dodge(width = 0.9))+
         labs(title="SJG026_DK114 ADT pathway scores", y = "Score")+
         geom_violin(trim=FALSE)+
         scale_fill_brewer(palette="RdBu") + theme_bw() +
         theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank()), width = 15, filename = "ADT_pathway_scores_violinplot_SJG026_DK114.pdf", limitsize = FALSE)

ggsave(ggplot(plot_SJG026_DK114_pathway_scores, aes(Pathway, Score, fill = `ADT Cluster`)) + 
         geom_violin(position = position_dodge(width = 0.9))+
         labs(title="SJG026_DK114 ADT pathway scores", y = "Score")+
         geom_violin(trim=FALSE)+
         scale_fill_brewer(palette="RdBu") + theme_bw() +
         theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank()), width = 15, filename = "ADT_pathway_scores_violinplot_clusters_SJG026_DK114.pdf", limitsize = FALSE)

#Dot plots
dotplots_clusters <- dotplot_pathway_scores_clusters(pathway_scores = SJG026_DK114_pathway_scores, clusters = c(0, 1, 2, 3, 4), patient = "SJG026_DK114")
ggsave(filename = "ADT_Dotplot_pathwayscores_clusters_SJG026_DK114.pdf", plot = dotplots_clusters, dpi = 700, width = 4.4)

dotplots_conditions <- dotplot_pathway_scores_conditions(pathway_scores = SJG026_DK114_pathway_scores, patient = "SJG026_DK114")
ggsave(filename = "ADT_Dotplot_pathwayscores_conditions_SJG026_DK114.pdf", plot = dotplots_conditions, dpi = 700, width = 4.4)

save.image("SJG026_DK114_ADT.RData")

#Adding metadata referring to the pathway scores
for (pathway in colnames(SJG026_DK114_pathway_scores)[1:9]){
  DK114 <- AddMetaData(DK114, metadata = SJG026_DK114_pathway_scores[, pathway], col.name = pathway)
}
#Feature plot
plot_data = function (data) {
  FeaturePlot(DK114, 
              features = data, 
              pt.size = 0.5,
              combine = TRUE, cols = c("lightgrey", "red"), min.cutoff = 0)
}
myplots <- lapply(colnames(DK114@meta.data)[36:ncol(DK114@meta.data)], plot_data)
ggsave(filename = "DK114_FeaturePlot_pathway_scores.pdf", dpi = 700, plot = wrap_plots(myplots), width = 17, height = 15)
