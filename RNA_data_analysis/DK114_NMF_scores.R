#################################
#       NMF data analysis       #   
#     single-cell NMF robust    #
#        programs scores        #
#    only programs from same    #
#    experimental conditions    #
#       no negative values      #
#           transformed         #
#         SJG-026 DK-114        #
#################################

#Loading the required packages
library(ggplot2)
library(viridis)
library(reshape)
library(Seurat)
library(patchwork)
library(readxl)

source("NMF_single_cell_robust_program_scores.R")
source("DotplotProgramsScores.R")

#Reading programs
NMF_programs <- as.list(read_xlsx(path = "Robust_programs_KM.xlsx", col_names = TRUE))

#Omit NA
NMF_programs <- lapply(NMF_programs, function(x) as.vector(na.omit(x)))

#DK114 ####
#Loading the Seurat object
DK114 <- readRDS("DK114_Seurat_split.integrated.RDS")
std_DK114 <- as.data.frame(GetAssayData(DK114, assay = "SCT", slot = "data"))
std_DK114 <- std_DK114 - rowMeans(std_DK114)

#Creating directory for output files
dir.create(file.path("SJG026_DK114_NMF_scores"), showWarnings = FALSE)
setwd(file.path("SJG026_DK114_NMF_scores"))

#Computing the scores for each of the NMF robust programs
sink("DK114_NMF_robust_programs_scores_outcome.txt")
DK114_NMF_programs_scores <- sc_NMF_robust_programs_scores(std_DK114, NMF_programs)
sink()

#Preparing the data frame for the plots
DK114_NMF_programs_scores <- as.data.frame(DK114_NMF_programs_scores)
DK114_NMF_programs_scores$CellID <- rownames(DK114_NMF_programs_scores)
DK114_NMF_programs_scores$Condition <- DK114[["orig.ident"]]
DK114_NMF_programs_scores$ADT_Cluster <- DK114[["adtClusterID"]]
DK114_NMF_programs_scores$RNA_Cluster <- DK114[["rnaClusterID"]]

#Saving the scores to csv
write.csv(x = as.matrix(DK114_NMF_programs_scores), 
          file = "SJG026_DK114_NMF_robust_programs_scores.csv")

#Saving the scores as rds
saveRDS(object = DK114_NMF_programs_scores, 
        file = "SJG026_DK114_NMF_robust_programs_scores.RDS")

#RNA Clusters dot plot
ggsave(dotplot_programs_scores_RNAclusters(programs = NMF_programs,
                                           RNAclusters = c(0, 1, 2, 3, 4, 5), 
                                           dataset_programs_scores = DK114_NMF_programs_scores, 
                                           no.programs = 6, 
                                           patient_dataset = "SJG026 DK114"), 
       filename = "DK114_NMF_robust_programs_scores_RNAclusters.pdf", dpi = 700, width = 4, height = 8)

#ADT Clusters dot plot
ggsave(dotplot_programs_scores_ADTclusters(programs = NMF_programs,
                                           ADTclusters = c(0, 1, 2, 3, 4), 
                                           dataset_programs_scores = DK114_NMF_programs_scores, 
                                           no.programs = 6, 
                                           patient_dataset = "SJG026 DK114"), 
       filename = "DK114_NMF_robust_programs_scores_ADTclusters.pdf", dpi = 700, width = 4, height = 8)

#Conditions dot plot
ggsave(dotplot_programs_scores_conditions(programs = NMF_programs,
                                          conditions = c("NT", "AG", "R", "RAG"), 
                                          dataset_programs_scores = DK114_NMF_programs_scores, 
                                          no.programs = 6, 
                                          patient_dataset = "SJG026 DK114"), 
       filename = "DK114_NMF_robust_programs_scores_conditions.pdf", dpi = 700, width = 4, height = 8)

#Adding metadata referring to the NMF robust programs scores
for (program in colnames(DK114_NMF_programs_scores)[1:6]){
  DK114 <- AddMetaData(DK114, metadata = DK114_NMF_programs_scores[, program], col.name = program)
}

#Feature plot
plot_data = function (data) {
  FeaturePlot(DK114, 
              features = data, 
              pt.size = 0.5,
              combine = TRUE, cols = c("lightgrey", "red"), min.cutoff = 0)
}

plot_data2 = function (data) {
  FeaturePlot(DK114, 
              features = data, 
              pt.size = 0.5,
              combine = TRUE, cols = c("lightgrey", "red"), min.cutoff = 0, reduction = "umap_tmm")
}

myplots <- lapply(colnames(DK114@meta.data)[36:ncol(DK114@meta.data)], plot_data)
ggsave(filename = "DK114_FeaturePlot_NMF_programs_scores.pdf", dpi = 700, plot = wrap_plots(myplots), width = 17, height = 12)
myplots <- lapply(colnames(DK114@meta.data)[36:ncol(DK114@meta.data)], plot_data2)
ggsave(filename = "DK114_FeaturePlot_NMF_programs_scores_TMMumap.pdf", dpi = 700, plot = wrap_plots(myplots), width = 17, height = 12)


matrix_for_programs <- GetAssayData(DK114, assay = "SCT", slot = "data")
matrix_for_programs <- matrix_for_programs[intersect(as.vector(unlist(NMF_programs)), rownames(matrix_for_programs)),]
write.csv(matrix_for_programs, file = "DK114_SCT_data_NMF_NMF_programs_genes.csv")

matrix_for_programs <- GetAssayData(DK114, assay = "SCT", slot = "scale.data")
matrix_for_programs <- matrix_for_programs[intersect(as.vector(unlist(NMF_programs)), rownames(matrix_for_programs)),]
write.csv(matrix_for_programs, file = "DK114_SCT_scale.data_NMF_NMF_programs_genes.csv")
