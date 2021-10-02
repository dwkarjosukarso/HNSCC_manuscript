#########################################
#          RNA data analysis            #
#           Markers & Plots             #
#               dataset                 #
#########################################

#This source contains the functions used for the Wilcoxon Rank Sum test-based 
#differential expression analysis of RNA data.
source("RNA_DE_Wilcoxon.R")

#DK114
#Loading the Seurat object
DK114 <- readRDS("DK114_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG026_DK114"), showWarnings = FALSE)
setwd(file.path("SJG026_DK114"))

#RNA Cluster-based and condition-based DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK114, assay = "SCT", slot = "scale.data", patient = "SJG026_DK114")
scRNA_DE_analysis(type = "conditions", seurat_object = DK114, assay = "SCT", slot = "scale.data", patient = "SJG026_DK114")

setwd("..")
rm(list = ls())

source("RNA_DE_Wilcoxon.R")

#DK104
#Loading the Seurat object
DK104 <- readRDS("DK104_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG026_DK104"), showWarnings = FALSE)
setwd(file.path("SJG026_DK104"))

#RNA Cluster-based and condition-based DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK104, assay = "SCT", slot = "scale.data", patient = "SJG026_DK104")
scRNA_DE_analysis(type = "conditions", seurat_object = DK104, assay = "SCT", slot = "scale.data", patient = "SJG026_DK104")

setwd("..")
rm(list = ls())

source("RNA_DE_Wilcoxon.R")

#DK118
#Loading the Seurat object
DK118 <- readRDS("DK118_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG017_DK118"), showWarnings = FALSE)
setwd(file.path("SJG017_DK118"))

#Cluster-based and condition-based DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK118, assay = "SCT", slot = "scale.data", patient = "SJG017_DK118")
scRNA_DE_analysis(type = "conditions", seurat_object = DK118, assay = "SCT", slot = "scale.data", patient = "SJG017_DK118")

setwd("..")
rm(list = ls())

source("RNA_DE_Wilcoxon.R")

#DK126
#Loading the Seurat object
DK126 <- readRDS("DK126_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG017_DK126"), showWarnings = FALSE)
setwd(file.path("SJG017_DK126"))

#Cluster-based and condition-based DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK126, assay = "SCT", slot = "scale.data", patient = "SJG017_DK126")
scRNA_DE_analysis(type = "conditions", seurat_object = DK126, assay = "SCT", slot = "scale.data", patient = "SJG017_DK126")

setwd("..")
rm(list = ls())

#RNA markers as computed in the ADT clusters

source("RNA_DE_Wilcoxon.R")

#DK114
#Loading the Seurat object
DK114 <- readRDS("DK114_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG026_ADT_DK114"), showWarnings = FALSE)
setwd(file.path("SJG026_ADT_DK114"))

#ADT Cluster-based RNA data DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK114, assay = "SCT", slot = "scale.data", idents_1 = "adtClusterID", patient = "adtID_SJG026_DK114")

setwd("..")
rm(list = ls())

source("RNA_DE_Wilcoxon.R")

#DK104
#Loading the Seurat object
DK104 <- readRDS("DK104_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG026_ADT_DK104"), showWarnings = FALSE)
setwd(file.path("SJG026_ADT_DK104"))

#ADT Cluster-based RNA data DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK104, assay = "SCT", slot = "scale.data", idents_1 = "adtClusterID", patient = "adtID_SJG026_DK104")

setwd("..")
rm(list = ls())

source("RNA_DE_Wilcoxon.R")

#DK118
#Loading the Seurat object
DK118 <- readRDS("DK118_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG017_ADT_DK118"), showWarnings = FALSE)
setwd(file.path("SJG017_ADT_DK118"))

#ADT Cluster-based RNA data DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK118, assay = "SCT", slot = "scale.data", idents_1 = "adtClusterID", patient = "adtID_SJG017_DK118")

setwd("..")
rm(list = ls())

source("RNA_DE_Wilcoxon.R")

#DK126
#Loading the Seurat object
DK126 <- readRDS("DK126_Seurat_split.integrated.RDS")

#Creating directory for output files
dir.create(file.path("SJG017_ADT_DK126"), showWarnings = FALSE)
setwd(file.path("SJG017_ADT_DK126"))

#ADT Cluster-based RNA data DE analysis: the output is also found in the global environment
scRNA_DE_analysis(type = "clusters", seurat_object = DK126, assay = "SCT", slot = "scale.data", idents_1 = "adtClusterID", patient = "adtID_SJG017_DK126")





