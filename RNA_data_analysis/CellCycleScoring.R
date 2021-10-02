#########################################
#          RNA data analysis            #
#          cell cycle scoring           #
#########################################

library(readxl)
library(Seurat)

#Loading Seurat objects of the two patients' datasets
DK114 <- readRDS(file = 'DK114_Seurat_split.integrated.RDS')
DK118 <- readRDS(file = 'DK118_Seurat_split.integrated.RDS')
DK104 <- readRDS(file = 'DK104_Seurat_split.integrated.RDS')
DK126 <- readRDS(file = 'DK126_Seurat_split.integrated.RDS')

#List of G1/S and G2/M phase genes from Dominguez et al. 2016 (https://www.nature.com/articles/cr201684#MOESM19)
cellcycle_genes_dominguez <- as.data.frame(read_xlsx(path = 'Dominguez et al 2016.xlsx', skip = 4))

#Core cell cycle genes from Dominguez et al.
core_cellcycle_genes_dominguez <- cellcycle_genes_dominguez[cellcycle_genes_dominguez$`Core 67`=='YES',]
G1S_genes_Dominguez <- core_cellcycle_genes_dominguez[core_cellcycle_genes_dominguez$Stage=='G1-S',]$`Gene Name`
G2M_genes_Dominguez <- core_cellcycle_genes_dominguez[core_cellcycle_genes_dominguez$Stage=="G2-M",]$`Gene Name`

#Puram et al.programs
puram_data <- as.data.frame(read_xlsx('Puram et al.xlsx'))
colnames(puram_data) <- c('Cluster', 'Inferred Annotation', 'Patient (MEEI)', seq(50))

#Selecting the genes in the G1/S cell cycle program from Puram et al.
G1S_genes_Puram <- as.vector(as.matrix(subset(na.omit(puram_data[puram_data$`Inferred Annotation`=='Cell cycle(G1/S)',]), select = -c(1:3))))

#Selecting the genes in the G2/M cell cycle program from Puram et al.
G2M_genes_Puram <- as.vector(as.matrix(subset(na.omit(puram_data[puram_data$`Inferred Annotation`=='Cell cycle(G2/m)',]), select = -c(1:3))))

#Selecting the genes in the G1/S + G2/M cell cycle program from Puram et al.
G1S_G2M_genes_Puram <- as.vector(as.matrix(subset(na.omit(puram_data[puram_data$`Inferred Annotation`=='Cell cycle(G1/S+G2/M)',]), select = -c(1:3))))

#Kinker et al. programs: RHPs from supplementary table 4
kinker_data <- as.data.frame(read_xlsx('Kinker et al.xlsx', sheet = 'Table S4', skip = 3))
G1S_genes_Kinker <- as.vector(kinker_data$`Cell Cycle - G1/S`)[!is.na(as.vector(kinker_data$`Cell Cycle - G1/S`))]
G2M_genes_Kinker <- as.vector(kinker_data$`Cell Cycle - G2/M`)[!is.na(as.vector(kinker_data$`Cell Cycle - G2/M`))]

G1S_genes <- unique(c(G1S_genes_Dominguez, G1S_genes_Kinker, G1S_genes_Puram, G1S_G2M_genes_Puram))
G2M_genes <- unique(c(G2M_genes_Dominguez, G2M_genes_Kinker, G2M_genes_Puram, G1S_G2M_genes_Puram))

DefaultAssay(DK114) <- "SCT"
DK114 <- CellCycleScoring(DK114, s.features = G1S_genes, g2m.features = G2M_genes, set.ident = TRUE)

DefaultAssay(DK118) <- "SCT"
DK118 <- CellCycleScoring(DK118, s.features = G1S_genes, g2m.features = G2M_genes, set.ident = TRUE)

DefaultAssay(DK104) <- "SCT"
DK104 <- CellCycleScoring(DK104, s.features = G1S_genes, g2m.features = G2M_genes, set.ident = TRUE)

DefaultAssay(DK126) <- "SCT"
DK126 <- CellCycleScoring(DK126, s.features = G1S_genes, g2m.features = G2M_genes, set.ident = TRUE)

library(cowplot)

levels(DK114) <- levels(DK118) <- levels(DK104) <- levels(DK126) <- c("G1", "G2M", "S")

p1 <- DimPlot(DK114, pt.size = 0.5)
p2 <- DimPlot(DK118, pt.size = 0.5)
p3 <- DimPlot(DK104)
p4 <- DimPlot(DK126)

ggsave(plot_grid(p1, p2, p3, p4), dpi = 300, filename = "CellCycleScoring_UMAP.pdf", width = 10, height = 8)
