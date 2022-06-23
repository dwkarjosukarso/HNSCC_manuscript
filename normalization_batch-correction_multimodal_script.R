# set working directory
setwd(DIR)

# load RNA data
# load data
RNA.counts <- read.csv(RNA_COUNT_TABLE, header= TRUE, row.names=1)
BC.counts <- read.csv(BC_COUNT_TABLE, header= TRUE, row.names=1)

# check number of cells in both RNA and Ab-barcode count table
ncol(RNA.counts)
ncol(BC.counts)

# if both ncol results are identical, no need to match the number of cells in both datasets
# if not, match them
# matching number of cells in both datasets and fill in with 0 on empty cells
# adjust colnames of both count tables
colnames(RNA.counts) <- substr(colnames(RNA.counts), 1, nchar(colnames(RNA.counts))-14)
# fill in the missing cell with 0
RNA.counts[setdiff(names(BC.counts), names(RNA.counts))] <- 0
BC.counts[setdiff(names(RNA.counts), names(BC.counts))] <- 0
write.csv(RNA.counts, FILENAME)
write.csv(BC.counts, FILENAME)

#Seurat
library(Seurat)

# Create Seurat object with minimal filtering to explore the data
Seurat <- CreateSeuratObject(counts =RNA.counts, min.cells = 0, min.features  = 0, project = "DK114")

# Add ADT data to Seurat object
Seurat[["ADT"]] <- CreateAssayObject(counts= BC.counts)

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  
# For non-UMI data, nUMI represents the sum of the non-normalized values within a cell. 
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and non-log-normalized counts. 
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
mito.genes <- grep(pattern= "MT-", x= rownames(Seurat@assays[["RNA"]]), value= TRUE)
percent.mito <- Matrix::colSums(Seurat@assays[["RNA"]][mito.genes, ])/Matrix::colSums(Seurat@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
# Seurat v2 function, but shows compatibility in Seurat v3
Seurat <- AddMetaData(object = Seurat, metadata = percent.mito, col.name = "percent.mito") 

# QC plot 1
pdf(FILENAME)
VlnPlot(object = Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# Feature Scatter
pdf(FILENAME)
FeatureScatter(object = Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf(FILENAME)
VlnPlot(object = Seurat, features = c("nFeature_ADT", "nCount_ADT"), ncol = 2)
dev.off()

pdf(FILENAME)
FeatureScatter(object = Seurat, feature1 = "nCount_ADT", feature2 = "nFeature_ADT")
dev.off()

# subset data to remove cells that do not meet the following criteria
# has less than 500 genes detected (keep nFeature_RNA > 500)
# has MT genes more than 50% (keep percent.mito < 0.5)
# has less than 20 ADT detected (keep nFeature_ADT > 20)
Seurat_filtered <- subset(Seurat, subset= nFeature_RNA >500 & nFeature_RNA <6000 & percent.mito > -Inf & percent.mito < 0.5 & nFeature_ADT > 20)
RNA.filter <- as.data.frame(GetAssayData(Seurat_filtered, assay= "RNA", slot= "counts"))
BC.filter  <- as.data.frame(GetAssayData(Seurat_filtered, assay= "ADT", slot= "counts"))
write.csv(RNA.filter, FILENAME)
write.csv(BC.filter, FILENAME)

# QC plot 2
pdf(FILENAME)
VlnPlot(object = Seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")
dev.off()

# Feature Scatter
pdf(FILENAME)
FeatureScatter(object = Seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
dev.off()

pdf(FILENAME)
VlnPlot(object = Seurat_filtered, features = c("nFeature_ADT", "nCount_ADT"), ncol = 2, group.by = "orig.ident")
dev.off()

pdf(FILENAME)
FeatureScatter(object = Seurat_filtered, feature1 = "nCount_ADT", feature2 = "nFeature_ADT", group.by = "orig.ident")
dev.off()

# count number of cells passed QC from each condition
NT <- length(colnames(Seurat_filtered)[grepl("^NT_", colnames(Seurat_filtered))])
AG <- length(colnames(Seurat_filtered)[grepl("^AG_", colnames(Seurat_filtered))])
R <- length(colnames(Seurat_filtered)[grepl("^R_", colnames(Seurat_filtered))])
RAG <- length(colnames(Seurat_filtered)[grepl("^RAG_", colnames(Seurat_filtered))])

number_of_cells <- data.frame(condition = c("NT", "AG", "R", "RAG"), n= c(NT, AG, R, RAG))
library(ggplot2)

number_of_cells_plot <- ggplot(number_of_cells, aes(x= condition, y= n)) + geom_col()
ggsave(FILENAME, number_of_cells_plot)

# RNA analysis
# RNA data normalization and batch effect correction
# define batches
Seurat_filtered@meta.data$plate <- rownames(Seurat_filtered@meta.data)
Seurat_filtered@meta.data$plate <- substr(Seurat_filtered@meta.data$plate,1,nchar(Seurat_filtered@meta.data$plate)-5)
Seurat_filtered@meta.data$plate <- gsub("^.*\\_", "plate", Seurat_filtered@meta.data$plate)
# batch correction with Seurat Integration & SCTrasnform normalization
Seurat_split <- SplitObject(Seurat_filtered, split.by = "plate")
for (i in 1:length(Seurat_split)) {
  Seurat_split[[i]] <-SCTransform(Seurat_split[[i]], assay= "RNA", verbose = FALSE)
}
Seurat_split.features <- SelectIntegrationFeatures(object.list = Seurat_split, nfeatures = 3000)
Seurat_split <- PrepSCTIntegration(object.list = Seurat_split, anchor.features = Seurat_split.features, 
                                   verbose = FALSE)
Seurat_split.anchors <- FindIntegrationAnchors(object.list = Seurat_split, normalization.method = "SCT", 
                                               anchor.features = Seurat_split.features, verbose = FALSE)
Seurat_split.integrated <- IntegrateData(anchorset = Seurat_split.anchors, normalization.method = "SCT", 
                                         verbose = FALSE)

# PCA of batch corrected edata
Seurat_split.integrated <- RunPCA(Seurat_split.integrated, verbose = FALSE)
PCAPlot(Seurat_split.integrated, group.by="orig.ident")
ElbowPlot(Seurat_split.integrated, ndims= 50)

# Clustering (use ndims based on Elbowplot results)
Seurat_split.integrated <- FindNeighbors(Seurat_split.integrated, assay= "SCT", dims= 1:25)
resolution <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (i in resolution){
  Seurat_split.integrated <- FindClusters(Seurat_split.integrated, assay= "SCT", resolution= i)
}

# check clustering stability
library(clustree)
clustree(Seurat_split.integrated, prefix="integrated_snn_res.")

# UMAP
# choose an appropriate resolution based on clustree results
Seurat_split.integrated <- FindClusters(Seurat_split.integrated, assay= "SCT", resolution= 0.3)
Seurat_split.integrated <- RunUMAP(Seurat_split.integrated, assay= "SCT", dims= 1:25)
DimPlot(Seurat_split.integrated, reduction= "umap", label= TRUE)
DimPlot(Seurat_split.integrated, reduction= "umap", label= TRUE, group.by= "orig.ident")
DimPlot(Seurat_split.integrated, reduction= "umap", label= TRUE, group.by= "plate", split.by= "orig.ident")
DimPlot(Seurat_split.integrated, reduction= "umap", label= TRUE, split.by= "orig.ident")

# stacked bargraph
# percent condition per cluster, change resolution based on the one used for generating UMAP
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(wesanderson)
stacked <- apply(table(Seurat_split.integrated@meta.data$integrated_snn_res.0.3, Seurat_split.integrated@meta.data$orig.ident), 
                 1, function(x){x*100/sum(x,na.rm=T)})
stacked <- as.data.frame(stacked)
stacked$condition <- rownames(stacked)
stacked2 <- melt(stacked)
stacked_bar <- ggplot(stacked2, aes(x= variable, y= value, fill= factor(condition, levels=c("NT", "AG", "R", "RAG")))) + 
  geom_bar(stat= "identity") + theme_classic() + scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")              
stacked_bar

# percent plate per cluster
Seurat_split.integrated@meta.data$replicate <- paste0(Seurat_split.integrated@meta.data$orig.ident, "_", Seurat_split.integrated@meta.data$plate)
stacked_plate <- apply(table(Seurat_split.integrated@meta.data$integrated_snn_res.0.3, Seurat_split.integrated@meta.data$replicate), 1, function(x){x*100/sum(x,na.rm=T)})
stacked_plate <- as.data.frame(stacked_plate)
stacked_plate$condition <- rownames(stacked_plate)
stacked_plate2 <- melt(stacked_plate)
stacked_plate_bar <- ggplot(stacked_plate2, aes(x= variable, y= value, 
                                                fill= factor(condition, levels=c("NT_plate1", "NT_plate2", "NT_plate3", "AG_plate1", "AG_plate2", "AG_plate3", 
                                                                                 "R_plate1", "R_plate2", "R_plate3", "RAG_plate1", "RAG_plate2", "RAG_plate3")))) + 
  geom_bar(stat= "identity") + theme_classic() +  scale_fill_brewer(palette= "Set3") + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")
stacked_plate_bar

# save cluster id
Seurat_split.integrated[["rnaClusterID"]] <- Idents(Seurat_split.integrated)


# Analysis of ADT data
# TMM normalization of ADT data
# generate meta data
library(dplyr)
library(tidyr)
library(tibble)
rownames(BC.filter) <- substr(rownames(BC.filter), 1, nchar(rownames(BC.filter))-11)

BC.long <- BC.filter %>% 
  rownames_to_column(var = 'ab_name') %>% 
  pivot_longer(-ab_name, names_to = 'sample_id', values_to = 'ab_count_raw') %>% 
  mutate(plate_number = case_when(grepl("^NT_1", sample_id) ~ "plate_1",
                                  grepl("^NT_2", sample_id) ~ "plate_2",
                                  grepl("^NT_3", sample_id) ~ "plate_3",
                                  grepl("^AG_1", sample_id) ~ "plate_1",
                                  grepl("^AG_2", sample_id) ~ "plate_2",
                                  grepl("^AG_3", sample_id) ~ "plate_3",
                                  grepl("^R_1", sample_id) ~ "plate_1",
                                  grepl("^R_2", sample_id) ~ "plate_2",
                                  grepl("^R_3", sample_id) ~ "plate_3",
                                  grepl("^RAG_1", sample_id) ~ "plate_1",
                                  grepl("^RAG_2", sample_id) ~ "plate_2",
                                  grepl("^RAG_3", sample_id) ~ "plate_3")) %>%
  mutate(treatment = case_when(grepl("^NT_", sample_id) ~ "NT",
                               grepl("^AG_", sample_id) ~ "AG",
                               grepl("^R_", sample_id) ~ "R",
                               grepl("^RAG_", sample_id) ~ "RAG"))


# To separate count data from metadata a file with antibody type information (ab.data), metadata (metadata), and a file with library size information (sample.annot) were created based on the sample_id. 
ab.data <- BC.long %>% 
  select(ab_name) %>% 
  mutate(ab_type = case_when(grepl("Phospho", ab_name) ~ 'phospho')) %>% 
  distinct()

ab.data$ab_type[is.na(ab.data$ab_type)] <- "total"

metadata <- BC.long %>% 
  select("sample_id", "plate_number", "treatment") %>% 
  distinct(sample_id, plate_number, treatment)

sample.annot <- BC.long %>% 
  group_by(sample_id) %>% 
  dplyr::mutate(lib_size = sum(ab_count_raw)) %>% 
  select(sample_id, plate_number, treatment, lib_size) %>% 
  ungroup() %>% 
  distinct()

# Calculate normfactors followed by effective library size (lib_size*normfactor). Divide raw counts by effective library size.
# BC.tmm contains normalized counts

normfactors <- edgeR::calcNormFactors(BC.filter, sumTrim = 0.05, logratioTrim = 0) %>% 
  enframe("sample_id", "normfactor")

BC.tmm <- BC.long %>%
  left_join(normfactors, by = "sample_id") %>%
  left_join(sample.annot, by = c('sample_id', 'plate_number', 'treatment')) %>%
  mutate(ab_count_normalized = ab_count_raw/(lib_size*normfactor))

## log transformation of TMM nromalized counts
## df.tmm.log contains log transformed normalized counts
BC.tmm.log <- BC.tmm %>% mutate(ab_count_log = log10(ab_count_normalized+1)) %>%
  select(-ab_count_normalized)

## df.tmm.log.batch contains normalized and log transformed counts after batch effect removal
batch.corrected.tmm.log <- select(BC.tmm.log, sample_id, ab_name, ab_count_log) %>% 
  pivot_wider(names_from = ab_name, values_from = ab_count_log) %>% 
  column_to_rownames("sample_id") %>% 
  as.matrix()

batch = column_to_rownames(sample.annot, var = "sample_id")[
  rownames(batch.corrected.tmm.log), 
  "plate_number"
  ]

batch.corrected.tmm.log <- sva::ComBat(t(batch.corrected.tmm.log), batch = batch) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  pivot_longer(-sample_id, names_to = "ab_name", values_to = "ab_count_abr")

BC.tmm.log.batch <- BC.tmm.log %>%
  select(-ab_count_log) %>%
  left_join(batch.corrected.tmm.log, by = c('ab_name', 'sample_id')) %>% 
  left_join(ab.data, by = 'ab_name')

# re-format to count tables to enter into Seurat
count.tmm <- BC.tmm.log.batch %>% select(ab_name, sample_id, ab_count_abr) %>%
  pivot_wider(names_from = sample_id, values_from = ab_count_abr) %>% 
  column_to_rownames("ab_name")
write.csv(count.tmm, FILENAME)

# enter TMM normalized counts to Seurat
Seurat_split.integrated[["TMM"]] <- CreateAssayObject(counts= count.tmm)

# Find variable features
Seurat_split.integrated <- FindVariableFeatures(Seurat_split.integrated, assay= "TMM")

# PCA
Seurat_split.integrated <- ScaleData(Seurat_split.integrated, assay = "TMM")
Seurat_split.integrated <- RunPCA(Seurat_split.integrated, assay= "TMM", reduction.name = "pca_tmm", reduction.key = "pca_tmm_")
PCAPlot(Seurat_split.integrated, reduction = "pca_tmm", group.by="orig.ident")
ElbowPlot(Seurat_split.integrated, reduction= "pca_tmm", ndims= 50)

# Euclidean distance matrix
# Visualization and clustering using a custom distance matrix in Seurat
tmm.data <- GetAssayData(Seurat_split.integrated, assay= "TMM", slot="data")
tmm.dist <- dist(t(as.matrix(tmm.data)))

# ADT-based clustering 
Seurat_split.integrated[["tmm_snn"]] <- FindNeighbors(tmm.dist)$snn
for (i in resolution){
  Seurat_split.integrated <- FindClusters(Seurat_split.integrated, assay= "TMM", resolution= i, graph.name= "tmm_snn")
}

# check clustering stability
library(clustree)
clustree(Seurat_split.integrated, prefix="tmm_snn_res.")

# UMAP
# choose resoltuion based on clustree
Seurat_split.integrated <- FindClusters(Seurat_split.integrated, assay= "TMM", resolution= 0.7, graph.name= "tmm_snn")
Seurat_split.integrated[["umap_tmm"]] <- RunUMAP(tmm.dist, assay= "TMM", reduction.key= "tmmUMAP_")
DimPlot(Seurat_split.integrated, reduction= "umap_tmm", label= TRUE)
DimPlot(Seurat_split.integrated, reduction= "umap_tmm", label= TRUE, group.by= "orig.ident")
DimPlot(Seurat_split.integrated, reduction= "umap_tmm", label= TRUE, group.by= "plate", split.by = "orig.ident")
DimPlot(Seurat_split.integrated, reduction= "umap_tmm", label= TRUE, split.by = "orig.ident")

# stacked bargraph
# percent condition per cluster, change resolution based on the one used for generating UMAP
BCstacked <- apply(table(Seurat_split.integrated@meta.data$tmm_snn_res.0.7, Seurat_split.integrated@meta.data$orig.ident), 1, function(x){x*100/sum(x,na.rm=T)})
BCstacked <- as.data.frame(BCstacked)
BCstacked$condition <- rownames(BCstacked)
BCstacked2 <- melt(BCstacked)
BCstacked_bar <- ggplot(BCstacked2, aes(x= variable, y= value, fill= factor(condition, levels=c("NT", "AG", "R", "RAG")))) + 
  geom_bar(stat= "identity") + theme_classic() + scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")               
BCstacked_bar

# percent plate per cluster
BCstacked_plate <- apply(table(Seurat_split.integrated@meta.data$tmm_snn_res.0.7, Seurat_split.integrated@meta.data$replicate), 1, function(x){x*100/sum(x,na.rm=T)})
BCstacked_plate <- as.data.frame(BCstacked_plate)
BCstacked_plate$condition <- rownames(BCstacked_plate)
BCstacked_plate2 <- melt(BCstacked_plate)
BCstacked_plate_bar <- ggplot(BCstacked_plate2, aes(x= variable, y= value, 
                                                    fill= factor(condition, levels=c("NT_plate1", "NT_plate2", "NT_plate3", "AG_plate1", "AG_plate2", "AG_plate3", 
                                                                                     "R_plate1", "R_plate2", "R_plate3", "RAG_plate1", "RAG_plate2", "RAG_plate3")))) + 
  geom_bar(stat= "identity") + theme_classic() +  scale_fill_brewer(palette= "Set3") + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")
BCstacked_plate_bar

# save ADT cluster ID
Seurat_split.integrated[["adtClusterID"]] <- Idents(Seurat_split.integrated)

# export RDS file
saveRDS(Seurat_split.integrated, FILENAME)

# CNV-based clustering
CNV_matrix <- read.table(CNV_COUNT_TABLE, row.names=1, header= TRUE)
CNV_matrix <- 2 * CNV_matrix
Seurat_split.integrated[["CNV"]] <- CreateAssayObject(counts= CNV_matrix)
Seurat_split.integrated <- FindVariableFeatures(Seurat_split.integrated, assay="CNV")
Seurat_split.integrated <-ScaleData(Seurat_split.integrated, assay= "CNV")
Seurat_split.integrated <- RunPCA(Seurat_split.integrated, assay= "CNV", reduction.name = "pca_cnv")
PCAPlot(Seurat_split.integrated, reduction = "pca_cnv", group.by= "orig.ident")
ElbowPlot(Seurat_split.integrated, ndims= 50)

# Add cell groupings determined by inferCNV to metadata
# adding inferCNV results to Seurat
cnv_group <- read.table("/scratch/karjosukarso/DK-114/inferCNV/new_infercnv_nogroup_subclusters/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings", sep="\t", header= T)
cnv_group <- cnv_group[grepl("all_observations", cnv_group$cell_group_name),]
rownames(cnv_group) <- cnv_group$cell
cnv_group$cell <- NULL
# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
Seurat_split.integrated <- AddMetaData(object = Seurat_split.integrated, metadata = cnv_group, col.name = "cnv") 
# remove a part of the string in cnv groupings
Seurat_split.integrated$cnv <- substr(Seurat_split.integrated$cnv, 35, 41)

# Euclidean distance matrix
# Visualization and clustering using a custom distance matrix in Seurat
CNV.data <- GetAssayData(Seurat_split.integrated, assay= "CNV", slot="counts")
CNV.dist <- dist(t(as.matrix(CNV.data)), method= "manhattan")

# CNV-based clustering 
Seurat_split.integrated[["cnv_snn"]] <- FindNeighbors(CNV.dist)$snn
for (i in resolution){
  Seurat_split.integrated <- FindClusters(Seurat_split.integrated, assay= "CNV", resolution= i, graph.name= "cnv_snn")
}

# check clustering stability
library(clustree)
clustree(Seurat_split.integrated, prefix="cnv_snn_res.")

# UMAP
# choose resoltuion based on clustree
Seurat_split.integrated <- FindClusters(Seurat_split.integrated, assay= "CNV", resolution= 0.1, graph.name= "cnv_snn")
Seurat_split.integrated[["umap_cnv"]] <- RunUMAP(CNV.dist, assay= "CNV", reduction.key= "cnvUMAP_")
DimPlot(Seurat_split.integrated, reduction= "umap_cnv", label= TRUE)
DimPlot(Seurat_split.integrated, reduction= "umap_cnv", label= TRUE, group.by= "orig.ident")
DimPlot(Seurat_split.integrated, reduction= "umap_cnv", label= TRUE, group.by= "cnv")
DimPlot(Seurat_split.integrated, reduction= "umap_cnv", label= TRUE, split.by= "orig.ident")
DimPlot(Seurat_split.integrated, reduction= "umap_cnv", label= TRUE, group.by= "plate", split.by = "orig.ident")
DimPlot(Seurat_split.integrated, reduction= "umap_cnv", label= TRUE, group.by= "cnv", split.by = "orig.ident")

# stacked bargraph
# percent condition per cluster, change resolution based on the one used for generating UMAP
CNVstacked <- apply(table(Seurat_split.integrated@meta.data$cnv_snn_res.0.1, Seurat_split.integrated@meta.data$orig.ident), 1, function(x){x*100/sum(x,na.rm=T)})
CNVstacked <- as.data.frame(CNVstacked)
CNVstacked$condition <- rownames(CNVstacked)
CNVstacked2 <- melt(CNVstacked)
CNVstacked_bar <- ggplot(CNVstacked2, aes(x= variable, y= value, fill= factor(condition, levels=c("NT", "AG", "R", "RAG")))) + 
  geom_bar(stat= "identity") + theme_classic() +  scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")               
CNVstacked_bar

# percent plate per cluster
CNVstacked_plate <- apply(table(Seurat_split.integrated@meta.data$cnv_snn_res.0.1, Seurat_split.integrated@meta.data$replicate), 1, function(x){x*100/sum(x,na.rm=T)})
CNVstacked_plate <- as.data.frame(CNVstacked_plate)
CNVstacked_plate$condition <- rownames(CNVstacked_plate)
CNVstacked_plate2 <- melt(CNVstacked_plate)
CNVstacked_plate_bar <- ggplot(CNVstacked_plate2, aes(x= variable, y= value, 
                                                      fill= factor(condition, levels=c("NT_plate1", "NT_plate2", "NT_plate3", "AG_plate1", "AG_plate2", "AG_plate3", 
                                                                                       "R_plate1", "R_plate2", "R_plate3", "RAG_plate1", "RAG_plate2", "RAG_plate3")))) + 
  geom_bar(stat= "identity") + theme_classic() +  scale_fill_brewer(palette= "Set3") + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")
CNVstacked_plate_bar

# save CNV clusterID
Seurat_split.integrated[["cnvClusterID"]] <- Idents(Seurat_split.integrated)


# compare RNA vs ADT vs CNV clusters
# based on ADT clustering
umap_rnaClusters <- DimPlot(Seurat_split.integrated, reduction = "umap_tmm", group.by = "rnaClusterID") + NoLegend()
umap_rnaClusters <- umap_rnaClusters + ggtitle("Clustering based on RNA") + theme(plot.title = element_text(hjust = 0.5))
umap_rnaClusters <- LabelClusters(plot = umap_rnaClusters, id = "rnaClusterID", size = 4)

umap_adtClusters <- DimPlot(Seurat_split.integrated, reduction = "umap_tmm", pt.size = 0.5, group.by = "adtClusterID") + NoLegend()
umap_adtClusters <- umap_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
umap_adtClusters <- LabelClusters(plot = umap_adtClusters, id = "adtClusterID", size = 4)

umap_cnvClusters <- DimPlot(Seurat_split.integrated, reduction = "umap_tmm", pt.size = 0.5, group.by = "cnvClusterID") + NoLegend()
umap_cnvClusters <- umap_cnvClusters + ggtitle("Clustering based on CNV signal") + theme(plot.title = element_text(hjust = 0.5))
umap_cnvClusters <- LabelClusters(plot = umap_cnvClusters, id = "cnvClusterID", size = 4)
# Note: for this comparison, RNA, protein, and CNV clustering are visualized on a umap
# generated using the ADT distance matrix.
library(patchwork)
wrap_plots(list(umap_rnaClusters, umap_adtClusters, umap_cnvClusters), ncol = 3)

# based on RNA clustering
umap_rnaClusters2 <- DimPlot(Seurat_split.integrated, reduction = "umap", group.by= "rnaClusterID") + NoLegend()
umap_rnaClusters2 <- umap_rnaClusters2+ ggtitle("Clustering based on RNA") + theme(plot.title = element_text(hjust = 0.5))
umap_rnaClusters2 <- LabelClusters(plot = umap_rnaClusters2, id = "rnaClusterID", size = 4)

umap_adtClusters2 <- DimPlot(Seurat_split.integrated, reduction = "umap", group.by = "adtClusterID", pt.size = 0.5) + NoLegend()
umap_adtClusters2 <- umap_adtClusters2 + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
umap_adtClusters2 <- LabelClusters(plot = umap_adtClusters2, id = "adtClusterID", size = 4)

umap_cnvClusters2 <- DimPlot(Seurat_split.integrated, reduction = "umap", pt.size = 0.5, group.by = "cnvClusterID") + NoLegend()
umap_cnvClusters2 <- umap_cnvClusters2 + ggtitle("Clustering based on CNV signal") + theme(plot.title = element_text(hjust = 0.5))
umap_cnvClusters2 <- LabelClusters(plot = umap_cnvClusters2, id = "cnvClusterID", size = 4)
# Note: for this comparison, both the RNA, protein, and CNV clustering are visualized on a umap
# generated using the RNA clustering
wrap_plots(list(umap_rnaClusters2, umap_adtClusters2, umap_cnvClusters2), ncol = 3)

# based on CNV clustering
umap_rnaClusters3 <- DimPlot(Seurat_split.integrated, reduction = "umap_cnv", group.by= "rnaClusterID") + NoLegend()
umap_rnaClusters3 <- umap_rnaClusters3+ ggtitle("Clustering based on RNA") + theme(plot.title = element_text(hjust = 0.5))
umap_rnaClusters3 <- LabelClusters(plot = umap_rnaClusters3, id = "rnaClusterID", size = 4)

umap_adtClusters3 <- DimPlot(Seurat_split.integrated, reduction = "umap_cnv", group.by = "adtClusterID", pt.size = 0.5) + NoLegend()
umap_adtClusters3 <- umap_adtClusters3 + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
umap_adtClusters3 <- LabelClusters(plot = umap_adtClusters3, id = "adtClusterID", size = 4)

umap_cnvClusters3 <- DimPlot(Seurat_split.integrated, reduction = "umap_cnv", pt.size = 0.5, group.by = "cnvClusterID") + NoLegend()
umap_cnvClusters3 <- umap_cnvClusters3 + ggtitle("Clustering based on CNV signal") + theme(plot.title = element_text(hjust = 0.5))
umap_cnvClusters3 <- LabelClusters(plot = umap_cnvClusters3, id = "cnvClusterID", size = 4)
# Note: for this comparison, both the RNA, protein, and CNV clustering are visualized on a umap
# generated using the RNA clustering
wrap_plots(list(umap_rnaClusters3, umap_adtClusters3, umap_cnvClusters3), ncol = 3)


# WNN clustering
Seurat_split.integrated <- FindMultiModalNeighbors(
  Seurat_split.integrated, reduction.list = list("pca", "pca_tmm", "pca_cnv"), 
  dims.list = list(1:25, 1:25, 1:25), modality.weight.name = "WNN.weight"
)

for (i in resolution){
  Seurat_split.integrated <- FindClusters(Seurat_split.integrated, resolution= i, graph.name= "wsnn", algorithm =1)
}

# check clustering stability
library(clustree)
clustree(Seurat_split.integrated, prefix="wsnn_res.")

# UMAP
Seurat_split.integrated <- RunUMAP(Seurat_split.integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Seurat_split.integrated<- FindClusters(Seurat_split.integrated, graph.name = "wsnn", algorithm = 1, resolution = 0.7, verbose = FALSE)
DimPlot(Seurat_split.integrated, reduction = 'wnn.umap', label = TRUE)
DimPlot(Seurat_split.integrated, reduction = 'wnn.umap', label = TRUE, group.by = "orig.ident")
DimPlot(Seurat_split.integrated, reduction = 'wnn.umap', label = TRUE, split.by = "orig.ident")

# stacked bargraph
# percent condition per cluster, change resolution based on the one used for generating UMAP
WNNstacked <- apply(table(Seurat_split.integrated@meta.data$wsnn_res.0.3, Seurat_split.integrated@meta.data$orig.ident), 1, function(x){x*100/sum(x,na.rm=T)})
WNNstacked <- as.data.frame(WNNstacked)
WNNstacked$condition <- rownames(WNNstacked)
WNNstacked2 <- melt(WNNstacked)
WNNstacked_bar <- ggplot(WNNstacked2, aes(x= variable, y= value, fill= factor(condition, levels=c("NT", "AG", "R", "RAG")))) + 
  geom_bar(stat= "identity") + theme_classic() + scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")                
WNNstacked_bar

# percent plate per cluster
WNNstacked_plate <- apply(table(Seurat_split.integrated@meta.data$wsnn_res.0.3, Seurat_split.integrated@meta.data$replicate), 1, function(x){x*100/sum(x,na.rm=T)})
WNNstacked_plate <- as.data.frame(WNNstacked_plate)
WNNstacked_plate$condition <- rownames(WNNstacked_plate)
WNNstacked_plate2 <- melt(WNNstacked_plate)
WNNstacked_plate_bar <- ggplot(WNNstacked_plate2, aes(x= variable, y= value, 
                                                      fill= factor(condition, levels=c("NT_plate1", "NT_plate2", "NT_plate3", "AG_plate1", "AG_plate2", "AG_plate3", 
                                                                                       "R_plate1", "R_plate2", "R_plate3", "RAG_plate1", "RAG_plate2", "RAG_plate3")))) + 
  geom_bar(stat= "identity") + theme_classic() + scale_fill_brewer(palette= "Set3") + labs(x= "Cluster", y="Percentage of cells", fill= "Conditions")
WNNstacked_plate_bar

# save WNN clusterID
Seurat_split.integrated[["wnnClusterID"]] <- Idents(Seurat_split.integrated)

# export metadata
metadata <- as.data.frame(Seurat_split.integrated[[]])
write.csv(metadata, FILENAME)