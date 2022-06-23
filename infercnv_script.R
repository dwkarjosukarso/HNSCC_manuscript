# set working directory
setwd(DIR)

# load library
library(infercnv)

# load reference count table
ref <- read.csv(REF_COUNT_TABLE, header= TRUE, row.names = 1)
# adjust colnames 
library(stringr)
colnames(ref) <- substr(colnames(ref), 1, nchar(colnames(ref))-14)

# load observation count table
obs <- read.csv(RNA_COUNT_TABLE, header= T, row.names= 1)

# combine reference and observation count tables
counts <- cbind(ref, obs)


# generate annotation dataframe
annot <- data.frame(x= colnames(counts))
library(dplyr)
annot<- annot %>% mutate(y= ifelse(str_detect(annot$x, "^NT_"), "malignant_NT",
                                   ifelse(str_detect(annot$x, "^AG_"), "malignant_AG",
                                          ifelse(str_detect(annot$x, "^R_"), "malignant_R",
                                                 ifelse(str_detect(annot$x, "^RAG_"), "malignant_RAG",
                                                        ifelse(str_detect(annot$x, "^CTRL"), "normal_KNP",
                                                               ifelse(str_detect(annot$x, "^AG1"), "normal_KNP",
                                                                      ifelse(str_detect(annot$x, "^AG2"), "normal_KNP",
                                                                             ifelse(str_detect(annot$x, "^AG3"), "normal_KNP", "NA")))))))))

rownames(annot) <- annot$x
annot$x <- NULL
names(annot) <- NULL

# load gene order file
gene_order <- read.table(GENE_ORDER, sep="\t")
rownames(gene_order) <- gene_order$V1
gene_order$V1 <- NULL
names(gene_order) <- NULL

# create infercnv object
infercnv <- CreateInfercnvObject(raw_counts_matrix= counts, annotations_file = annot, gene_order_file = gene_order, ref_group_names = "normal_KNP")

# run inferCNV
infercnv_nogroup_subclusters <- infercnv::run(infercnv, cutoff=0.1, out_dir = "new_infercnv_nogroup_subclusters", 
                                              cluster_by_groups = F, num_threads = 5, HMM= T, analysis_mode = "subclusters", tumor_subcluster_partition_method= "random_trees",
                                              denoise= T, BayesMaxPNormal = 0.5, plot_probabilities = T, no_prelim_plot = T,
)



