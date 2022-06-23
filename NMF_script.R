# set working directory
setwd(DIR)

# load library
library(Seurat)

#Loading Seurat object
Seurat <- readRDS(SEURAT_OBJECT)

#Gets scaled data
expr <- as.data.frame(GetAssayData(Seurat, assay= "SCT", slot= "data"))
expr <- expr - rowMeans(expr)
expr[expr < 0] <- 0

# calculate NMF scores based on Kinker, et al Method
library(NMF)

# define function
nmf_programs <- function(expr, is.log=F, rank, method="snmf/r", seed=1) {
  
  nmf_programs <- nmf(expr, rank=rank, method=method, seed=seed, .options= 'p1')
  
  nmf_programs_scores <- list(w_basis=basis(nmf_programs), h_coef=t(coef(nmf_programs)))
  
  return(nmf_programs_scores)
}

# run on 5:10 ranges of rank
nmf_program_5 <- nmf_programs(expr, rank=5)
nmf_program_6 <- nmf_programs(expr, rank=6)
nmf_program_7 <- nmf_programs(expr, rank=7)
nmf_program_8 <- nmf_programs(expr, rank=8)
nmf_program_9 <- nmf_programs(expr, rank=9)
nmf_program_10 <- nmf_programs(expr, rank=10)

# save output
saveRDS(nmf_program_5, "nmf_5.rds")
saveRDS(nmf_program_6, "nmf_6.rds")
saveRDS(nmf_program_7, "nmf_7.rds")
saveRDS(nmf_program_8, "nmf_8.rds")
saveRDS(nmf_program_9, "nmf_9.rds")
saveRDS(nmf_program_10, "nmf_10.rds")