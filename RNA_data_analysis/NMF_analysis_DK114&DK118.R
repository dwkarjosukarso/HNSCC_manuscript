#########################################
#          RNA data analysis            #
#        NMF robust expression          #
#     programs on DK114 & DK118         #
#             datasets                  #
#########################################

#This function extracts the top 50 genes by NMF basis score
#Which correspond to the metagenes.
source("NMF_metagenes_function.R")

#Kinker et al. robust NMF programs function.
source("robust_nmf_programs.R")

#Metagenes in each dataset for the ranks 5 to 10
#Each object is a list whose elements correspond to the ranks
SJG026_DK114_metagenes <- NMFmetagenesSeurat(path_nmf_rds = "./NMF_SJG026_DK114/", dataset_patient = "SJG026_DK114")
SJG017_DK118_metagenes <- NMFmetagenesSeurat(path_nmf_rds = "./NMF_SJG017_DK118/", dataset_patient = "SJG017_DK118")

SJG026_DK114_metagenes <- do.call(cbind, SJG026_DK114_metagenes)
SJG017_DK118_metagenes <- do.call(cbind, SJG017_DK118_metagenes)

colnames(SJG026_DK114_metagenes) <- unlist(lapply(colnames(SJG026_DK114_metagenes), function(x) paste0("SJG026_DK114_", x)))
colnames(SJG017_DK118_metagenes) <- unlist(lapply(colnames(SJG017_DK118_metagenes), function(x) paste0("SJG017_DK118_", x)))

#Assigning each matrix to an element in the list standing for the dataset of origin
nmf_programs<- list(SJG026_DK114_metagenes,
                    SJG017_DK118_metagenes)
                   
names(nmf_programs) <- c("SJG026_DK114",
                         "SJG017_DK118")

#Robust programs
RB <- robust_nmf_programs(nmf_programs)

for (robust_program in RB[grep("DK114", RB)]){
  print(robust_program)
  print(nmf_programs$SJG026_DK114[, robust_program])
  cat("\n")
  write.csv(nmf_programs$SJG026_DK114[, robust_program], file = paste0(robust_program, '.csv'))
}

for (robust_program in RB[grep("DK118", RB)]){
  print(robust_program)
  print(nmf_programs$SJG017_DK118[, robust_program])
  cat("\n")
  write.csv(nmf_programs$SJG017_DK118[, robust_program], file = paste0(robust_program, '.csv'))
}
