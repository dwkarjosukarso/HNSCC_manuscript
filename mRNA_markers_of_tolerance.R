##########################################################
#       Identification of the markers of tolerance       # 
#     using Random Forest & validation with stepwise     #  
#              regression-based linear model             #  
##########################################################

library(Seurat)
library(readxl)
library(writexl)
library(randomForest)
library(rfPermute)
library(caret)
library(MASS)

#Loading the robust expression programs
robust_expression_programs <- lapply(as.list(read_xlsx(path = "Robust_programs_KM.xlsx", col_names = TRUE)), function(x) as.vector(na.omit(x)))

#DK114
#Loading the Seurat object
DK114 <- readRDS("DK114_Seurat_split.integrated.RDS") 

#Retrieving the data slot as a matrix 
DK114_SCT_data <- as.matrix(GetAssayData(object = DK114, assay = "SCT", slot = "data"))

#Transformation of the matrix into a vector
DK114_SCT_data_vector <- as.vector(DK114_SCT_data)

#Select the expression level for which 5% of the genes have their average expression among all cells 
#at that level or higher. This expression level corresponds to the 95th percentile of this distribution
DK114_top_level <- quantile(sort(DK114_SCT_data_vector), probs = c(0.95))
DK114_genes_means <- rowMeans(DK114_SCT_data)

#Filter the SCT data slot for the genes with expression level => than the 95th percentile
DK114_high_expr_genes_means <- DK114_genes_means[DK114_genes_means >= DK114_top_level]
DK114_high_expr_genes <- names(DK114_high_expr_genes_means)

#Saving data frame based on the expression of these highly expressed genes only
DK114_all_high_expr_mRNAs <- as.data.frame(DK114[["SCT"]]@data[DK114_high_expr_genes,])
DK114_all_high_expr_mRNAs_tolerance <- as.data.frame(t(DK114_all_high_expr_mRNAs))

#Filtering the robust expression programs for the genes filtered by the expression level 
DK114_high_expr_genes_robust_expression_programs <- lapply(robust_expression_programs, function(x) intersect(DK114_high_expr_genes, x))

#DK118
#Loading the Seurat object
DK118 <- readRDS("DK118_Seurat_split.integrated.RDS") 

#Retrieving the data slot as a matrix 
DK118_SCT_data <- as.matrix(GetAssayData(object = DK118, assay = "SCT", slot = "data"))

#Transformation of the matrix into a vector
DK118_SCT_data_vector <- as.vector(DK118_SCT_data)

#Since we want to select the expression level for which 5% of the genes have their average expression among all cells at that level
#Thus, this expression level corresponds to the third quartile of this distribution
DK118_top_level <- quantile(sort(DK118_SCT_data_vector), probs = c(0.95))
DK118_genes_means <- rowMeans(DK118_SCT_data)

#Filter the SCT data slot for the genes with expression level greater or equal than the 95th percentile
DK118_high_expr_genes_means <- DK118_genes_means[DK118_genes_means >= DK118_top_level]
DK118_high_expr_genes <- names(DK118_high_expr_genes_means)

#Saving data frame based on the expression of these highly expressed genes only
DK118_all_high_expr_mRNAs <- as.data.frame(DK118[["SCT"]]@data[DK118_high_expr_genes,])
DK118_all_high_expr_mRNAs_tolerance <- as.data.frame(t(DK118_all_high_expr_mRNAs))

#Filtering the robust expression programs for the genes filtered by the expression level such that 5% of the genes have >= that expression level
DK118_high_expr_genes_robust_expression_programs <- lapply(robust_expression_programs, function(x) intersect(DK118_high_expr_genes, x))

#Taking the union of the mRNAs selected in the two patients' datasets
#Keeping the info on what genes are patient-specific
#First we check that for both patients we have the same programs as elements in the respective lists of mRNAs

if (length(DK114_high_expr_genes_robust_expression_programs) == length(DK118_high_expr_genes_robust_expression_programs) && names(DK114_high_expr_genes_robust_expression_programs) == names(DK118_high_expr_genes_robust_expression_programs)){
  mRNA_tolerance_markers <- list()
  DK114_specific <- list()
  DK118_specific <- list()
  for (program in 1:length(DK114_high_expr_genes_robust_expression_programs)){
    mRNA_tolerance_markers[[program]] <- unique(union(DK114_high_expr_genes_robust_expression_programs[[program]], DK118_high_expr_genes_robust_expression_programs[[program]]))
    DK114_specific[[program]]  <- setdiff(DK114_high_expr_genes_robust_expression_programs[[program]], DK118_high_expr_genes_robust_expression_programs[[program]])
    DK118_specific[[program]]  <- setdiff(DK118_high_expr_genes_robust_expression_programs[[program]], DK114_high_expr_genes_robust_expression_programs[[program]])
  }
}

#Assigning names to the programs
#In this way, we know in which programs - from either of the patients- the mRNAs are included
names(mRNA_tolerance_markers) <- names(DK114_high_expr_genes_robust_expression_programs)

#Creating a data frame out of each list, and saving it as excel file ready for KM plot
for (program in names(mRNA_tolerance_markers)){
  onecolumn_dataframe <- data.frame(mRNA_tolerance_markers[[program]], stringsAsFactors = TRUE)
  write_xlsx(onecolumn_dataframe, path = paste0("mRNA_tolerance_markers_union_",gsub("/", "_", program), ".xlsx"))
}

#Loading the tolerance score from both the patients
DK114_tolerance_score <- read.delim(file = "SJG026_DK114_tolerance_score.csv", sep = ",", row.names = 1)
DK118_tolerance_score <- read.delim(file = "SJG017_DK118_tolerance_score.csv", sep = ",", row.names = 1)

#Retrieving mRNA expression for both DK114 and DK118 of the mRNAs that can be candidate for tolerance markers from robust expression programs
DK114_candidate_mRNAs_expression <- as.data.frame(DK114[["SCT"]]@data[unique(unlist(mRNA_tolerance_markers)),])
DK118_candidate_mRNAs_expression <- as.data.frame(DK118[["SCT"]]@data[unique(unlist(mRNA_tolerance_markers)),])

#Creating dataframes for random forest
DK114_tol_mRNA <- as.data.frame(t(DK114_candidate_mRNAs_expression))
DK114_tol_mRNA$Tolerance_score <- DK114_tolerance_score$tolerance_score_DK114
DK118_tol_mRNA <- as.data.frame(t(DK118_candidate_mRNAs_expression))
DK118_tol_mRNA$Tolerance_score <- DK118_tolerance_score$tolerance_score_DK118
rownames(DK114_tol_mRNA) <- unlist(lapply(rownames(DK114_tol_mRNA), function(x) paste0("DK114_", x)))
rownames(DK118_tol_mRNA) <- unlist(lapply(rownames(DK118_tol_mRNA), function(x) paste0("DK118_", x)))
all_data_tol_mRNA <- rbind(DK114_tol_mRNA, DK118_tol_mRNA)
colnames(all_data_tol_mRNA) <- make.names(colnames(all_data_tol_mRNA))

#Finding best mtry for random forest based on the highly expressed genes present in the robust expression programs

#1. With tuneRF we find the best range
set.seed(1)
tuneRF_res = lapply(1:10,function(i){
  tr = tuneRF(x = all_data_tol_mRNA[, -ncol(all_data_tol_mRNA)], 
              y = all_data_tol_mRNA[, ncol(all_data_tol_mRNA)],
              mtryStart=2, step=0.9, ntreeTry = 100, trace = TRUE, improve = 1e-5)
  tr = data.frame(tr)
  tr$RMSE = sqrt(tr[,2])
  tr
})

tuneRF_res = do.call(rbind,tuneRF_res)

#2. The min and max mtry identified will be our range for 10-fold cross validation to identify the best mtry
set.seed(1)
control <- trainControl(method = "cv", number = 10, returnResamp = "all", allowParallel = TRUE)
tunegrid <- expand.grid(.mtry = c(min(tuneRF_res$mtry):max(tuneRF_res$mtry)+5))
caret_res <- train(Tolerance_score ~., data = all_data_tol_mRNA, method = "rf", metric = "RMSE", 
                   tuneGrid = tunegrid, ntree = 100, trControl=control)

#The best mtry is
caret_res$bestTune

#Random forest with permutations: in this way, we obtain a p-value of the importance values associated to each variable.
set.seed(1)
tolerance.rf.permut <- rfPermute(Tolerance_score~., 
                                 data = all_data_tol_mRNA,
                                 nrep = 100, 
                                 importance = TRUE, num.cores = 32, 
                                 mtry = caret_res$bestTune$mtry, ntree = 1000)

rp.importance(tolerance.rf.permut)
varImpPlot(tolerance.rf.permut, type = 1, n.var = nrow(rp.importance(tolerance.rf.permut)))

#Identification of the mRNA markers of tolerance
#1. Selecting mRNAs with p-value of importance < 0.05
mRNA_imp <- as.data.frame(rp.importance(tolerance.rf.permut))
write.csv(mRNA_imp, file = "mRNA_importance_RF_union_results.csv")
#Retaining just the %IncMSE column and p-val, which correspond to the two first columns
mRNA_imp <- mRNA_imp[,c(1,2)]
#Retaining only most significant mRNAs
significant_mRNAs <- mRNA_imp[mRNA_imp$`%IncMSE.pval` < 0.05,]

#2. Select all the variables whose importance values is > than the mean of the importance values
mean_imp <- mean(significant_mRNAs$`%IncMSE`)
significant_top_mRNAs <- significant_mRNAs[significant_mRNAs$`%IncMSE` > mean_imp,] 
write.csv(significant_top_mRNAs, file = "significant_and_important_mRNAs_RF_union.csv")

#3. Linear model based on the same set of candidate mRNAs and stepwise regression based on AIC with both the patients' data
set.seed(1)
tolerance_lm <- lm(Tolerance_score ~ ., data = all_data_tol_mRNA)
tolerance_lm_step <- stepAIC(tolerance_lm, direction = "both", trace = 0)

#10-fold cross-validation of the AIC-based stepwise regression lm
set.seed(1)
#Specify the cross-validation method
ctrl_10 <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
#Fit a regression model and use k-fold CV to evaluate performance
tolerance_model_10 <- train(tolerance_lm_step$call[[2]], 
                            data = all_data_tol_mRNA, 
                            method = "lm", trControl = ctrl_10)

#View summary
print(tolerance_model_10)

#Now we compare the predictors that have been retained and whether they satisfied the significance criterion in the random forest model
#Total candidates mRNAs 
length(unique(unlist(mRNA_tolerance_markers)))

#The coefficients that are retained are (the first one is the intercept):
length(tolerance_lm_step$coefficients)-1

#Whereas the coefficients with high significance in the random forest model are
nrow(significant_mRNAs)

#What is the intersection among these two lists of candidate markers of tolerance?
intersect(rownames(significant_mRNAs), names(tolerance_lm_step$coefficients)[-1])
length(intersect(rownames(significant_mRNAs), names(tolerance_lm_step$coefficients)[-1]))

#This means that the % of genes in the stepwise regression linear model are shared with the significant mRNAs according to the random forest model:
length(intersect(rownames(significant_mRNAs), names(tolerance_lm_step$coefficients)[-1]))/length(names(tolerance_lm_step$coefficients)[-1])*100

#On the other hand, the % of genes in the random forest model that are shared with the stepwise regression linear model are:
length(intersect(rownames(significant_mRNAs), names(tolerance_lm_step$coefficients)[-1]))/length(rownames(significant_mRNAs))*100

#The differences in between the two sets:
setdiff(rownames(significant_mRNAs), names(tolerance_lm_step$coefficients)[-1])
setdiff(names(tolerance_lm_step$coefficients)[-1], rownames(significant_mRNAs))

#Using the TOP mRNAs selected with RF - satisfying also the importance criterion
intersect(rownames(significant_top_mRNAs), names(tolerance_lm_step$coefficients)[-1])
length(intersect(rownames(significant_top_mRNAs), names(tolerance_lm_step$coefficients)[-1]))

#This means that the % of genes in the stepwise regression linear model are shared with the significant mRNAs according to the random forest model:
length(intersect(rownames(significant_top_mRNAs), names(tolerance_lm_step$coefficients)[-1]))/length(names(tolerance_lm_step$coefficients)[-1])*100

#On the other hand, the % of genes in the random forest model that are shared with the stepwise regression linear model are:
length(intersect(rownames(significant_top_mRNAs), names(tolerance_lm_step$coefficients)[-1]))/length(rownames(significant_top_mRNAs))*100

#The differences in between the two sets:
setdiff(rownames(significant_mRNAs), names(tolerance_lm_step$coefficients)[-1])
setdiff(names(tolerance_lm_step$coefficients)[-1], rownames(significant_mRNAs))
#Summary of the two models:
tolerance.rf.permut
summary(tolerance_lm_step)

#Running patient-specific models
DK114_specific_candidate_markers_mRNAs_expression <- as.data.frame(DK114[["SCT"]]@data[unique(unlist(DK114_high_expr_genes_robust_expression_programs)),])
DK118_specific_candidate_markers_mRNAs_expression <- as.data.frame(DK118[["SCT"]]@data[unique(unlist(DK118_high_expr_genes_robust_expression_programs)),])
DK114_specific_tol_mRNA <- t(DK114_specific_candidate_markers_mRNAs_expression)
DK118_specific_tol_mRNA <- t(DK118_specific_candidate_markers_mRNAs_expression)

DK114_specific_tol_mRNA$Tolerance_score <- DK114_tolerance_score$tolerance_score_DK114
DK118_specific_tol_mRNA$Tolerance_score <- DK118_tolerance_score$tolerance_score_DK118

#In order to find the patient-specific mRNAs involved in tolerance, we run patient-specific RF models
#1. With tuneRF we find the best range
#DK114
colnames(DK114_tol_mRNA) <- make.names(colnames(DK114_tol_mRNA))
set.seed(1)
DK114_tuneRF_res = lapply(1:10,function(i){
  tr = tuneRF(x = DK114_tol_mRNA[, -ncol(DK114_tol_mRNA)], 
              y = DK114_tol_mRNA[, ncol(DK114_tol_mRNA)],
              mtryStart=2, step=0.9, ntreeTry = 100, trace = TRUE, improve = 1e-5)
  tr = data.frame(tr)
  tr$RMSE = sqrt(tr[,2])
  tr
})

DK114_tuneRF_res = do.call(rbind,DK114_tuneRF_res)

#2. The min and max mtry identified will be our range for 10-fold cross validation to identify the best mtry
set.seed(1)
DK114_control <- trainControl(method = "cv", number = 10, returnResamp = "all", allowParallel = TRUE)
DK114_tunegrid <- expand.grid(.mtry = c(min(DK114_tuneRF_res$mtry):max(DK114_tuneRF_res$mtry)+5))
DK114_caret_res <- train(Tolerance_score ~., data = DK114_tol_mRNA, method = "rf", metric = "RMSE", 
                         tuneGrid = DK114_tunegrid, ntree = 100, trControl=DK114_control)

#The best mtry is
DK114_caret_res$bestTune

#Random forest with permutations: in this way, we obtain a p-value of the importance values associated to each variable.
set.seed(1)
DK114_tolerance.rf.permut <- rfPermute(Tolerance_score~., 
                                       data = DK114_tol_mRNA,
                                       nrep = 100, 
                                       importance = TRUE, num.cores = 32, 
                                       mtry = DK114_caret_res$bestTune$mtry, ntree = 1000)

#Identification of the mRNA markers of tolerance
#1. Selecting mRNAs with p-value of importance < 0.05
DK114_mRNA_imp <- as.data.frame(rp.importance(DK114_tolerance.rf.permut))
write.csv(DK114_mRNA_imp, file = "DK114_mRNA_importance_RF_union_results.csv")
#Retaining just the %IncMSE column and p-val, which correspond to the two first columns
DK114_mRNA_imp <- DK114_mRNA_imp[,c(1,2)]
#Retaining only most significant mRNAs
DK114_significant_mRNAs <- DK114_mRNA_imp[DK114_mRNA_imp$`%IncMSE.pval` < 0.05,]

#2. Select all the variables whose importance values is > than the mean of the importance values
DK114_mean_imp <- mean(DK114_significant_mRNAs$`%IncMSE`)
DK114_significant_top_mRNAs <- DK114_significant_mRNAs[DK114_significant_mRNAs$`%IncMSE` > DK114_mean_imp,] 
write.csv(DK114_significant_top_mRNAs, file = "DK114_significant_and_important_mRNAs_RF_union.csv")

#DK118
#1. With tuneRF we find the best range
colnames(DK118_tol_mRNA) <- make.names(colnames(DK118_tol_mRNA))
set.seed(1)
DK118_tuneRF_res = lapply(1:10,function(i){
  tr = tuneRF(x = DK118_tol_mRNA[, -ncol(DK118_tol_mRNA)], 
              y = DK118_tol_mRNA[, ncol(DK118_tol_mRNA)],
              mtryStart=2, step=0.9, ntreeTry = 100, trace = TRUE, improve = 1e-5)
  tr = data.frame(tr)
  tr$RMSE = sqrt(tr[,2])
  tr
})

DK118_tuneRF_res = do.call(rbind,DK118_tuneRF_res)

#2. The min and max mtry identified will be our range for 10-fold cross validation to identify the best mtry
set.seed(1)
DK118_control <- trainControl(method = "cv", number = 10, returnResamp = "all", allowParallel = TRUE)
DK118_tunegrid <- expand.grid(.mtry = c(min(DK118_tuneRF_res$mtry):max(DK118_tuneRF_res$mtry)+5))
DK118_caret_res <- train(Tolerance_score ~., data = DK118_tol_mRNA, method = "rf", metric = "RMSE", 
                         tuneGrid = DK118_tunegrid, ntree = 100, trControl=DK118_control)

#The best mtry is
DK118_caret_res$bestTune

#Random forest with permutations: in this way, we obtain a p-value of the importance values associated to each variable.
set.seed(1)
DK118_tolerance.rf.permut <- rfPermute(Tolerance_score~., 
                                       data = DK118_tol_mRNA,
                                       nrep = 100, 
                                       importance = TRUE, num.cores = 32, 
                                       mtry = DK118_caret_res$bestTune$mtry, ntree = 1000)

#Identification of the mRNA markers of tolerance
#1. Selecting mRNAs with p-value of importance < 0.05
DK118_mRNA_imp <- as.data.frame(rp.importance(DK118_tolerance.rf.permut))
write.csv(DK118_mRNA_imp, file = "DK118_mRNA_importance_RF_union_results.csv")
#Retaining just the %IncMSE column and p-val, which correspond to the two first columns
DK118_mRNA_imp <- DK118_mRNA_imp[,c(1,2)]
#Retaining only most significant mRNAs
DK118_significant_mRNAs <- DK118_mRNA_imp[DK118_mRNA_imp$`%IncMSE.pval` < 0.05,]

#2. Select all the variables whose importance values is > than the mean of the importance values
DK118_mean_imp <- mean(DK118_significant_mRNAs$`%IncMSE`)
DK118_significant_top_mRNAs <- DK118_significant_mRNAs[DK118_significant_mRNAs$`%IncMSE` > DK118_mean_imp,] 
write.csv(DK118_significant_top_mRNAs, file = "DK118_significant_and_important_mRNAs_RF_union.csv")

save.image(file = "mRNAs_markers_of_tolerance.RData")
