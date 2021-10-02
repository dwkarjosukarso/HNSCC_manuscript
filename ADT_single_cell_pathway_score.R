#################################
#       ADT data analysis       #   
#   single-cell pathway scores  #
#################################

#Single cell pathway scores ####
sc_pathway_scores <- function(Seurat_split.integrated){
  
  std_ADT <- Seurat_split.integrated[["TMM"]]@scale.data
  
  #WNT pathway
  WNT_score <- (std_ADT["LRP6-Phospho",] + std_ADT["GSK3B-Phospho",] + std_ADT["Frizzled-3",] + std_ADT["PLC-G2-Phospho",] + std_ADT["PLC-G1 -Phospho",])/5
  names(WNT_score) <- colnames(std_ADT)
  
  #BMP pathway
  BMP_score <- (std_ADT["BMP-2/BMP-4",] + std_ADT["Smad1/5/9-Phospho",] + std_ADT["Smad1/5-Phospho",])/3
  
  #TNFa pathway
  TNFa_score <- (std_ADT["MAPKp38-Phospho-RD",] + std_ADT["NFKBp65-Phospho",] + std_ADT["IKB-A-Phospho",] + std_ADT["CREB-Phospho",] + std_ADT["MAPKp38-Phospho-CST",] + std_ADT["JNK",] + std_ADT["IKB-A",] + std_ADT["NFKBp65/RelA-Phospho",] + std_ADT["MKK3/MKK6-Phospho",] + std_ADT["NFKBp65/RelA",] + std_ADT["MAPK-APK2-Phospho",] + std_ADT["JNK-Phospho",])/12
  
  #TGFb pathway
  TGFb_score <- (std_ADT["Smad3",] + std_ADT["Smad2/3-Phospho",])/2
  
  #Notch pathway
  Notch_score <- (std_ADT["Notch-1", ] + std_ADT["Notch-2",] + std_ADT["Notch-3",] + std_ADT["Notch1-Cleaved",] + std_ADT["Jagged1",])/5
  
  #EGF pathway
  EGF_score <- (std_ADT["mTOR-Phospho",] + std_ADT["EGFR-Phospho-Y1173",] + std_ADT["EGFR-Phospho-Y1045",] + std_ADT["Akt-Phospho",] + std_ADT["Akt1/2/3-Phospho",] + std_ADT["RPS6-Phospho",] + std_ADT["c-Fos-Phospho",] + std_ADT["ERK1/ERK2",] + std_ADT["Src-Phospho",] + std_ADT["c-Jun-Phospho",])/10
  
  #JAK-STAT pathway
  JAKSTAT_score <- (std_ADT["STAT1p91",] + std_ADT["STAT1-Phospho",] + std_ADT["STAT3-Phospho",] + std_ADT["STAT5-Phospho",] + std_ADT["Jak1",] + std_ADT["Jak1-Phospho",])/6
  
  #Adhesion pathway score
  adhesion_score <- (std_ADT["IntegrinB1",] + std_ADT["IntegrinA6",] + std_ADT["FAK",] + std_ADT["FAK-Phospho",] + std_ADT["Src",] + std_ADT["Src-Phospho",])/6
  
  #Differentiation pathway score
  differentiation_score <- (std_ADT["KLF4",] + std_ADT["Notch-1",] + std_ADT["Notch1-Cleaved",] + std_ADT["Notch-2",] + std_ADT["Notch-3",] + std_ADT["Jagged1",] + std_ADT["BMP-2/BMP-4",] + std_ADT["Smad1/5-Phospho",] + std_ADT["Smad1/5/9-Phospho",] + std_ADT["BMPR-II",] + std_ADT["p63/TP73L",] + std_ADT["Transglutaminase1/BC1",] + std_ADT["IntegrinA6",] + std_ADT["IntegrinB1",] + std_ADT["EGFR-Phospho-Y1173",] + std_ADT["EGFR-Phospho-Y1045",])/16
  
  pathway_scores <- data.frame(WNT_score = WNT_score,
                               BMP_score = BMP_score,
                               TNFa_score = TNFa_score,
                               TGFb_score = TGFb_score,
                               Notch_score = Notch_score,
                               EGF_score = EGF_score,
                               JAKSTAT_score = JAKSTAT_score,
                               adhesion_score = adhesion_score,
                               differentiation_score = differentiation_score,
                               condition = Seurat_split.integrated[["orig.ident"]],
                               ADT_cluster = Seurat_split.integrated[["adtClusterID"]],
                               sample = colnames(std_ADT))
  return(pathway_scores)
}