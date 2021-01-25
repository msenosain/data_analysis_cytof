# Load functions
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/50_supervised_analysis.R")
# Load data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellclusters.RData")
ref <- ref[-42,] # sample 13376 has been collected twice and will be merged, only one ref is needed
CDE <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2021JAN13_SA_MF.csv')
CDE <- CDE[match(ref$pt_ID, CDE$pt_ID),] # Select only pts with CyTOF data

# Add new class and SILA to annot_df
pts = unique(annot_df$pt_ID)
annot_df['SILA'] <- NA
annot_df['n_op1'] <- NA
annot_df['n_op2'] <- NA
for (i in pts) {
    k <- which(annot_df$pt_ID == i)
    annot_df[k, 'SILA'] <- CDE[which(CDE$pt_ID == i), 'SILA']
    annot_df[k, 'n_op1'] <- CDE[which(CDE$pt_ID == i), 'n_op1']
    annot_df[k, 'n_op2'] <- CDE[which(CDE$pt_ID == i), 'n_op2']
}
annot_df$n_op1 <- as.factor(annot_df$n_op1)
annot_df$n_op2 <- as.factor(annot_df$n_op2)
annot_df$subtype_A <- as.factor(annot_df$subtype_A)
annot_df$subtype_B <- as.factor(annot_df$subtype_B)
annot_df$clusters_A <- as.factor(annot_df$clusters_A)
annot_df$clusters_B <- as.factor(annot_df$clusters_B)


# # Equal sampling by n_op1
# seed = 65
# set.seed(seed)
# smpl <- annot_df[1,][-1,]
# k = 35000
# o1 = unique(annot_df$n_op1)
# for(i in o1){
#   x <- dplyr::sample_n(annot_df[which(annot_df$n_op1 %in% i),], size=k, replace=F)
#   smpl <- rbind(smpl, x)
# }
# # Equal sampling by n_op2
# seed = 65
# set.seed(seed)
# smpl <- annot_df[1,][-1,]
# k = 24000
# o2 = unique(annot_df$n_op2)
# for(i in o2){
#   x <- dplyr::sample_n(annot_df[which(annot_df$n_op2 %in% i),], size=k, replace=F)
#   smpl <- rbind(smpl, x)
# }


# Entire data set UMAP
seed = 65
set.seed(seed)
#print(as.matrix(colnames(annot_df)))
# sampling k cells per patient
smpl <- annot_df[1,][-1,]
k = 1000
pts = unique(annot_df$pt_ID)
for(i in pts){
  x <- dplyr::sample_n(annot_df[which(annot_df$pt_ID %in% i),], size=k, replace=F)
  smpl <- rbind(smpl, x)
}


smpl <- cbind(denoisingCTF::t_asinh(smpl[,1:77]), smpl[,78:89])
col_id <- c(15,20,29:31,34,37,40,44,46,47)

#UMAP
umap_smp <- umap::umap(smpl[,col_id])
umap_smp <- as.data.frame(cbind(smpl,
                                    UMAP1 = umap_smp$layout[,1],
                                    UMAP2 = umap_smp$layout[,2]))



save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_allcelltypes.RData")