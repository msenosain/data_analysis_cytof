
# Load functions
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/51_supervised_analysis_viz.R")
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


###############################################################################
# GENERAL UMAP
###############################################################################
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

umap_f <- function(dt, sample_by_patient = TRUE, sample_k = 1000, col_id, 
    subset = FALSE, subset_col, subset_nm, seed = 65){
    if(subset){
        k <- which(dt[,subset_col] == subset_nm)
        dt <- dt[k,]
    }

    set.seed(seed)
    #print(as.matrix(colnames(annot_df)))
    if(sample_by_patient){
        # sampling sample_k cells per patient
        smpl <- dt[1,][-1,]
        pts = unique(dt$pt_ID)
        replacement = FALSE
        for(i in pts){
            len <- length(which(dt$pt_ID %in% i))
            if(len<sample_k){
                replacement = TRUE
            }
            x <- dplyr::sample_n(dt[which(dt$pt_ID %in% i),], size=sample_k, replace=replacement)
            smpl <- rbind(smpl, x)
        }        
    } else {
        smpl <- dplyr::sample_n(dt, size=sample_k, replace=F)
    }

    smpl <- cbind(denoisingCTF::t_asinh(smpl[,1:77]), smpl[,78:89])

    #UMAP
    umap_smp <- umap::umap(smpl[,col_id])
    umap_smp <- as.data.frame(cbind(smpl,
                                        UMAP1 = umap_smp$layout[,1],
                                        UMAP2 = umap_smp$layout[,2]))
    umap_smp$n_op1 <- as.factor(umap_smp$n_op1)
    umap_smp$n_op2 <- as.factor(umap_smp$n_op2)
    umap_smp$subtype_A <- as.factor(umap_smp$subtype_A)
    umap_smp$subtype_B <- as.factor(umap_smp$subtype_B)
    umap_smp$clusters_A <- as.factor(umap_smp$clusters_A)
    umap_smp$clusters_B <- as.factor(umap_smp$clusters_B)

    return(umap_smp)
}


###############################################################################
# ALL CELLS UMAP
###############################################################################

ft_cols <- col_id <- c(15,20,29:31,34,37,40,44,46,47)
umap_smp <- umap_f(annot_df, sample_k = 1000, sample_by_patient = TRUE, 
    col_id = ft_cols, subset = FALSE, seed = 65)

save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_allcelltypes.RData")

###############################################################################
# Epi UMAP
###############################################################################
ft_cols <- c(15,18,19,21,22,23,26,28,31,33,35,37,38,40,41,50,51)
umap_smp <- umap_f(annot_df, sample_k = 25000, sample_by_patient = FALSE,
    col_id = ft_cols, subset = TRUE, subset_col='subtype_B', 
    subset_nm='Epithelial', seed = 50)
save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_epi.RData")

###############################################################################
# Endo UMAP
###############################################################################
ft_cols <- c(18, 19, 21, 22, 26, 28, 31, 33, 35, 38, 41, 42, 50, 51)
umap_smp <- umap_f(annot_df, sample_k = 25000, sample_by_patient = FALSE, 
    col_id = ft_cols, subset = TRUE, 
    subset_col='subtype_B', subset_nm='Endothelial', seed = 50)
save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_endo.RData")

###############################################################################
# FMes UMAP
###############################################################################
ft_cols <- c(18, 19, 21, 22, 26, 27, 28, 33, 35, 38, 41, 42, 43, 48, 50, 51)
umap_smp <- umap_f(annot_df, sample_k = 25000, sample_by_patient = FALSE, 
    col_id = ft_cols, subset = TRUE, 
    subset_col='subtype_B', subset_nm='Fib_Mesenchymal', seed = 50)
save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_fmes.RData")

###############################################################################
# CD4 UMAP
###############################################################################
ft_cols <- c(18:22, 26, 28, 31, 35, 42, 50)
umap_smp <- umap_f(annot_df, sample_k = 25000, sample_by_patient = FALSE, 
    col_id = ft_cols, subset = TRUE, 
    subset_col='subtype_B', subset_nm='Th_cells', seed = 50)
save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_cd4t.RData")

###############################################################################
# CD8 UMAP
###############################################################################
ft_cols <- c(18:22, 26, 28, 31, 35, 42, 50)
umap_smp <- umap_f(annot_df, sample_k = 25000, sample_by_patient = FALSE, 
    col_id = ft_cols, subset = TRUE, 
    subset_col='subtype_B', subset_nm='Tc_cells', seed = 50)
save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_cd8t.RData")

###############################################################################
# DNT UMAP
###############################################################################
ft_cols <- c(18:22, 26, 28, 31, 35, 42, 50)
umap_smp <- umap_f(annot_df, sample_k = 25000, sample_by_patient = FALSE,
    col_id = ft_cols, subset = TRUE, 
    subset_col='subtype_B', subset_nm='DNT_cells', seed = 50)
save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_dnt.RData")

###############################################################################
# Myeloid UMAP
###############################################################################
ft_cols <- c(18:22, 26, 28, 31, 35, 42, 50)
umap_smp <- umap_f(annot_df, sample_k = 25000, sample_by_patient = FALSE, 
    col_id = ft_cols, subset = TRUE, 
    subset_col='subtype_B', subset_nm='Myeloid', seed = 50)
save(umap_smp, file = "/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/supervised/UMAP_mye.RData")
