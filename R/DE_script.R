load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")

source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ML_model.R")


med_endo <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Endothelial', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))

med_epi <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Epithelial', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))

med_fibmes <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Fib_Mesenchymal', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))

med_tcells <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B2',
    celltype_name='T_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))

med_CD8T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Tc_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))

med_CD4T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Th_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))

med_mye <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Myeloid', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))

med_NK <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='NK_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', groups=c('G', 'P'))


# Save RData
save(med_endo, med_epi, med_fibmes, med_tcells, med_CD8T, med_CD4T, med_mye, 
    med_NK, file ='DE_bygroup.RData')