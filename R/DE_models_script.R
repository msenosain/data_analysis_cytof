load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")

source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ML_model.R")


df <- equal_sampling(annot_df, group_col = 'CANARY', ptID_col = 'pt_ID', 
    groups = c('G', 'P'))

# 15 18 19 21 22 23 24 25 26 27 28 29 31 33 35 37 38 40 41 42 43 45 48 50 51

# Endothelial cells
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'cell_type_B', celltype_name = 'Endothelial',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'Endothelial', allowParallel = TRUE, workers = 12)

# Fib_Mesenchymal
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'cell_type_B', celltype_name = 'Fib_Mesenchymal',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'Fib_Mesenchymal', allowParallel = TRUE, workers = 12)

# T cells (one)
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'subtype_B2', celltype_name = 'T_cells',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'T_cells', allowParallel = TRUE, workers = 12)

# T cells CD8+
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'subtype_B', celltype_name = 'Tc_cells',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'Tc_cells', allowParallel = TRUE, workers = 12)

# T cells CD4+
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'subtype_B', celltype_name = 'Th_cells',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'Th_cells', allowParallel = TRUE, workers = 12)

# Myeloid
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'subtype_B', celltype_name = 'Myeloid',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'Myeloid', allowParallel = TRUE, workers = 12)

# NK
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'subtype_B', celltype_name = 'NK_cells',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'NK_cells', allowParallel = TRUE, workers = 12)

# Epithelial
ML_model(df, balance = T,train_size = 0.75, subset_celltype = T, 
    celltype_col = 'cell_type_B', celltype_name = 'Epithelial',
    alg = 'RF', class_col = 'CANARY', seed = 45, ask_features = F,
    ft_idxs = c(15, 17:31, 33:35, 37:48, 50:51),
    label = 'Epithelial', allowParallel = TRUE, workers = 12)





# col1 <- brewer.pal(12, "Set3")
# heatmap.2(as.matrix(prcnt_by_pt[x,]), col=greenred(75),
#               trace="none",
#               keysize=1,
#               margins=c(8,6),
#               scale="none",
#               dendrogram="none",
#               Colv = FALSE,
#               Rowv = FALSE,
#               cexRow=0.5 + 1/log10(dim(prcnt_by_pt)[1]),
#               cexCol=1.25,
#               main="Genes grouped by categories",
#               RowSideColors=col1[as.numeric(ref2$CANARY[x])]
# )










