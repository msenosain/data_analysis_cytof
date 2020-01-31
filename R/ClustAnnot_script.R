#########################################################
# Data without CD90
#########################################################

# Load data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/woCD90.RData")

# Source code
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ClustAnnot_functions.R")

# Clustering
dt_cl <- clustering(big_df, n_clusters = 10, iterations = 200, seed = 45) 
# 15 20 29 31 37 40

# Clustering Evaluation
hm_data <- ClusterEval_data(dt_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

tsne_data <- ClusterEval_data(dt_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(dt_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

table(dt_cl$cluster)/nrow(dt_cl)*100

# Cluster annotation
epithelial <- c(6,7,9,10)
endothelial <- c(8)
mesenchymal <- c(4)
immune <- c(1,2,3)
nothing <- c(5)
#unknown <- c(6)

ct_ls <- list('Epithelial'= epithelial, 'Endothelial'=endothelial, 
                   'Mesenchymal'=mesenchymal, 'Immune'=immune, 'Nothing'=nothing)

annot_df <- ClusterAnnotation(data = big_df, df_cluster = dt_cl, 
    ls_annotation = ct_ls, annotation_col = 'cell_type', cl_delete = T, 
    cl_delete_name = 'nothing')

prcnt_by_pt <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'cell_type')

save(hm_data, tsne_data, umap_data, annot_df, prcnt_by_pt, file = 'majorcelltypes.RData')


# Do this for both, wCD90 woCD90 separately and then merge.
# Subset cell types for further subtyping

#########################################################
# Data with CD90
#########################################################

# Load data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/withCD90.RData")

# Source code
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ClustAnnot_functions.R")


# Clustering
dt_cl <- clustering(big_df, n_clusters = 12, iterations = 200, seed = 45)
# 15 20 29 31 37 40 49

# Clustering Evaluation
hm_data <- ClusterEval_data(dt_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

tsne_data <- ClusterEval_data(dt_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(dt_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

table(dt_cl$cluster)/nrow(dt_cl)*100

# Cluster annotation
epithelial <- c(1,2,4)
endothelial <- c(7)
fibroblasts <- c(9)
mesenchymal <- c(11)
immune <- c(3,5,6,8,12)
nothing <- c(10)

ct_ls <- list('Epithelial'= epithelial, 'Endothelial'=endothelial, 
                   'Fibroblasts' = fibroblasts, 'Mesenchymal'=mesenchymal, 
                   'Immune'=immune, 'Nothing'=nothing)

annot_df <- ClusterAnnotation(data = big_df, df_cluster = dt_cl, 
    ls_annotation = ct_ls, annotation_col = 'cell_type', cl_delete = T, 
    cl_delete_name = 'Nothing')

prcnt_by_pt <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'cell_type')

save(hm_data, tsne_data, umap_data, annot_df, prcnt_by_pt, file = 'majorcelltypes.RData')


#########################################################
# Combining both w and wo CD90
#########################################################

#save(annot_df, ref, annot_df_wo, ref_wo, file = 'majorcelltypes.RData')

# Source code
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ClustAnnot_functions.R")

annot_df <- change_colname(annot_df, 'cell_type', 'cell_type_A')
annot_df_wo <- change_colname(annot_df_wo, 'cell_type', 'cell_type_A')

annot_df['cell_type_B'] <- annot_df$cell_type_A
annot_df_wo['cell_type_B'] <- annot_df_wo$cell_type_A

annot_df <- edit_names(annot_df, col_name = 'cell_type_B', var_oldname = 'Fibroblasts',
                var_newname = 'Fib_Mesenchymal')
annot_df <- edit_names(annot_df, col_name = 'cell_type_B', var_oldname = 'Mesenchymal',
                var_newname = 'Fib_Mesenchymal')

annot_df_wo <- edit_names(annot_df_wo, col_name = 'cell_type_B', var_oldname = 'Mesenchymal',
                var_newname = 'Fib_Mesenchymal')

annot_df <- rbind(annot_df_wo, annot_df)
ref <- rbind(ref_wo, ref)

save(annot_df, ref, file = 'majorcelltypes_merged.RData')


# Immune cells
##########################################

# Clustering for Immune cells
ct_immune <- subset(annot_df, cell_type_A == 'Immune')

imm_cl <- clustering(ct_immune, n_clusters = 8, iterations = 200, seed = 55)
# 30 34 44 46 47

# Clustering Evaluation
hm_data <- ClusterEval_data(imm_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')
table(imm_cl$cluster)/nrow(imm_cl)*100

tsne_data <- ClusterEval_data(imm_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(imm_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
myeloid <- c(3,4)
th_cells <- c(5)
tc_cells <- c(6)
t_cells <- c(2,8)
nk_cells <- c(1)
other_immune <- c(7)

im_ls <- list('Myeloid'= myeloid, 'Th_cells'=th_cells, 'Tc_cells'=tc_cells, 
    'T_cells'= t_cells, 'NK_cells'=nk_cells, 'Other_immune'=other_immune)

annot_imm <- ClusterAnnotation(data = ct_immune, df_cluster = imm_cl, 
    ls_annotation = im_ls, annotation_col = 'subtype')


save(hm_data, tsne_data, umap_data, annot_imm, file = 'immunesubtypes.RData')


# Epithelial cells
##########################################

# Clustering for Epithelial cells
ct_epi <- subset(annot_df, cell_type_A == 'Epithelial')

DetermineNumberOfClusters(ct_epi, k_max=15, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

epi_cl <- clustering(ct_epi, n_clusters = 10, iterations = 200, seed = 50)
# 15 18 19 21 22 23 24 25 26 27 28 29 31 33 35 37 38 40 41 42 43 45 48 50 51

# Clustering Evaluation
hm_data <- ClusterEval_data(epi_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(epi_cl$cluster)/nrow(epi_cl)*100

tsne_data <- ClusterEval_data(epi_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(epi_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
epi_ls <- list()
for(i in 1:10){
    epi_ls[i] <- c(i)
}
names(epi_ls) <- paste0('Epithelial_',rep(1:10))

annot_epi <- ClusterAnnotation(data = ct_epi, df_cluster = epi_cl, 
    ls_annotation = epi_ls, annotation_col = 'subtype')

save(hm_data, tsne_data, umap_data, annot_epi, file = 'episubtypes.RData')


# Other cell types
##########################################

annot_imm <- change_colname(annot_imm, 'subtype', 'subtype_A')
annot_epi <- change_colname(annot_epi, 'subtype', 'subtype_A')

annot_imm['subtype_B'] <- annot_imm$subtype_A
annot_epi['subtype_B'] <- annot_epi$subtype_A

ct_mes <- subset(annot_df, cell_type_A == 'Mesenchymal')
ct_mes[,'subtype_A'] <- ct_mes$cell_type_A

ct_fib <- subset(annot_df, cell_type_A == 'Fibroblasts')
ct_fib[,'subtype_A'] <- ct_fib$cell_type_A

ct_fibm <- rbind(ct_mes, ct_fib)
ct_fibm[,'subtype_B'] <- ct_fibm$cell_type_B

ct_endo <- subset(annot_df, cell_type_A == 'Endothelial')
ct_endo[,'subtype_A'] <- ct_endo$cell_type_A
ct_endo[,'subtype_B'] <- ct_endo$cell_type_A

annot_df <- rbind(annot_imm, annot_epi, ct_fibm, ct_endo)



# Create another col in which all t cell types are labeled as t cells
annot_df['subtype_A2'] <- annot_df$subtype_A
annot_df$subtype_A2 <- as.character(annot_df$subtype_A2)
k <- which(annot_df$subtype_A %in% c('Tc_cells','Th_cells'))
annot_df[k,'subtype_A2'] <- 'T_cells'
annot_df$subtype_A2 <- factor(annot_df$subtype_A2)

annot_df['subtype_B2'] <- annot_df$subtype_B
annot_df$subtype_B2 <- as.character(annot_df$subtype_B2)
k <- which(annot_df$subtype_B %in% c('Tc_cells','Th_cells'))
annot_df[k,'subtype_B2'] <- 'T_cells'
annot_df$subtype_B2 <- factor(annot_df$subtype_B2)

# Create another col in which all epi types are labeled as epi (w/immune subtypes)
annot_df['subtype_A3'] <- annot_df$subtype_A2
annot_df$subtype_A3 <- as.character(annot_df$subtype_A3)
k <- which(annot_df$subtype_A2 %in% c(paste0('Epithelial_',rep(1:10))))
annot_df[k,'subtype_A3'] <- 'Epithelial'
annot_df$subtype_A3 <- factor(annot_df$subtype_A3)

annot_df['subtype_B3'] <- annot_df$subtype_B2
annot_df$subtype_B3 <- as.character(annot_df$subtype_B3)
k <- which(annot_df$subtype_B2 %in% c(paste0('Epithelial_',rep(1:10))))
annot_df[k,'subtype_B3'] <- 'Epithelial'
annot_df$subtype_B3 <- factor(annot_df$subtype_B3)


save(ref, annot_df, file = 'cellsubtypes.RData')




##########################################
# ***Only for first 11 tumors***
##########################################

# Immune cells
##########################################

# Clustering for Immune cells
ct_immune <- subset(annot_df, cell_type == 'immune')

imm_cl <- clustering(ct_immune, n_clusters = 8, iterations = 200, seed = 55) 
# 30 34 44 46 47

table(imm_cl$cluster)/nrow(imm_cl)*100

# Clustering Evaluation
hm_data <- ClusterEval_data(imm_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

tsne_data <- ClusterEval_data(imm_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(imm_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
myeloid <- c(4,7)
th_cells <- c(6)
tc_cells <- c(5)
t_cells <- c(1,8)
nk_cells <- c(3)
other_immune <- c(2)

im_ls <- list('Myeloid'= myeloid, 'Th_cells'=th_cells, 'Tc_cells'=tc_cells, 
    'T_cells'= t_cells, 'NK_cells'=nk_cells, 'Other_immune'=other_immune)

annot_imm <- ClusterAnnotation(data = ct_immune, df_cluster = imm_cl, 
    ls_annotation = im_ls, annotation_col = 'subtype')

prcnt_by_pt <- ClassAbundanceByPt(data=annot_imm, ptID_col = 'pt_ID', 
    class_col = 'subtype')

save(hm_data, tsne_data, umap_data, annot_imm, prcnt_by_pt, file = 'immunesubtypes.RData')


# Epithelial cells
##########################################

# Clustering for Epithelial cells
ct_epi <- subset(annot_df, cell_type == 'epithelial')

DetermineNumberOfClusters(ct_epi, k_max=15, plot=T,smooth=0.2,
                                      iter.max=200, seed = 45)

epi_cl <- clustering(ct_epi, n_clusters = 10, iterations = 200, seed = 45)
# 15 18 19 21 22 23 24 25 26 27 28 29 31 33 35 37 38 40 41 42 43 45 48 50 51
# 15 18 19 21 22 23 25 26 29 31 33 35 37 38 40 41 50

# Clustering Evaluation
hm_data <- ClusterEval_data(epi_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

table(epi_cl$cluster)/nrow(epi_cl)*100

tsne_data <- ClusterEval_data(epi_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(epi_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

# Cluster annotation
epi_ls <- list()
for(i in 1:10){
    epi_ls[i] <- c(i)
}
names(epi_ls) <- paste0('Epithelial_',rep(1:10))

annot_epi <- ClusterAnnotation(data = ct_epi, df_cluster = epi_cl, 
    ls_annotation = epi_ls, annotation_col = 'subtype')

prcnt_by_pt <- ClassAbundanceByPt(data=annot_epi, ptID_col = 'pt_ID', 
    class_col = 'subtype')

save(hm_data, tsne_data, umap_data, annot_epi, prcnt_by_pt, file = 'episubtypes.RData')


# Other cell types

ct_mes <- subset(annot_df, cell_type == 'mesenchymal')
ct_mes[,'subtype'] <- ct_mes$cell_type

ct_endo <- subset(annot_df, cell_type == 'endothelial')
ct_endo[,'subtype'] <- ct_endo$cell_type

annot_df <- rbind(annot_imm, annot_epi, ct_mes, ct_endo)



# Create another col in which all t cell types are labeled as t cells
annot_df['subtype2'] <- annot_df$subtype
k <- which(annot_df$subtype %in% c('Tc_cells','Th_cells'))
annot_df[k,'subtype2'] <- 'T_cells'

# Create another col in which all epi types are labeled as epi (w/immune subtypes)
annot_df['subtype3'] <- annot_df$subtype2
annot_df$subtype3 <- as.character(annot_df$subtype3)
k <- which(annot_df$subtype2 %in% c(paste0('Epithelial_',rep(1:10))))
annot_df[k,'subtype3'] <- 'Epithelial'



annot_df$cell_type <- factor(annot_df$cell_type)
annot_df$subtype <- factor(annot_df$subtype)
annot_df$subtype2 <- factor(annot_df$subtype2)
annot_df$subtype3 <- factor(annot_df$subtype3)

#levels(droplevels(annot_df$subtype3))

save(annot_df, file = 'cellsubtypes.RData')

load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/woCD90.RData")

save(ref, annot_df, file = 'cellsubtypes.RData')


#####################################################
# Controls clustering
#####################################################

load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/CyTOF_ADC_controls.RData")

source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ClustAnnot_functions.R")


# Clustering
dt_cl <- clustering(big_df, n_clusters = 2, iterations = 200, seed = 45) 
# 15 20 29 31 37 40
# 15, 17:31, 33:35, 37:48, 50, 51
# 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 33 34 35 37 38 39 40 41 42 43 44 45 46 47 48 50 51

# Clustering Evaluation
hm_data <- ClusterEval_data(dt_cl, eval_type = 'heatmap')
ClusterEval_plot(hm_data, data_type = 'heatmap')

tsne_data <- ClusterEval_data(dt_cl, eval_type = 'tSNE', sample_size = 25000)
ClusterEval_plot(tsne_data, data_type = 'tSNE')

umap_data <- ClusterEval_data(dt_cl, eval_type = 'UMAP', sample_size = 25000)
ClusterEval_plot(umap_data, data_type = 'UMAP')

table(dt_cl$cluster)/nrow(dt_cl)*100

# Cluster annotation
epithelial <- c(6,7,9,10)
endothelial <- c(8)
mesenchymal <- c(4)
immune <- c(1,2,3)
nothing <- c(5)
#unknown <- c(6)

ct_ls <- list('Epithelial'= epithelial, 'Endothelial'=endothelial, 
                   'Mesenchymal'=mesenchymal, 'Immune'=immune, 'Nothing'=nothing)

annot_df <- ClusterAnnotation(data = big_df, df_cluster = dt_cl, 
    ls_annotation = ct_ls, annotation_col = 'cell_type', cl_delete = T, 
    cl_delete_name = 'nothing')

prcnt_by_pt <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'cell_type')

save(hm_data, tsne_data, umap_data, annot_df, prcnt_by_pt, file = 'majorcelltypes.RData')







