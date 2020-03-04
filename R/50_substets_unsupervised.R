load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")

source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/20_ClustAnnot_functions.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/40_DE_functions.R")

# Get percentage
# ctb <- ClassAbundanceByPt(annot_df, ptID_col = 'pt_ID', class_col = 'cell_type_B')
# stb2 <- ClassAbundanceByPt(annot_df, ptID_col = 'pt_ID', class_col = 'subtype_B2')
# stb <- ClassAbundanceByPt(annot_df, ptID_col = 'pt_ID', class_col = 'subtype_B')

# prcnt_pt <- cbind(Endothelial = ctb$Endothelial, Fib_Mesenchymal = ctb$Fib_Mesenchymal,
#     T_cells = stb2$T_cells, Tc_cells = stb$Tc_cells, Th_cells = stb$Th_cells, 
#     Myeloid = stb$Myeloid, NK_cells = stb$NK_cells, stb[,grep('Epithelial', colnames(stb))])


# get expression

med_all <- median_by_pt(annot_df, ref, subset_celltype=F,
    ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)

match_all <- function(med_subtype, med_all){
    med_subtype <- med_subtype[match(med_all$pt_ID,med_subtype$pt_ID),]
    med_subtype$pt_ID <- med_all$pt_ID # assuming med_fibmes has no NAs
    med_subtype$CANARY <- med_all$CANARY # assuming med_fibmes has no NAs

    med_subtype
}

med_endo <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Endothelial', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_endo)[3:ncol(med_endo)] <- paste0('Endo_',colnames(med_endo)[3:ncol(med_endo)])
med_endo <- match_all(med_endo, med_all)
# med_endo_w <- cbind(med_endo[,1:2], med_endo[,3:ncol(med_endo)]*prcnt_pt$Endothelial)


med_fibmes <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Fib_Mesenchymal', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_fibmes)[3:ncol(med_fibmes)] <- paste0('FMes_',colnames(med_fibmes)[3:ncol(med_fibmes)])
med_fibmes <- match_all(med_fibmes, med_all)


med_tcells <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B2',
    celltype_name='T_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_tcells)[3:ncol(med_tcells)] <- paste0('Tcells_',colnames(med_tcells)[3:ncol(med_tcells)])
med_tcells <- match_all(med_tcells, med_all)


med_CD8T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Tc_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_CD8T)[3:ncol(med_CD8T)] <- paste0('Tc_',colnames(med_CD8T)[3:ncol(med_CD8T)])
med_CD8T <- match_all(med_CD8T, med_all)



med_CD4T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Th_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_CD4T)[3:ncol(med_CD4T)] <- paste0('Th_',colnames(med_CD4T)[3:ncol(med_CD4T)])
med_CD4T <- match_all(med_CD4T, med_all)


med_mye <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Myeloid', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_mye)[3:ncol(med_mye)] <- paste0('Mye_',colnames(med_mye)[3:ncol(med_mye)])
med_mye <- match_all(med_mye, med_all)


med_NK <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='NK_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_NK)[3:ncol(med_NK)] <- paste0('NK_',colnames(med_NK)[3:ncol(med_NK)])
med_NK <- match_all(med_NK, med_all)


# Epithelial By subset
med_epi <- list()
# med_epi_w <- list()
e <- length(grep('Epithelial', unique(annot_df$subtype_B)))
for (i in 1:e){
    med_epi[[i]] <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name=paste0('Epithelial_', i), ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
    med_epi[[i]] <- match_all(med_epi[[i]], med_all)
}
names(med_epi) <- paste0('Epi_', rep(1:e))
x <- do.call("cbind", med_epi)
k <- c(grep('CANARY', colnames(x)), grep('pt_ID', colnames(x)))
med_epi <- x[,-sort(k)[3:length(k)]]


# Epithelial as one
med_epi <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Epithelial', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_epi)[3:ncol(med_epi)] <- paste0('Epi_',colnames(med_epi)[3:ncol(med_epi)])
med_epi <- match_all(med_epi, med_all)


# Merging all
raw_exp <- do.call('cbind', 
    list(med_endo, med_fibmes, med_tcells, med_CD8T, med_CD4T, med_mye, med_NK, med_epi))

k <- c(grep('CANARY', colnames(raw_exp)), grep('pt_ID', colnames(raw_exp)))
raw_exp <- raw_exp[,-sort(k)[3:length(k)]]
raw_exp[is.na(raw_exp)] <- 0





#############################################################################################################
#############################################################################################################
excl_i <- grep("*EpCAM*|*Cytokeratin*|*CK7*|*CD31*|*CD45*|*CD56*|*CD8*|*CD3*|*CD11b*|*CD90*|*CD4*|*Vimentin|*TP63*", colnames(raw_exp))

raw_exp <- raw_exp[,-excl_i]

set.seed(53)

DetermineNumberOfClusters(raw_exp[,3:ncol(raw_exp)], k_max = 20, plot = T,
    ask_ft = F, arcsn_tr = F) #4

cl_4 <- kmeans(raw_exp[,3:ncol(raw_exp)], centers=4, iter.max = 100)$cluster

tsne_smp <- Rtsne::Rtsne(raw_exp[,3:ncol(raw_exp)], 
    check_duplicates = FALSE, perplexity =20, pca_scale = F)

pca_smp <- prcomp(raw_exp[,3:ncol(raw_exp)], center = TRUE,scale. = TRUE)

plot(tsne_smp$Y[,1], tsne_smp$Y[,2], col=cl_4, pch=19)
plot(pca_smp$x[,1], pca_smp$x[,2], col=cl_4, pch=19)


library(ComplexHeatmap)
data <- as.matrix(raw_exp[,3:ncol(raw_exp)])
scale_max = max(data)
heat_palette_med <- c("Darkblue", "white", "red")
pairs.breaks_med <- c(0,1.5,scale_max)
Heatmap(data, name = "mat", row_km = 4, column_km = 5,
    col = circlize::colorRamp2(pairs.breaks_med, heat_palette_med),
  heatmap_legend_param = list(color_bar = "continuous"))






    

