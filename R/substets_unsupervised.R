load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")

source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ClustAnnot_functions.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/DE_functions.R")

# Get percentage
ctb <- ClassAbundanceByPt(annot_df, ptID_col = 'pt_ID', class_col = 'cell_type_B')
stb2 <- ClassAbundanceByPt(annot_df, ptID_col = 'pt_ID', class_col = 'subtype_B2')
stb <- ClassAbundanceByPt(annot_df, ptID_col = 'pt_ID', class_col = 'subtype_B')

prcnt_pt <- cbind(Endothelial = ctb$Endothelial, Fib_Mesenchymal = ctb$Fib_Mesenchymal,
    T_cells = stb2$T_cells, Tc_cells = stb$Tc_cells, Th_cells = stb$Th_cells, 
    Myeloid = stb$Myeloid, NK_cells = stb$NK_cells, stb[,grep('Epithelial', colnames(stb))])


# get expression 

med_endo <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Endothelial', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_endo)[3:ncol(med_endo)] <- paste0('Endo_',colnames(med_endo)[3:ncol(med_endo)])
med_endo <- med_endo[match(rownames(ctb),med_endo$pt_ID),]
med_endo_w <- cbind(med_endo[,1:2], med_endo[,3:ncol(med_endo)]*prcnt_pt$Endothelial)


med_fibmes <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='cell_type_B',
    celltype_name='Fib_Mesenchymal', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_fibmes)[3:ncol(med_fibmes)] <- paste0('FMes_',colnames(med_fibmes)[3:ncol(med_fibmes)])
med_fibmes <- med_fibmes[match(rownames(ctb),med_fibmes$pt_ID),]
med_fibmes_w <- cbind(med_fibmes[,1:2], med_fibmes[,3:ncol(med_fibmes)]*prcnt_pt$Fib_Mesenchymal)


med_tcells <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B2',
    celltype_name='T_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_tcells)[3:ncol(med_tcells)] <- paste0('Tcells_',colnames(med_tcells)[3:ncol(med_tcells)])
med_tcells <- med_tcells[match(rownames(ctb),med_tcells$pt_ID),]
med_tcells_w <- cbind(med_tcells[,1:2], med_tcells[,3:ncol(med_tcells)]*prcnt_pt$T_cells)


med_CD8T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Tc_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_CD8T)[3:ncol(med_CD8T)] <- paste0('TCD8_',colnames(med_CD8T)[3:ncol(med_CD8T)])
med_CD8T <- med_CD8T[match(rownames(ctb),med_CD8T$pt_ID),]
med_CD8T$pt_ID <- med_fibmes$pt_ID # assuming med_fibmes has no NAs
med_CD8T$CANARY <- med_fibmes$CANARY # assuming med_fibmes has no NAs
med_CD8T_w <- cbind(med_CD8T[,1:2], med_CD8T[,3:ncol(med_CD8T)]*prcnt_pt$Tc_cells)



med_CD4T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Th_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_CD4T)[3:ncol(med_CD4T)] <- paste0('TCD4_',colnames(med_CD4T)[3:ncol(med_CD4T)])
med_CD4T <- med_CD4T[match(rownames(ctb),med_CD4T$pt_ID),]
med_CD4T_w <- cbind(med_CD4T[,1:2], med_CD4T[,3:ncol(med_CD4T)]*prcnt_pt$Th_cells)


med_mye <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='Myeloid', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_mye)[3:ncol(med_mye)] <- paste0('Mye_',colnames(med_mye)[3:ncol(med_mye)])
med_mye <- med_mye[match(rownames(ctb),med_mye$pt_ID),]
med_mye_w <- cbind(med_mye[,1:2], med_mye[,3:ncol(med_mye)]*prcnt_pt$Myeloid)


med_NK <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name='NK_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_NK)[3:ncol(med_NK)] <- paste0('NK_',colnames(med_NK)[3:ncol(med_NK)])
med_NK <- med_NK[match(rownames(ctb),med_NK$pt_ID),]
med_NK_w <- cbind(med_NK[,1:2], med_NK[,3:ncol(med_NK)]*prcnt_pt$NK_cells)


med_epi <- list()
med_epi_w <- list()
e <- length(grep('Epithelial', unique(annot_df$subtype_B)))
for (i in 1:e){
    med_epi[[i]] <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B',
    celltype_name=paste0('Epithelial_', i), ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
    med_epi[[i]] <- med_epi[[i]][match(rownames(ctb),med_epi[[i]]$pt_ID),]
    med_epi[[i]]$pt_ID <- med_fibmes$pt_ID # assuming med_fibmes has no NAs
    med_epi[[i]]$CANARY <- med_fibmes$CANARY # assuming med_fibmes has no NAs
    med_epi_w[[i]] <- cbind(med_epi[[i]][,1:2], med_epi[[i]][,3:ncol(med_epi[[i]])]*prcnt_pt[,paste0('Epithelial_', i)])
}
names(med_epi) <- paste0('Epi_', rep(1:e))
x <- do.call("cbind", med_epi)
k <- c(grep('CANARY', colnames(x)), grep('pt_ID', colnames(x)))
med_epi <- x[,-sort(k)[3:length(k)]]

names(med_epi_w) <- paste0('Epi_', rep(1:e))
x <- do.call("cbind", med_epi_w)
k <- c(grep('CANARY', colnames(x)), grep('pt_ID', colnames(x)))
med_epi_w <- x[,-sort(k)[3:length(k)]]



# Merging all
raw_exp <- do.call('cbind', 
    list(med_endo, med_fibmes, med_tcells, med_CD8T, med_CD4T, med_mye, med_NK, med_epi))

k <- c(grep('CANARY', colnames(raw_exp)), grep('pt_ID', colnames(raw_exp)))
raw_exp <- raw_exp[,-sort(k)[3:length(k)]]
raw_exp[is.na(raw_exp)] <- 0


weighted_exp <- do.call('cbind', 
    list(med_endo_w, med_fibmes_w, med_tcells_w, med_CD8T_w, med_CD4T_w, med_mye_w, med_NK_w, med_epi_w))

k <- c(grep('CANARY', colnames(weighted_exp)), grep('pt_ID', colnames(weighted_exp)))
weighted_exp <- weighted_exp[,-sort(k)[3:length(k)]]
weighted_exp[is.na(weighted_exp)] <- 0



#############################################################################################################
#############################################################################################################


DetermineNumberOfClusters(raw_exp[,3:ncol(raw_exp)], k_max = 20, plot = T,
    ask_ft = F, arcsn_tr = F) #4

set.seed(55)
cl_4 <- kmeans(raw_exp[,3:ncol(raw_exp)], centers=4, iter.max = 100)$cluster
cl_8 <- kmeans(raw_exp[,3:ncol(raw_exp)], centers=8, iter.max = 100)$cluster

tsne_smp <- Rtsne::Rtsne(raw_exp[,3:ncol(raw_exp)], 
    check_duplicates = FALSE, perplexity =14)

tsne_smp <- as.data.frame(cbind(cluster4 = cl_4, cluster8 = cl_8,
                                tSNE1 = tsne_smp$Y[,1],
                                tSNE2 = tsne_smp$Y[,2]))

plot(tsne_smp$tSNE1, tsne_smp$tSNE2, col=tsne_smp$cluster4, pch=19)


DetermineNumberOfClusters(weighted_exp[,3:ncol(weighted_exp)], k_max = 20, plot = T,
    ask_ft = F, arcsn_tr = F) 

set.seed(55)
cl_4 <- kmeans(weighted_exp[,3:ncol(weighted_exp)], centers=4, iter.max = 100)$cluster
cl_8 <- kmeans(weighted_exp[,3:ncol(weighted_exp)], centers=8, iter.max = 100)$cluster

tsne_smp <- Rtsne::Rtsne(weighted_exp[,3:ncol(weighted_exp)], 
    check_duplicates = FALSE, perplexity =14)

tsne_smp <- as.data.frame(cbind(cluster4 = cl_4, cluster8 = cl_8,
                                tSNE1 = tsne_smp$Y[,1],
                                tSNE2 = tsne_smp$Y[,2]))

plot(tsne_smp$tSNE1, tsne_smp$tSNE2, col=tsne_smp$cluster4, pch=19)




stb_pr <- stb/100

DetermineNumberOfClusters(stb_pr, k_max = 20, plot = T,
    ask_ft = F, arcsn_tr = F) #8

set.seed(55)
cl_4 <- kmeans(stb_pr, centers=4, iter.max = 100)$cluster
cl_8 <- kmeans(stb_pr, centers=8, iter.max = 100)$cluster

tsne_smp <- Rtsne::Rtsne(stb_pr, 
    check_duplicates = FALSE, perplexity =20)

tsne_smp <- as.data.frame(cbind(cluster4 = cl_4, cluster8 = cl_8,
                                tSNE1 = tsne_smp$Y[,1],
                                tSNE2 = tsne_smp$Y[,2]))

plot(tsne_smp$tSNE1, tsne_smp$tSNE2, col=tsne_smp$cluster4, pch=19)



DetermineNumberOfClusters(test[,3:ncol(test)], k_max = 10, plot = T,
    ask_ft = F, arcsn_tr = F) #8

set.seed(55)
cl_4 <- kmeans(test[,3:ncol(test)], centers=3, iter.max = 100)$cluster
cl_8 <- kmeans(test[,3:ncol(test)], centers=8, iter.max = 100)$cluster

tsne_smp <- Rtsne::Rtsne(stb_pr, 
    check_duplicates = FALSE, perplexity =20)

tsne_smp <- as.data.frame(cbind(cluster4 = cl_4, cluster8 = cl_8,
                                tSNE1 = tsne_smp$Y[,1],
                                tSNE2 = tsne_smp$Y[,2]))

plot(tsne_smp$tSNE1, tsne_smp$tSNE2, col=tsne_smp$cluster4, pch=19)




library(plotly)
p <- plot_ly(x = colnames(weighted_exp)[3:ncol(weighted_exp)], y = as.factor(weighted_exp$ptID_col),
    z = as.matrix(weighted_exp[,3:ncol(weighted_exp)]), type = "heatmap")
p

test <- cbind(weighted_exp[,1:2], scale(weighted_exp[3:ncol(weighted_exp)]))

heatmap(as.matrix(test[,3:ncol(test)]))

heatmap(as.matrix(weighted_exp[,3:ncol(weighted_exp)]), scale='none')

heatmap(as.matrix(raw_exp[,3:ncol(raw_exp)]),  scale='none')

heatmap(as.matrix(prcnt_pt), scale='col')



x <- compositions::clr(stb)

heatmap(as.matrix(x), scale='none')

heatmap(as.matrix(stb), scale='none')



    

