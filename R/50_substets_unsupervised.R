load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")

source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/20_ClustAnnot_functions.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/40_DE_functions.R")


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
sbst_exp <- do.call('cbind', 
    list(med_endo, med_fibmes, med_tcells, med_CD8T, med_CD4T, med_mye, med_NK, med_epi))

k <- c(grep('CANARY', colnames(sbst_exp)), grep('pt_ID', colnames(sbst_exp)))
sbst_exp <- sbst_exp[,-sort(k)[3:length(k)]]
sbst_exp[is.na(sbst_exp)] <- 0


save(sbst_exp, file= 'subset_exprmat.RData')


#############################################################################################################
# Integration of Clinical data

load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/subset_exprmat.RData")

library(readxl)
CDE_TMA36 <- read_excel("~/Documents/Massion_lab/CDE/CDE_TMA36_2020FEB25_SA.xlsx", 
    sheet = "ADC_mafe_processed", col_types = c("date", 
        "text", "text", "date", "text", "text", 
        "numeric", "numeric", "text", "text", 
        "text", "numeric", "text", "numeric", 
        "numeric", "numeric", "text", "text", 
        "text", "text", "date", "numeric", 
        "text", "text", "text", "text", "date", 
        "text", "date", "text", "text", "date", 
        "text", "text", "text", "text", "text", 
        "text", "numeric", "text", "text", 
        "text", "date", "text", "date", "text", 
        "date", "text", "date", "date", "text", 
        "text", "text", "date", "text", "text", 
        "text", "text", "date", "text", "text", 
        "text", "date", "text", "text", "date", 
        "text"))

x <- match(sbst_exp$pt_ID, CDE_TMA36$Patient_ID)
pData_cytof <- CDE_TMA36[x,]
colnames(pData_cytof)[2] <- 'pt_ID'
pData_cytof[9,3] <- 'N'
pData_cytof$Family_History_Cancer_Type[9] <- 'Unknown'

# Option 1: RData object with 2 data_frames
save(sbst_exp, pData_cytof, file= 'subset_exprmat_clinical_annotations.RData')


#2 Option 2: ExpressionSet Object containing assay data and clinical annotations
rownames(pData_cytof) <- paste0(rep('smp'), '_', 1:nrow(pData_cytof))
sbst_exp$CANARY <- NULL
sbst_exp$pt_ID <- NULL
rownames(sbst_exp) <- paste0(rep('smp'), '_', 1:nrow(sbst_exp))
aData_cytof <- t(sbst_exp)

library(Biobase)

eSet_cytof <- ExpressionSet(assayData = aData_cytof, phenoData = AnnotatedDataFrame(pData_cytof))

save(eSet_cytof, file= 'eSet_cytof.RData')

#############################################################################################################

# Unsupervised analysis of patients proteocmic profile

load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/subset_exprmat_clinical_annotations.RData")

excl_i <- grep("*EpCAM*|*Cytokeratin*|*CK7*|*CD31*|*CD45*|*CD56*|*CD8*|*CD3*|*CD11b*|*CD90*|*CD4*|*Vimentin|*TP63*", colnames(sbst_exp))

sbst_exp <- sbst_exp[,-excl_i]

set.seed(53)

DetermineNumberOfClusters(sbst_exp[,3:ncol(sbst_exp)], k_max = 20, plot = T,
    ask_ft = F, arcsn_tr = F) #4


library(ComplexHeatmap)
library(RColorBrewer)

data <- as.matrix(sbst_exp[,3:ncol(sbst_exp)])
scale_max = max(data)
heat_palette_med <- c("Darkblue", "white", "red")
pairs.breaks_med <- c(0,1.5,scale_max)

library(circlize)

ha = rowAnnotation(
    foo = runif(71), 
    bar = sample(letters[1:3], 71, replace = TRUE),
    col = list(foo = col_fun,
               bar = c("a" = "red", "b" = "green", "c" = "blue")
    ),
    #gp = gpar(col = "black"),
    simple_anno_size = unit(0.5, "cm")
)


    canary = as.factor(pData_cytof$CANARY), 
    gender = as.factor(pData_cytof$Gender),
    smoking = as.factor(pData_cytof$Smoking_Status),
    stage = as.factor(pData_cytof$`8th_edition_path_stage`),
    age = pData_cytof$Age_at_collection,
    life_status = as.factor(pData_cytof$Living_Status),
    BMI = pData_cytof$BMI,
    family_cancer = as.factor(pData_cytof$Family_History_Cancer_Type),
    prior_cancer = as.factor(pData_cytof$Prior_Cancer)
    prior_cancer = as.factor(pData_cytof$Prior_Cancer_Type)
    fev1 = as.numeric(pData_cytof$`FEV1 (% Pred)`)



ha = rowAnnotation(
    #fev1 = as.numeric(pData_cytof$`FEV1 (% Pred)`),
    smoking = as.factor(pData_cytof$Smoking_Status),
    stage = as.factor(pData_cytof$`8th_edition_path_stage`),
    life_status = as.factor(pData_cytof$Living_Status),
    simple_anno_size = unit(0.5, "cm")
)

Heatmap(data, name = "mat", row_km = 4, column_km = 5,
    col = circlize::colorRamp2(pairs.breaks_med, heat_palette_med),
  heatmap_legend_param = list(color_bar = "continuous"), right_annotation = ha)








cl_4 <- kmeans(sbst_exp[,3:ncol(sbst_exp)], centers=4, iter.max = 100)$cluster

tsne_smp <- Rtsne::Rtsne(sbst_exp[,3:ncol(sbst_exp)], 
    check_duplicates = FALSE, perplexity =20, pca_scale = F)

pca_smp <- prcomp(sbst_exp[,3:ncol(sbst_exp)], center = TRUE,scale. = TRUE)

plot(tsne_smp$Y[,1], tsne_smp$Y[,2], col=cl_4, pch=19)
plot(pca_smp$x[,1], pca_smp$x[,2], col=cl_4, pch=19)