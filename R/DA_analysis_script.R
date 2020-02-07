# Load data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")

# Source functions
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/DA_analysis.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/ClustAnnot_functions.R")


library(gplots)
library(RColorBrewer)
library(pheatmap)


######################################################
# Good and Poor 
######################################################

## All celltypes

ctB <- DA_analysis(annot_df, ref, class_col='cell_type_B', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')

stB <- DA_analysis(annot_df, ref, class_col='subtype_B', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')

stB2 <- DA_analysis(annot_df, ref, class_col='subtype_B2', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')

stB3 <- DA_analysis(annot_df, ref, class_col='subtype_B3', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')


## Immune only

k <- which(annot_df$cell_type_B == 'Immune')
sbst <- annot_df[k,]
sbst$subtype_B <- factor(sbst$subtype_B)
sbst$subtype_B2 <- factor(sbst$subtype_B2)

stB_imm <- DA_analysis(sbst, ref, class_col='subtype_B', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')

stB2_imm <- DA_analysis(sbst, ref, class_col='subtype_B2', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')

## Epithelial only
k <- which(annot_df$cell_type_B == 'Epithelial')
sbst <- annot_df[k,]
sbst$subtype_B <- factor(sbst$subtype_B)

stB_epi <- DA_analysis(sbst, ref, class_col='subtype_B', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')


## Stromal only
sbst <- annot_df[-k,]
sbst$subtype_B <- factor(sbst$subtype_B)
sbst$subtype_B2 <- factor(sbst$subtype_B2)

stB_str <- DA_analysis(sbst, ref, class_col='subtype_B', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')

stB2_str <- DA_analysis(sbst, ref, class_col='subtype_B2', group_levels=c('G','P'),
     group_col = 'CANARY', ptID_col = 'pt_ID')



# ctB
pheatmap(ctB_test$RA_clrt, cluster_cols = T, cluster_rows = F, annotation_row = ctB_test$ref_t, 
    gaps_row = c(16), cellheight = 5, cellwidth = 20)

# stB_epi
pheatmap(stB_epi$RA_clrt, cluster_cols = T, cluster_rows = F, annotation_row = stB_epi$ref_t, 
    gaps_row = c(7), cellheight = 5, cellwidth = 20)






