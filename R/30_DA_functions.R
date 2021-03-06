    
DA_analysis <- function(df, ref, class_col, group_levels, group_col, 
    ptID_col = 'pt_ID', transformation = FALSE){


    if(length(group_levels) != 2){
        stop('Differential abundance analysis requieres a vector of length 2 for group_levels')
    }

    library(operators)
    library(dplyr)

    ref <- ref[-42,] # remove extra sample from 13376 patient
    df <- df[df[,group_col] %in% group_levels,]
    ref <- ref[ref[,group_col] %in% group_levels,]

    prcnt_by_pt <- ClassAbundanceByPt(data=df, ptID_col = ptID_col, 
        class_col = class_col)
    ref <- data.frame(row.names = ref[,ptID_col], 
        CANARY = ref[,group_col], nothing = rep(0, nrow(prcnt_by_pt)))

    # Making sure that order matches
    ref <- ref[order(rownames(ref)),]
    prcnt_by_pt <- prcnt_by_pt[order(rownames(prcnt_by_pt)),]

    # Separating matrices by group
    g_i <- which(ref[,group_col] == group_levels[1])
    p_i <- which(ref[,group_col] == group_levels[2])

    if (transformation){
        g_t <- compositions::clr(prcnt_by_pt[g_i,])
        p_t <- compositions::clr(prcnt_by_pt[p_i,])
    } else {
        g_t <- prcnt_by_pt[g_i,]
        p_t <- prcnt_by_pt[p_i,]
    }

    RA_raw <- rbind(prcnt_by_pt[g_i,], prcnt_by_pt[p_i,])

    RA_clrt <- rbind(g_t, p_t)
    ref_t <- rbind(ref[g_i,], ref[p_i,])
    ref_t$nothing <- NULL

    pv <- c()
    for(i in 1:ncol(RA_clrt)){
        pv <- c(pv, wilcox.test(g_t[,i], p_t[,i])$p.value)
    }

    pvals <- data.frame('cell_type' = colnames(RA_clrt), 'p.value' = pv)
    corrected_pvals <- p.adjust(pvals$p.value, method = 'BH')
    pvals['FDR_corrected'] <- corrected_pvals

    pv_method <- "Wilcoxon rank sum test"
    correction_method <- "Benjamini & Hochberg / FDR"

    DA_results <- list('RA_raw'=RA_raw, 'RA_clrt'=RA_clrt, 'ref_t'=ref_t, 
        'pvals'=pvals, 'pv_method'=pv_method, 'correction_method' = correction_method)

    return(DA_results)
}



#pheatmap(prcnt_t, cluster_cols = T, cluster_rows = F, annotation_row = ref2_t, 
#    gaps_row = c(7), cellheight = 5, cellwidth = 20)