

DA_analysis <- function(df, ref, class_col, CANARY_groups=c('G', 'I', 'P')){

    library(operators)
    library(dplyr)

    if('I' %!in% CANARY_groups){
        ref <- ref[-42,] # remove extra sample from 13376 patient
        df <- subset(df, CANARY %in% CANARY_groups)
        ref <- subset(ref, CANARY %in% CANARY_groups)

        prcnt_by_pt <- ClassAbundanceByPt(data=df, ptID_col = 'pt_ID', 
            class_col = class_col)
        ref <- data.frame(row.names = ref[,'pt_ID'], 
            CANARY = ref[,'CANARY'], nothing = rep(0, nrow(prcnt_by_pt)))

        # Making sure that order matches
        ref <- ref[order(rownames(ref)),]
        prcnt_by_pt <- prcnt_by_pt[order(rownames(prcnt_by_pt)),]

        # Separating matrices by group
        g_i <- which(ref$CANARY == 'G')     
        p_i <- which(ref$CANARY == 'P')

        g_t <- compositions::clr(prcnt_by_pt[g_i,])
        p_t <- compositions::clr(prcnt_by_pt[p_i,])

        RA_raw <- rbind(prcnt_by_pt[g_i,], prcnt_by_pt[p_i,])

        RA_clrt <- rbind(g_t, p_t)
        ref_t <- rbind(ref[g_i,], ref[p_i,])
        ref_t$nothing <- NULL

        pv <- c()
        for(i in 1:ncol(RA_clrt)){
            pv <- c(pv, wilcox.test(g_t[,i], p_t[,i])$p.value)
        }

        pvals <- data.frame('cell_type' = colnames(RA_clrt), 'p.value' = pv)
        pv_method <- "Wilcoxon rank sum test"

    } else {
        df <- subset(df, CANARY %in% CANARY_groups)
        ref <- subset(ref, CANARY %in% CANARY_groups)

        prcnt_by_pt <- ClassAbundanceByPt(data=df, ptID_col = 'pt_ID', 
            class_col = class_col)
        ref <- data.frame(row.names = rownames(prcnt_by_pt), 
            CANARY = ref[-42,'CANARY'], nothing = rep(0, nrow(prcnt_by_pt)))

        g_i <- which(ref$CANARY == 'G')
        i_i <- which(ref$CANARY == 'I') 
        p_i <- which(ref$CANARY == 'P')

        g_t <- compositions::clr(prcnt_by_pt[g_i,])
        i_t <- compositions::clr(prcnt_by_pt[i_i,])
        p_t <- compositions::clr(prcnt_by_pt[p_i,])

        RA_raw <- rbind(prcnt_by_pt[g_i,], prcnt_by_pt[p_i,])

        RA_clrt <- rbind(g_t, i_t, p_t)
        ref_t <- rbind(ref[g_i,], ref[i_i,], ref[p_i,])
        ref_t$nothing <- NULL

        pv <- c()
        for(i in 1:ncol(RA_clrt)){
            pv <- c(pv, kruskal.test(list(g_t[,i], i_t[,i], p_t[,i]))$p.value)
        }

        pvals <- data.frame('cell_type' = colnames(RA_clrt), 'p.value' = pv)
        pv_method <- "Kruskal-Wallis rank sum test"
    }

    DA_results <- list('RA_raw'=RA_raw, 'RA_clrt'=RA_clrt, 'ref_t'=ref_t, 
        'pvals'=pvals, 'pv_method'=pv_method)
    return(DA_results)

}


#pheatmap(prcnt_t, cluster_cols = T, cluster_rows = F, annotation_row = ref2_t, 
#    gaps_row = c(7), cellheight = 5, cellwidth = 20)