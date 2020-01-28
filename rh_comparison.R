source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/clustering_f.R")

load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/woCD90/cellsubtypes.RData")

rd <- denoisingCTF::t_asinh(annot_df[,1:5])
k <- which(rd$`103Rh`< asinh(12/5))
rh_neg_df <- annot_df[k,]
rh_pos_df <- annot_df[-k,]
markers <- colnames(rh_pos_df)[c(4, 15, 17:31, 33:35, 37:48, 50, 51)]

# Rhodium positive
dt <- denoisingCTF::t_asinh(rh_pos_df[markers])
hm_data <- aggregate(dt, list(as.factor(rh_pos_df$subtype)), median)[,-1]
hm_data <- round(as.matrix(hm_data[,1:ncol(hm_data)-1]), digits = 3)
rownames(hm_data) <- levels(rh_pos_df$subtype)
ClusterEval_plot(hm_data, data_type = 'heatmap')

# Rhodium negative
dt <- denoisingCTF::t_asinh(rh_neg_df[markers])
hm_data <- aggregate(dt, list(as.factor(rh_neg_df$subtype)), median)[,-1]
hm_data <- round(as.matrix(hm_data[,1:ncol(hm_data)-1]), digits = 3)
rownames(hm_data) <- unique(rh_neg_df$subtype)
ClusterEval_plot(hm_data, data_type = 'heatmap')

# All dt points
dt <- denoisingCTF::t_asinh(annot_df[markers])
hm_data <- aggregate(dt, list(as.factor(annot_df$subtype)), median)[,-1]
hm_data <- round(as.matrix(hm_data[,1:ncol(hm_data)-1]), digits = 3)
rownames(hm_data) <- unique(annot_df$subtype)
ClusterEval_plot(hm_data, data_type = 'heatmap')




######### Rh neg

markers <- colnames(big_df)[c(4, 15, 17:31, 33:35, 37:48, 50, 51)]

dt <- denoisingCTF::t_asinh(big_df[markers])
hm_data <- aggregate(dt, list(as.factor(dt_cl$cluster)), median)[,-1]
hm_data <- round(as.matrix(hm_data[,1:ncol(hm_data)-1]), digits = 3)
ClusterEval_plot(hm_data, data_type = 'heatmap')


rd <- denoisingCTF::t_asinh(big_df[,1:5])
k <- which(rd$`103Rh`< asinh(12/5))

rhneg <- big_df[k,]
dt_lost <- big_df[-k,]
#markers <- colnames(big_df)[c(15, 17:31, 33:35, 37:48, 50, 51)]

markers <- colnames(big_df)[c(4, 15, 17:31, 33:35, 37:48, 50, 51)]

dt <- denoisingCTF::t_asinh(rhneg[markers])
hm_data <- aggregate(dt, list(as.factor(dt_cl$cluster[k])), median)[,-1]
hm_data <- round(as.matrix(hm_data[,1:ncol(hm_data)-1]), digits = 3)
ClusterEval_plot(hm_data, data_type = 'heatmap')

dt <- denoisingCTF::t_asinh(dt_lost[markers])
hm_data <- aggregate(dt, list(as.factor(dt_cl$cluster[-k])), median)[,-1]
hm_data <- round(as.matrix(hm_data[,1:ncol(hm_data)-1]), digits = 3)
ClusterEval_plot(hm_data, data_type = 'heatmap')

