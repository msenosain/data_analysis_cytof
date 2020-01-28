
# Selecting columns of interest
wo_cd90 <- colnames(big_df)[c(15, 17:31, 33:35, 37:48, 50, 51)]
w_cd90 <- colnames(big_df)[c(15, 17:31, 33:35, 37:51)]
few <- c("141Pr_EpCAM", "145Nd_CD31", "154Sm_CD45",  "156Gd_Vimentin", "161Dy_Cytokeratin", "164Dy_CK7")

c_cols <- few
scaled_smp <- denoisingCTF::t_asinh(big_df[,c_cols])

# k-means clustering
seed <- 101
n_clusters <- 10
set.seed(seed)
cl_data <- stats::kmeans(scaled_smp, centers = n_clusters, iter.max = 200)$cluster

# Computing median by cluster by column
smp_cl <- cbind(scaled_smp, 'cluster'= cl_data)
smp_median <- aggregate(smp_cl, list(smp_cl$cluster), median)[,-1]

# Plotting a heatmap
x <- round(as.matrix(smp_median[,1:ncol(smp_median)-1]), digits = 3)
library(gplots)
scale_max = max(x)
heat_palette_med <- colorRampPalette(c("black", "yellow","#FAF7C9"))
pairs.breaks_med <- c(seq(0, scale_max/6.6, by = 0.1),
                      seq(scale_max/6.6, scale_max/3.3, by = 0.1),
                      seq(scale_max/3.3, scale_max, by = 0.1))

heatmap.2(x,
          main = "Median protein expression",
          dendrogram = "both",
          Rowv = TRUE,
          Colv = TRUE,
          breaks = pairs.breaks_med,
          revC = FALSE,
          symkey = FALSE,
          symbreaks = FALSE,
          scale = "none",
          cexRow = 0.8,
          cexCol = 0.8,
          key = TRUE,
          col = heat_palette_med,
          trace = "none",
          density.info = 'none',
          sepcolor="#424242",
          margins = c(6,10),
          colsep=1:ncol(x),
          rowsep=1:nrow(x),
          sepwidth=c(0.005,0.005),
          keysize = 1,
          key.title = 'Intensity',
          key.xlab= "Arcsinh Transform",
          extrafun = box(lty = "solid"),
          cellnote = x,
          notecol = 'red',
          srtCol=45
          )
save(x, file= 'k10heatmapdata.RData')

# Generating a tSNE plot

# tsne
smpl <- dplyr::sample_n(smp_cl, size=100000, replace=F) # just a sample
tsne_smp <- Rtsne::Rtsne(smpl,  check_duplicates = FALSE)

# Editing df
scaled_smp_cl <- as.data.frame(cbind(smpl,
                                     tsne1 = tsne_smp$Y[,1],
                                     tsne2 = tsne_smp$Y[,2]))
# Save scaled_smp_cl
save(scaled_smp_cl, file = 'kmns_tsne_smpl.RData')

# Plotting
library(ggplot2)
library(gridExtra)
library(reshape2)
test <- melt(scaled_smp_cl, id= c('cluster', 'tsne1', 'tsne2'))
var_list <- unique(test$variable)
pl <- list()

# Defining cluster centers
edata <- cbind(scaled_smp_cl$tsne1, scaled_smp_cl$tsne2, scaled_smp_cl$cluster)
colnames(edata) <- c('x', "y", "z")
center <- aggregate(cbind(x,y) ~ z, data = edata, median)

p_cl <- ggplot(scaled_smp_cl, aes(x=tsne1, y=tsne2, colour = factor(cluster)))+
  geom_point(alpha=0.3) + theme_bw() + ggtitle('cluster')+
  annotate("text", label = center[,1], x=center[,2], y = center[,3],
           size = 6, colour = "black", fontface = 'bold')
pl[[1]] <- p_cl
for(i in seq_along(var_list)) {
  p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=tsne1, y=tsne2, colour = value)) +
    geom_point(alpha=0.3) + theme_bw() + ggtitle(var_list[i]) +
    scale_colour_gradient(low = "gray75", high = "blue")
  pl[[i+1]] <- p
}

grid.arrange(grobs=pl)

# Cluster annotation
epithelial <- c(4,5,6,8)
endothelial <- c(9)
mesenchymal <- c(7)
immune <- c(2,3,10)
nothing <- c(1)
#unknown <- c(6)

# Cluster merge
ct_ls <- list('epithelial'= epithelial, 'endothelial'=endothelial, 
                   'mesenchymal'=mesenchymal, 'immune'=immune, 'nothing'=nothing)

ct_k <- data.frame(cbind(cluster=cl_data, cell_type=rep(NA, length(cl_data))))
for (i in 1:length(ct_ls)){
  k <- which(smp_cl$cluster %in% ct_ls[[i]])
  ct_k[k,'cell_type'] <- names(ct_ls)[i]
}


# Translate to bigdf
bdf <- cbind(denoisingCTF::t_asinh(big_df[,w_cd90]), 
    pt_ID=big_df$pt_ID, CD90=big_df$CD90, ct_k)

# Check general abundances (plot)
table(bdf$cell_type)/nrow(bdf)*100

# Heatmap of all markers (to see how the weird clusters behave, 1 and 6)
x_cl <- aggregate(bdf[,1:33], list(bdf$cluster), median)
x_cl <- round(as.matrix(x_cl[,2:ncol(x_cl)]), digits = 3)
library(gplots)
scale_max = max(x_cl)
heat_palette_med <- colorRampPalette(c("black", "yellow","#FAF7C9"))
pairs.breaks_med <- c(seq(0, scale_max/6.6, by = 0.1),
                      seq(scale_max/6.6, scale_max/3.3, by = 0.1),
                      seq(scale_max/3.3, scale_max, by = 0.1))

heatmap.2(x_cl,
          main = "Median protein expression",
          dendrogram = "both",
          Rowv = TRUE,
          Colv = TRUE,
          breaks = pairs.breaks_med,
          revC = FALSE,
          symkey = FALSE,
          symbreaks = FALSE,
          scale = "none",
          cexRow = 0.8,
          cexCol = 0.8,
          key = TRUE,
          col = heat_palette_med,
          trace = "none",
          density.info = 'none',
          sepcolor="#424242",
          margins = c(6,10),
          colsep=1:ncol(x_cl),
          rowsep=1:nrow(x_cl),
          sepwidth=c(0.005,0.005),
          keysize = 1,
          key.title = 'Intensity',
          key.xlab= "Arcsinh Transform",
          extrafun = box(lty = "solid"),
          cellnote = x_cl,
          notecol = 'red',
          srtCol=45
          )

save(x_cl, file= 'k10heatmapdataALLMARKERS.RData')
save(bdf, file= 'woCD90dataclusters.RData')

# Remove cluster 1 ("nothing")
k <- which(bdf$cell_type == 'nothing')
bdf <- bdf[-k,]

# Compute cell type abundances by pt
ptids <- unique(bdf$pt_ID)
prcnt <- table(bdf$cell_type)

for(i in ptids){
    k <- subset(bdf, pt_ID==i)
    k <-table(k$cell_type)/nrow(k)*100
    prcnt <- rbind(prcnt,k)
}

prcnt <- data.frame(prcnt)[-1,]
rownames(prcnt) <- ptids
save(prcnt, file= 'celltypespercents.RData')

##########################################
# ***Only for first 11 tumors***
##########################################

CANARY <- c('I', 'G', 'P', 'G', 'P', 'G', 'P', 'P', 'G', 'P', 'P')
ref <- cbind(ref, CANARY)

# Create separate dfs for each cell type
ct_immune <- subset(bdf, cell_type == 'immune')
ct_epi <- subset(bdf, cell_type == 'epithelial')
ct_mes <- subset(bdf, cell_type == 'mesenchymal')
ct_endo <- subset(bdf, cell_type == 'endothelial')

# Clustering for Immune cells

# markers
c_cols <- c("155Gd_CD56", "159Tb_CD4", "168Er_CD8", "170Yb_CD3", "171Yb_CD11b")
# "144Nd_HLA-ABC", "174Yb_HLA-DR", "156Gd_Vimentin"
# "166Er_CD44", "169Tm_CD24", "154Sm_CD45", "156Gd_Vimentin"
ct_immune_m <- ct_immune[,c_cols]

# k-means clustering
seed <- 101
n_clusters <- 8
set.seed(seed)
cl_immune <- stats::kmeans(ct_immune_m, centers = n_clusters, iter.max = 500)$cluster

# Computing median by cluster by column
ct_immune_m <- cbind(ct_immune_m, 'cluster'= cl_immune)
ct_immune_m_agg <- aggregate(ct_immune_m, list(ct_immune_m$cluster), median)[,-1]

# Plotting a heatmap
x <- round(as.matrix(ct_immune_m_agg[,1:ncol(ct_immune_m_agg)-1]), digits = 3)
library(gplots)
scale_max = max(x)
heat_palette_med <- colorRampPalette(c("black", "yellow","#FAF7C9"))
pairs.breaks_med <- c(seq(0, scale_max/6.6, by = 0.1),
                      seq(scale_max/6.6, scale_max/3.3, by = 0.1),
                      seq(scale_max/3.3, scale_max, by = 0.1))

heatmap.2(x,
          main = "Median protein expression",
          dendrogram = "both",
          Rowv = TRUE,
          Colv = TRUE,
          breaks = pairs.breaks_med,
          revC = FALSE,
          symkey = FALSE,
          symbreaks = FALSE,
          scale = "none",
          cexRow = 0.8,
          cexCol = 0.8,
          key = TRUE,
          col = heat_palette_med,
          trace = "none",
          density.info = 'none',
          sepcolor="#424242",
          margins = c(6,10),
          colsep=1:ncol(x),
          rowsep=1:nrow(x),
          sepwidth=c(0.005,0.005),
          keysize = 1,
          key.title = 'Intensity',
          key.xlab= "Arcsinh Transform",
          extrafun = box(lty = "solid"),
          cellnote = x,
          notecol = 'red',
          srtCol=45
          )
save(x, file= 'k10heatmapdata_immune.RData')

# Generating a tSNE plot

# tsne
smpl_immune <- dplyr::sample_n(ct_immune_m, size=50000, replace=F) # just a sample
tsne_smp <- Rtsne::Rtsne(smpl_immune[,-11],  check_duplicates = FALSE)

# Editing df
tsne_immune <- as.data.frame(cbind(smpl_immune,
                                     tsne1 = tsne_smp$Y[,1],
                                     tsne2 = tsne_smp$Y[,2]))
# Save tsne_immune
save(tsne_immune, file = 'kmns_tsne_immune.RData')

# Plotting
library(ggplot2)
library(gridExtra)
library(reshape2)
test <- melt(tsne_immune, id= c('cluster', 'tsne1', 'tsne2'))
var_list <- unique(test$variable)
pl <- list()

# Defining cluster centers
edata <- cbind(tsne_immune$tsne1, tsne_immune$tsne2, tsne_immune$cluster)
colnames(edata) <- c('x', "y", "z")
center <- aggregate(cbind(x,y) ~ z, data = edata, median)

p_cl <- ggplot(tsne_immune, aes(x=tsne1, y=tsne2, colour = factor(cluster)))+
  geom_point(alpha=0.3) + theme_bw() + ggtitle('cluster')+
  annotate("text", label = center[,1], x=center[,2], y = center[,3],
           size = 6, colour = "black", fontface = 'bold')
pl[[1]] <- p_cl
for(i in seq_along(var_list)) {
  p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=tsne1, y=tsne2, colour = value)) +
    geom_point(alpha=0.3) + theme_bw() + ggtitle(var_list[i]) +
    scale_colour_gradient(low = "gray75", high = "blue")
  pl[[i+1]] <- p
}

grid.arrange(grobs=pl)

# Cluster annotation
myeloid <- c(1,2)
th_cells <- c(3)
tc_cells <- c(4)
t_cells <- c(6,8)
nk_cells <- c(5)
other_immune <- c(7)

# Cluster merge
im_ls <- list('Myeloid'= myeloid, 'Th_cells'=th_cells, 'Tc_cells'=tc_cells, 
    'T_cells'= t_cells, 'NK_cells'=nk_cells, 'Other_immune'=other_immune)

im_k <- data.frame(cbind(cluster=cl_immune, subtype=rep(NA, length(cl_immune))))
for (i in 1:length(im_ls)){
  k <- which(ct_immune_m$cluster %in% im_ls[[i]])
  im_k[k,'subtype'] <- names(im_ls)[i]
}

# Translate to ct_immune
bdf_immune <- cbind(ct_immune, subtype=im_k$subtype)

# Check general abundances (plot)
table(bdf$cell_type)/nrow(bdf)*100

# Add 'subtypes' column to the DF

# Clustering for Epithelial cells

# Clustering for Endothelial cells

# Clustering for Mesenchymal cells

# Rbind all in one big df


