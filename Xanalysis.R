#################################
# k-means clustering

# Reading file
file_name <- '180321_2_8356_G_normalized_denoised.FCS' 
smp <- as.data.frame(flowCore::read.FCS(file_name, transformation = FALSE)@exprs)
c_cols <- c(15, 20, 29, 31, 37, 40)

# tsne
scaled_smp <- denoisingCTF::t_asinh(smp[,c_cols])
tsne_smp <- Rtsne::Rtsne(scaled_smp,  check_duplicates = FALSE)

# clustering
n_clusters <- 3
cl_data <- stats::kmeans(scaled_smp, centers = n_clusters)$cluster

# Editing df
scaled_smp_cl <- as.data.frame(cbind(class = as.factor(cl_data),
                                     scaled_smp,
                                     tsne1 = tsne_smp$Y[,1],
                                     tsne2 = tsne_smp$Y[,2]))
vars <- c('EpCAM', 'CD31', 'CD45', 'Vimentin', 'CK', 'CK7')
cn <- c('class', vars, 'tsne1', 'tsne2')
colnames(scaled_smp_cl) <- cn

# Plotting
library(ggplot2)
library(gridExtra)
library(reshape2)
test <- melt(scaled_smp_cl, id= c('class', 'tsne1', 'tsne2'))
var_list <- unique(test$variable)
pl <- list()
p_cl <- ggplot(scaled_smp_cl, aes(x=tsne1, y=tsne2, colour = class))+
  geom_point(alpha=0.3) + theme_bw() + ggtitle('class')
pl[[1]] <- p_cl
for(i in seq_along(var_list)) {
  p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=tsne1, y=tsne2, colour = value)) +
    geom_point(alpha=0.3) + theme_bw() + ggtitle(var_list[i]) +
    scale_colour_gradient(low = "gray75", high = "blue")
  pl[[i+1]] <- p
}

grid.arrange(grobs=pl)

#################################

file_name <- '180321_2_8356_G_normalized_beads_neg.fcs' 
smp <- as.data.frame(flowCore::read.FCS(file_name, transformation = FALSE)@exprs)
c_cols <- c(15, 20, 29, 31, 37, 40)

# tsne
scaled_smp <- denoisingCTF::t_asinh(smp[,c_cols])
tsne_smp <- Rtsne::Rtsne(scaled_smp,  check_duplicates = FALSE)

# clustering
n_clusters <- 3
cl_data <- stats::kmeans(scaled_smp, centers = n_clusters)$cluster

# Editing df
scaled_smp_cl <- as.data.frame(cbind(class = as.factor(cl_data),
                                     scaled_smp,
                                     tsne1 = tsne_smp$Y[,1],
                                     tsne2 = tsne_smp$Y[,2]))
vars <- c('EpCAM', 'CD31', 'CD45', 'Vimentin', 'CK', 'CK7')
cn <- c('class', vars, 'tsne1', 'tsne2')
colnames(scaled_smp_cl) <- cn

# Plotting
library(ggplot2)
library(gridExtra)
library(reshape2)
test <- melt(scaled_smp_cl, id= c('class', 'tsne1', 'tsne2'))
var_list <- unique(test$variable)
pl <- list()
p_cl <- ggplot(scaled_smp_cl, aes(x=tsne1, y=tsne2, colour = class))+
  geom_point(alpha=0.3) + theme_bw() + ggtitle('class')
pl[[1]] <- p_cl
for(i in seq_along(var_list)) {
  p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=tsne1, y=tsne2, colour = value)) +
    geom_point(alpha=0.3) + theme_bw() + ggtitle(var_list[i]) +
    scale_colour_gradient(low = "gray75", high = "blue")
  pl[[i+1]] <- p
}

grid.arrange(grobs=pl)

