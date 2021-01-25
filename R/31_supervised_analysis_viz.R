# Visualization

dendrogram_barplot <- function(data, dist_method = 'euclidean', 
    hclust_method = 'ward.D', corr_mat = FALSE, coul, xlabl="Cell type (% of sample)",
    ylabl="Patient ID", legend_p ='none', BPOrderAsDendrogram=T, bp_order,
    cex_names=0.75, cex_axis=0.75, cex_lab=1, cex_sub=1){
    
    # Dendrogram
    if(corr_mat){
        corr_pt <- Hmisc::rcorr(t(as.matrix(data)), type = 'spearman')
        dist_mat <- dist(corr_pt$r, method = dist_method)
    } else {
        dist_mat <- dist(data, method = dist_method)
    }
    hclust_avg <- hclust(dist_mat, method = hclust_method)
    par(mar=c(2,7,4,2), lwd=2, mgp=c(3,1,0), las=1, font.lab=2)
    plot(hclust_avg,cex = 0.8, hang = -1, main = paste0(hclust_method, ' linkage'))
    
    # Barplot
    #coul = coul
    if(BPOrderAsDendrogram){
      data <- data[hclust_avg$order,] #dendrogram order
    }else{
      data <- data[bp_order,]
    }
    
    data <- t(as.matrix(data))
    
    par(mgp=c(1.5,0.25,0), las=1, mar = c(3, 5, 1, 7.5), xpd=TRUE, font.lab=2, lwd=1) # axis label locations
    barplot(data,
            col=coul ,
            border='white',
            horiz=TRUE,
            cex.names = cex_names,
            cex.axis = cex_axis,
            cex.lab = cex_lab,
            #cex.sub = cex_sub,
            legend = FALSE)
    title(xlab=xlabl, mgp = c(1.5, 0.5, 0))
    title(ylab=ylabl,  mgp = c(3, 0.1, 0))
    # legend('topright', legend = rownames(data), fill = coul, ncol = 1, inset=c(-0.3,0),
    #            cex = 0.75) 

}


ClusterUMAP_plot <- function(data, 
                             density = F,
                             color_by_cluster = T, 
                             color_by_protein = F,
                             color_by_continuous = F, 
                             cluster_col, ft_cols, 
                             pallete = F,
                             plot_clusnames = F,
                             plot_colors,
                             plot_title,
                             lg_names,
                             x_lim = c(-10,10),
                             y_lim = c(-10,10),
                             nrow_prot = 2,
                             heights_prot = unit(c(1.3,1.3), c("in", "in")),
                             title_size_prot = 10, epi_clus=F){
  
  if(density){
    test <- data[c('UMAP1', 'UMAP2')]
    p_cl <- ggplot(test, aes(x = UMAP1, y = UMAP2)) + 
      geom_point(shape = 20, size = 0.01) + 
      ggtitle(plot_title) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line =  element_blank(),
            legend.position = "right",
            legend.title=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
            legend.background = element_blank(),
            legend.spacing.y = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 8),
            plot.title = element_text(size=16, face="bold"),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5,"cm"),
            axis.title=element_text(size=10)) +
      stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", n=100) +
      #xlim(x_lim) +
      #ylim(y_lim) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
      scale_fill_viridis_c(option = 'magma') #https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
    return(p_cl)
  }
  
  if(color_by_cluster){
    test <- data[c(cluster_col, 'UMAP1', 'UMAP2')]
    test[,cluster_col] <- as.factor(test[,cluster_col])
    
    p_cl <- ggplot(test) + 
      geom_point(aes_string(x='UMAP1', y='UMAP2', colour = test[,cluster_col]), 
                 shape = 19, size = 0.01, alpha = 0.4) +
      ggtitle(plot_title) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line =  element_blank(),
            legend.position = "right",
            legend.title=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
            legend.background = element_blank(),
            legend.spacing.y = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 8),
            plot.title = element_text(size=16, face="bold"),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5,"cm"),
            axis.title=element_text(size=10)) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
      guides(colour = guide_legend(override.aes = list(size=3, alpha=1, pch=15))) #override legend features
    
    if(pallete){
      #p_cl <- p_cl + scale_color_brewer(palette="Set3", labels = lg_names) + scale_color_hue(l=60, c=65)
      p_cl <- p_cl + scale_color_jcolors(palette = 'pal8', labels = lg_names)
    } else {
      p_cl <- p_cl + scale_color_manual(values=plot_colors, labels = lg_names)
    }
    
    if (plot_clusnames){
      edata <- cbind(test$UMAP1, test$UMAP2, cluster_col=test[cluster_col])
      colnames(edata) <- c('x', "y", "z")
      center <- aggregate(cbind(x,y) ~ z, data = edata, median)
      if(epi_clus){
        lb <- sapply(strsplit(as.character(center[,1]), "_"), "[[", 2)
      }else{
          lb <- center[,1]
        }
      p_cl <- p_cl + annotate("text", label = lb, 
                              x=center[,2], y = center[,3], size = 6, colour = "black", fontface = 'bold')
    }
    return(p_cl)
  }
  
  if(color_by_protein){
    data2 <- data[,ft_cols]
    data <- cbind(data2, data[c('UMAP1', 'UMAP2')])
    #data <- data[c(colnames(data)[ft_cols], 'UMAP1', 'UMAP2')]
    test <- melt(data, id= c('UMAP1', 'UMAP2'))
    var_list <- unique(test$variable)
    
    #x <- sapply(strsplit(var_list, "_"), "[[", 2)
    pl <- list()
    for(i in seq_along(var_list)) {
      p <- ggplot(subset(test,test$variable==var_list[i]), aes(x=UMAP1, y=UMAP2, colour = value)) + 
        geom_point(shape = 20, alpha=0.4, size = 0.5) + 
        theme_bw() + 
        ggtitle(sapply(strsplit(as.character(var_list[i]), "_"), "[[", 2)) + 
        scale_colour_gradient(low = "grey72", high = "blue") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.background = element_blank(), axis.line =  element_blank(),
              #axis.ticks=element_blank(), 
              axis.title=element_blank(),
              axis.text=element_blank(), aspect.ratio = 1,
              legend.position = "none",
              plot.title = element_text(size=title_size_prot)) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0))
      
      pl[[i]] <- p
    }
    umap_pl <- grid.arrange(grobs=pl, nrow=nrow_prot, heights=heights_prot)
    return(umap_pl)
  }

  if(color_by_continuous){
    test <- data[c(cluster_col, 'UMAP1', 'UMAP2')]
    
    p_cl <- ggplot(test) + 
      geom_point(aes_string(x='UMAP1', y='UMAP2', colour = test[,cluster_col]), 
                 shape = 20, size = 0.2, alpha = 0.4) +
      ggtitle(plot_title) +
      theme_bw()+
      scale_colour_gradient2(low = "#3498DB", mid= 'white', high = "#EC7063", midpoint = 0.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line =  element_blank(),
            legend.position = "right",
            legend.title=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
            legend.background = element_blank(),
            legend.spacing.y = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 8),
            plot.title = element_text(size=16, face="bold"),
            legend.key.size = unit(0.5, "cm"),
            legend.key.width = unit(0.5,"cm"),
            axis.title=element_text(size=10)) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
      scale_y_continuous(limits = y_lim, expand = c(0, 0)) 
  
  return(p_cl)
  }  
}

