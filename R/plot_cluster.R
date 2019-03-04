#' Plot Cluster
#'
#' This function plots a desired cluster contained in the seurat object
#' @param seurat_object the object to retrieve the data from.
#' @param clusters cluster or cluster to plot cells from.
#' @param alpha_neg transparency value for negative points.
#' @param reduction dimensionality reduction, could be \code{c('pca','tsne','umap')}
#' @keywords positive cells
#' @export
#' @examples
#' plot_positives()

plot_cluster = function(seurat_object,clusters,alpha_neg = 0.1,reduction = 'umap'){
  
  reductions = list('pca' = c('PC_1','PC_2'),'tsne' = c('tSNE_1','tSNE_2'),"umap" = c('UMAP_1','UMAP_2'))
  
  
  plot_data = FetchData(seurat_object,var = c(unlist(reductions[reduction])))
  plot_data$cluster = as.character(Idents(seurat_object)[rownames(plot_data)])
  plot_data$cluster[!plot_data$cluster%in%clusters] = 'none'
  plot_data$transparency = ifelse(plot_data$cluster=='none','neg','pos')
  
  # set transparency values:
  alphas = c('pos' = 1,'neg'= alpha_neg)
  
  p = ggplot(plot_data,aes_string(x = colnames(plot_data)[1],y = colnames(plot_data)[2]))+
    geom_point(aes(alpha = transparency,col = cluster))+
    scale_alpha_manual(values = alphas, guide = F)+
    theme_bw()+
    ggtitle(paste0('Clusters: ',clusters,collapse = ','))
  return(p)
}