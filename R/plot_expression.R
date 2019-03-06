#' Plot Expression
#'
#' This function plots the expression of the feature across cells. It's analogous to the \code{Seurat::FeaturePlot()} but it changes the transparency of negative cells so all positive ones are visible.
#' @param seurat_object the object to retrieve the data from.
#' @param marker a gene to plot.
#' @param alpha_neg transparency value for negative points.
#' @param reduction dimensionality reduction, could be \code{c('pca','tsne','umap')}
#' @keywords expression
#' @export
#' @examples
#' cluster_data()

plot_expression = function(seurat_object,marker,alpha_neg = 0.1,reduction = 'umap'){
  
  reductions = list('pca' = c('PC_1','PC_2'),'tsne' = c('tSNE_1','tSNE_2'),"umap" = c('UMAP_1','UMAP_2'))
  
  
  plot_data = FetchData(seurat_object,var = c(unlist(reductions[reduction]),marker))
  colnames(plot_data)[3] = 'marker'
  plot_data$transparency = ifelse(plot_data$marker>0,'pos','neg')
  
  # set transparency values:
  alphas = c('pos' = 1,'neg'= alpha_neg)
  
  plot_data$marker[plot_data$marker==0] = NA
  
  p = ggplot(plot_data,aes_string(x = colnames(plot_data)[1],y = colnames(plot_data)[2]))+
    geom_point(aes(alpha = transparency,col = marker))+
    scale_alpha_manual(values = alphas, guide = F)+
    scale_color_gradient(na.value = 'grey',low = "red", high = "yellow")+
    theme_bw()+
    ggtitle(paste0(marker))
  return(p)
}
# dimenshionality reduction plot - the default Seurat function doesn't change the transparency of the points based on the expression - many 
# positive cells are overlapped by negative ones