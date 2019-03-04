#' Plot Cells
#'
#' This function plots the specific cells.
#' @param seurat_object the object to retrieve the data from.
#' @param cells cells names (could be retrieved using \code{Seurat::Cells()})
#' @param alpha_neg transparency value for negative points.
#' @param reduction dimensionality reduction, could be \code{c('pca','tsne','umap')}
#' @keywords positive cells
#' @export
#' @examples
#' plot_positives()


plot_cells = function(seurat_object,cells,alpha_neg = 0.1,reduction = 'umap'){
  
  reductions = list('pca' = c('PC_1','PC_2'),'tsne' = c('tSNE_1','tSNE_2'),"umap" = c('UMAP_1','UMAP_2'))
  
  
  plot_data = FetchData(seurat_object,var = c(unlist(reductions[reduction])))
  plot_data$cell = ifelse(rownames(plot_data) %in% cells,1,0)
  plot_data$transparency = ifelse(plot_data$cell>0,'pos','neg')
  
  # set transparency values:
  alphas = c('pos' = 1,'neg'= alpha_neg)
  
  p = ggplot(plot_data,aes_string(x = colnames(plot_data)[1],y = colnames(plot_data)[2]))+
    geom_point(aes(alpha = transparency,col = cell))+
    scale_alpha_manual(values = alphas, guide = F)+
    scale_color_gradient(low = "grey", high = "red",guide = F)+
    theme_bw()+ggtitle(paste0(length(cells),' cells'))
  return(p)
}