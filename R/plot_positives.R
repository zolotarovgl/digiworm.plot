#' Plot Multiple Positives
#'
#' This function plots the cells that are positive at multiple markers at the same time.
#' @param seurat_object the object to retrieve the data from.
#' @param markers positive markers of the cells of interest.
#' @param marker_names optional names for markers
#' @param alpha_neg transparency value for negative points.
#' @param reduction dimensionality reduction, could be \code{c('pca','tsne','umap')}
#' @keywords positive cells
#' @export
#' @examples
#' plot_positives()

plot_positives = function(seurat_object,markers,marker_names=NULL,alpha_neg = 0.1,reduction = 'umap'){
  
  reductions = list('pca' = c('PC_1','PC_2'),'tsne' = c('tSNE_1','tSNE_2'),"umap" = c('UMAP_1','UMAP_2'))
  
  
  plot_data = FetchData(seurat_object,var = c(unlist(reductions[reduction]),markers))
  colnames(plot_data)[3:length(colnames(plot_data))] = paste0('marker',1:(length(colnames(plot_data))-2))
  plot_data$col = ifelse(apply(plot_data[,grepl('marker',colnames(plot_data))],1,FUN = function(x) all(x>0)),1,0)
  plot_data$transparency = ifelse(plot_data$col>0,'pos','neg')
  
  if(is.null(marker_names)){
    marker_labels = gsub('-','',str_extract(markers,'-[0-9]+-'))
  }else{
    marker_labels = marker_names
  }
  # set transparency values:
  alphas = c('pos' = 1,'neg'= alpha_neg)
  
  p = ggplot(plot_data,aes_string(x = colnames(plot_data)[1],y = colnames(plot_data)[2]))+
    geom_point(aes(alpha = transparency,col = col))+
    scale_alpha_manual(values = alphas, guide = F)+
    scale_color_gradient(low = "grey", high = "red")+
    theme_bw()+
    labs(subtitle = paste0(sum(plot_data$col),' cells'))+
    ggtitle(paste0(marker_labels,collapse = '/','+'))
  return(p)
} 