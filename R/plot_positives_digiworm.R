#' Plot Multiple Positives in Digiworm
#'
#' This function plots the cells that are positive at multiple markers at the same time using tSNE coordinatess from digiworm.
#' @param seurat_object the object to retrieve the cell data from ( cells that don't appear in this object are not plotted!)
#' @param markers positive markers of the cells of interest.
#' @param marker_names optional names for markers
#' @param alpha_neg transparency value for negative points.
#' @keywords positive cells
#' @export
#' @examples
#' plot_positives_digiworm()

plot_positives_digiworm = function(seurat_object,markers,marker_names = NULL,alpha_neg = 0.1){
  
  cell_data  = FetchData(seurat_object,vars = markers)
  
  plot_subset = function(subset,alpha_neg = alpha_neg){
    plot_data = coordinates[[subset]]
    
    # positive cells 
    positive_cells = rownames(cell_data)[apply(cell_data,1,FUN = function(x) all(x>0))]
    
    plot_data$cell = ifelse(rownames(plot_data)%in%positive_cells,1,0)
    
    
    plot_data$transparency = ifelse(plot_data$cell>0,'pos','neg')
    alphas = c('pos' = 1,'neg'= alpha_neg)
    
    
    p = ggplot(plot_data,aes(x = tSNE_1,y = tSNE_2))+
      geom_point(aes(alpha = transparency,col = cell))+
      scale_alpha_manual(values = alphas, guide = F)+
      scale_color_gradient(low = "grey", high = "red",guide = F)+
      theme_bw()+ggtitle(paste0(sum(rownames(plot_data)%in%positive_cells),'/ ',length(rownames(plot_data)),' cells'))
    return(p)
  }
  
  out = lapply(1:4, FUN = function(i) plot_subset(i))
  
  if(is.null(marker_names)){
    marker_labels = gsub('-','',str_extract(markers,'-[0-9]+-'))
  }else{
    marker_labels = marker_names
  }
  
  return(do.call('grid.arrange',c(out,top=paste0(marker_labels,collapse = '/','+'))))
}
