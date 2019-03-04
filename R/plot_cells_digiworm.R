#' Plot Cells in Digiworm
#'
#' Plots specified cells across digiworm subsets (so far neuronal only)
#' @param cells cells names (could be retrieved using \code{Seurat::Cells()})
#' @param alpha_neg transparency value for negative points.
#' @keywords digiworm cells
#' @export
#' @examples
#' plot_cells_digiworm()

plot_cells_digiworm = function(cells,alpha_neg = 0.1){
  
  plot_subset = function(subset,cells,alpha_neg = 0.1){
    plot_data = coordinates[[subset]]
    plot_data$cell = ifelse(rownames(plot_data)%in%cells,1,0)
    
    
    plot_data$transparency = ifelse(plot_data$cell>0,'pos','neg')
    alphas = c('pos' = 1,'neg'= alpha_neg)
    
    
    p = ggplot(plot_data,aes(x = tSNE_1,y = tSNE_2))+
      geom_point(aes(alpha = transparency,col = cell))+
      scale_alpha_manual(values = alphas, guide = F)+
      scale_color_gradient(low = "grey", high = "red",guide = F)+
      theme_bw()+ggtitle(paste0(sum(rownames(plot_data)%in%cells),'/ ',length(rownames(plot_data)),' cells'))
    return(p)
  }
  
  out = lapply(1:4, FUN = function(i) plot_subset(i,cells = cells ))
  
  return(do.call('grid.arrange',out))
}