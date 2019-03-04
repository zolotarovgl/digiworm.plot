#' Plot tSNE from digiworm
#'
#' This function plots a tSNE image from digiworm for a given marker
#' @param subset the subset of cells from digiworm: 1 - all clusters; 2 - neural; 3 - ciliated; 4 - non-ciliated.
#' @param marker the marker to plot.
#' @param alpha_neg transparency value for negative points.
#' @keywords positive cells
#' @export
#' @examples
#' plot_tsne()

plot_tsne = function(subset, marker,alpha_neg){
  plot_data = coordinates[[subset]]
  gene = FetchData(all_clusters_small, marker ,cells = rownames(plot_data))
  plot_data = merge(plot_data,gene, by  = 0, all = T)
  colnames(plot_data)[4] = 'marker'
  summary(plot_data$marker)
  plot_data$transparency = ifelse(plot_data$marker>0,'pos','neg')
  alphas = c('pos' = 1,'neg'=alpha_neg)
  plot_data$marker = as.numeric(plot_data$marker)
  plot_data$marker[plot_data$marker<=0]=NA
  p = ggplot(plot_data,aes(x = tSNE_1,y = tSNE_2))+
    geom_point(aes(alpha = transparency,col = marker))+scale_alpha_manual(values = alphas,guide = F)+
    scale_color_gradient(na.value = 'grey',low = "red", high = "yellow")
  return(p)
}

#' Plot tSNE from digiworm
#'
#' This function plots a tSNE image for all subsets of cells in digiworm.
#' @param marker the marker to plot.
#' @param alpha_neg transparency value for negative points.
#' @keywords positive cells
#' @export
#' @examples
#' plot_digiworm()

plot_digiworm = function(marker,alpha_neg){
  p1 = plot_tsne(1,marker,alpha_neg)+ggtitle(paste0(marker,' all'))
  p2 = plot_tsne(2,marker,alpha_neg)+ggtitle(paste0(marker,' all_neurons'))
  p3 = plot_tsne(3,marker,alpha_neg)+ggtitle(paste0(marker,' ciliated'))
  p4 = plot_tsne(4,marker,alpha_neg)+ggtitle(paste0(marker,' non-ciliated'))
  grid.arrange(p1,p2,p3,p4, nrow =2)
}