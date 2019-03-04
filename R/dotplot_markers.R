#' Dot Plot of Markers
#'
#' This function plots the expression of the markers across clusters in the Seurat object. It's analogous to the \code{Seurat::DotPlot()} but it doesn't change the order of markers.
#' If the named list of markers provided, the function use the names as titles for plots.
#' @param seurat_object the object to retrieve the data from.
#' @param markers a genes to plot.
#' @keywords dotplot
#' @export
#' @examples
#' dotplot_markers()

dotplot_markers  = function(seurat_object, markers){
  seurat_subset = FetchData(seurat_object,vars = markers)
  seurat_subset$cluster = Idents(seurat_object)
  cluster_data = lapply(markers, FUN = function(gene) {
    gene_data = t(sapply(levels(seurat_subset$cluster),FUN = function(x) {
      out = rep(NA,3)
      out[1] = mean(subset(seurat_subset,cluster == x)[,gene])
      out[2] = sum(subset(seurat_subset,cluster == x)[,gene]>0)/length(subset(seurat_subset,cluster == x)[,gene]>0)
      out[3] = length(subset(seurat_subset,cluster == x)[,gene]>0)
      return(out)
    }))
    gene_data = as.data.frame(gene_data)
    gene_data$cluster = rownames(gene_data)
    gene_data$gene = gene
    colnames(gene_data)[c(1,2,3)] = c('mean','perc','n_cells')
    return(gene_data)
  })
  cluster_data = do.call('rbind',cluster_data)
  row.names(cluster_data) = NULL
  
  if(!is.null(names(markers))){
    cluster_data$gene = names(markers)[match(cluster_data$gene,markers)]
  }
  
  ggplot(cluster_data,aes(y = gene,x = cluster,group = cluster))+
    geom_point(aes(fill = mean, size = perc),shape = 21)+
    scale_fill_gradient(low = 'red',high = 'yellow')+
    scale_radius(range = c(0,20),breaks = c(0.1,0.5,1),limits = c(0,NA))+
    theme_bw()+
    ggtitle('Cluster comparison')
}