#' Cluster Data
#'
#' This function computes the mean expression and number of positive cells for markers in the Seurat object.
#' @param seurat_object the object to retrieve the data from.
#' @param markers markers .
#' @keywords cluster
#' @export
#' @examples
#' cluster_data()

cluster_data = function(seurat_object,markers){
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
  return(cluster_data)
}