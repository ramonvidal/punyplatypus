#' Next-neighbors in different clusters
#' @param seurat_object Dataframe with seurat marker predictions.
#' @param top_nn use mouse or human (default: "mouse")
#' @keywords single-cell next-neighbors clusters
#' @export
#' @return array with cells sharing neighbors in different clusters
#' @examples
#' nn_bw_clusters(pbmc, top_nn = 10)



nn_bw_clusters <- function(seurat_object, top_nn=10){
  cells_high<-""
  for (i in seq(1:nrow(seurat_object@graphs$SCT_snn))){
    Snn<-names(head(sort(seurat_object@graphs$SCT_snn[i,], decreasing = T), top_nn))
    if (all(seurat_object@meta.data[Snn,]$seurat_clusters == seurat_object@meta.data[Snn,]$seurat_clusters[1] )==FALSE){
      cells_high<-c(cells_high, Snn[1])
    }
  }
  return(cells_high)
}

