#' Random gene list
#' This function creates random gene lists.
#' @param n_genes How many random genes? Defaults to 10.
#' @import biomaRt
#' @keywords genes
#' @export
#' @return A 2 column table with gene names in the second column (gene).
#' @examples
#' randomGeneList(n_genes=10)


randomGeneList <- function(n_genes=10){
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  res <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol","chromosome_name",'strand','transcript_start','transcript_end'), mart = mouse)
  res<-as.data.frame(sample(res$mgi_symbol, n_genes))
  colnames(res)<-"gene"
  return(res)
}


