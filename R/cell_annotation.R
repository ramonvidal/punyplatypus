#' Cell type guess
#' This function is a R adaptation from Alona similar python function. It uses PanglaoDB database
#' @param markers Dataframe with seurat marker predictions.
#' @param species use mouse or human. Default is mouse
#' @keywords single-cell celltype
#' @export
#' @import dplyr
#' @return A list containing celltype predictions with p-values
#' @examples
#' celltypes(markers, species = "mouse")


#TODO :
# Make p-value a parameter
# remove the dplyr library call from here
#

#library(dplyr)

cellType <- function(markers, species="mouse"){
  markers<-markers[, c(2,6,7)]
  markers$gene<-toupper(markers$gene)
  split_list<- split(markers,markers$cluster)
  all_genes_dataset<-unique(markers$gene)
  markerSet<-markersSet(species=species)
  markerSet_freq<-table(markerSet$official.gene.symbol)
  weights = as.data.frame(1+sqrt(((max(markerSet_freq)-markerSet_freq)/(max(markerSet_freq)-min(markerSet_freq)))))


  split_list_marker<- split(markerSet,markerSet$cell.type)
  all_genes_markers<-unique(markerSet$official.gene.symbol)

  output<-c("", "", 0.1)
  output<-as.data.frame(output)
  colnames(output)<-"NA"
  rownames(output)<-c("cluster", "celltype", "pvalue")

  for (i in 1:length(split_list)){
    #results=NA
    #names(results)<-"None"
    cluster<-as.data.frame(split_list[[i]])
    cluster<-merge(cluster, weights[weights$Var1 %in% cluster$gene,], by.x="gene", by.y="Var1")
    for (j in 1:length(split_list_marker)){
      markers_cellT<-as.data.frame(split_list_marker[[j]])
      markers_cellT<-merge(markers_cellT, weights[weights$Var1 %in% markers_cellT$official.gene.symbol,], by.x="official.gene.symbol", by.y="Var1")
      #how many expressed genesets are found in the geneset?
      ct_exp<-(intersect(cluster$gene,markers_cellT$official.gene.symbol))
      # how many _non_ expressed genes are found in the geneset?
      diff_gene_dataset<-setdiff(all_genes_dataset, ct_exp)
      ct_non_exp<-(intersect(diff_gene_dataset,markers_cellT$official.gene.symbol))
      # how many expressed genes are NOT found in the geneset?
      ct_exp_not_found<-setdiff(cluster$gene,markers_cellT$official.gene.symbol)
      # how many _non_ expressed genes are NOT found in the geneset?
      not_exp_not_found_in_geneset<-(setdiff(diff_gene_dataset, markers_cellT$official.gene.symbol))
      odds<-fisher.test(matrix(c(length(ct_exp), length(ct_non_exp), length(ct_exp_not_found), length(not_exp_not_found_in_geneset)),  nrow = 2))
      if(odds$p.value<0.05){ #Make this a parameter
        fish<-odds$p.value
        results<-c(unique(as.character(cluster$cluster)), unique(as.character(markers_cellT$cell.type)), fish)
        results<-as.data.frame(results)
        rownames(results)<-c("cluster", "celltype", "pvalue")
        colnames(results)<-j
        output<-cbind(output, results)
      }
    }
    #print(unique(as.character(cluster$cluster)))
  }
  output<-output[,2:ncol(output)]
  output<-t(output)
  output<-as.data.frame(output)
  output$pvalue<-as.numeric(as.character(output$pvalue))
  output<-output[order(output$cluster, output$pvalue),]
  rownames(output)<-NULL
  return(output)
}

markersSet <- function(species="mouse"){
  ma<-read.csv(system.file("data", "markers.tsv", package = "punyplatypus"), sep="\t")
    if (species == 'mouse'){
      s = 'Mm'
    }else if (species == 'athaliana'){
      s = 'At'
    } else{
      s = 'Hs'
    }
  ma<-ma[grepl(s, ma$species),]
  ma<-ma[ma$ubiquitousness.index<0.05,]
  ma_ss<-ma[,2:3]
  return(ma_ss)
}

