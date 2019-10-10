#' Merge Bismark coverage files
#' This function merge a list of Bismark coverage files by adding methylation values
#' @param list list of bismark files
#' @import GenomicRanges
#' @keywords bismark merge union methylation
#' @export
#' @return A dataframe with merged genomic locations and added methylation values
#' @examples
#' concatBismark(list(bedA, bedB, bedC))



concatBismark <- function(listbismark=list){
  bed1 <- data.table::data.table(
    CHROM = listbismark[[1]]$V1,
    START = listbismark[[1]]$V2,
    STOP = listbismark[[1]]$V3,
    PERC = listbismark[[1]]$V4,
    MET= listbismark[[1]]$V5,
    UNMET= listbismark[[1]]$V6,
    key = c("CHROM","START","STOP", "PERC", "MET", "UNMET")
  )
  tmp<-bed1
  for (l in 2:length(listbismark)){
    bed2 <- data.table::data.table(
      CHROM = listbismark[[l]]$V1,
      START = listbismark[[l]]$V2,
      STOP = listbismark[[l]]$V3,
      PERC = listbismark[[l]]$V4,
      MET= listbismark[[l]]$V5,
      UNMET= listbismark[[l]]$V6,
      key = c("CHROM","START","STOP", "PERC", "MET", "UNMET")
    )
    tmp<-intersectBedFiles.GR(tmp, bed2)
  }
  return(tmp)
}

intersectBedFiles.GR <- function(bed1,bed2) {
  bed1 <- GenomicRanges::makeGRangesFromDataFrame(bed1,
                                   seqnames.field = "CHROM",start.field="START",end.field="STOP", keep.extra.columns=TRUE, ignore.strand=TRUE,
                                   seqinfo=NULL)
  bed2 <- GenomicRanges::makeGRangesFromDataFrame(bed2,
                                   seqnames.field = "CHROM",start.field="START",end.field="STOP", keep.extra.columns=TRUE, ignore.strand=TRUE,
                                   seqinfo=NULL)
  grUnion <- suppressWarnings(GenomicRanges::union(bed1,bed2))
  bed1<-GenomicRanges::merge(grUnion, bed1, all.x = TRUE)
  bed2<-GenomicRanges::merge(grUnion, bed2, all.x = TRUE)
  bed1[is.na(bed1$MET)]$MET<-0
  bed2[is.na(bed2$MET)]$MET<-0
  bed1[is.na(bed1$UNMET)]$UNMET<-0
  bed2[is.na(bed2$UNMET)]$UNMET<-0
  TMET=bed1$MET+bed2$MET
  TUNMET=bed1$UNMET+bed2$UNMET

  if(TMET+TUNMET==0){
    percent<-0
  }else{
    percent<-TMET/(TMET+TUNMET)*100
  }


  bedTot<-data.table::data.table(CHROM=as.character(seqnames(bed1)),START=start(bed1),STOP=start(bed1), PERC=percent , MET=TMET, UNMET=TUNMET )
  return(bedTot)
}
