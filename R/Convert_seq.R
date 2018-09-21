#' Convert_seq
#'
#'You can extract the targe sequence and convert it into dataframe
#'
#' @param Chr Chromosome information
#' @param Target.Start Starting location of target sequence
#' @param Target.End ending ocation of target sequence
#' @param Strand +/- strand
#' @param chromosome chromosome information
#' e.g chromosome <- BSgenome.Hsapiens.UCSC.hg19
#' @param NO. name of the each target sequence;default is NULL
#' @return a dataframe
#' @export
#'
#' @examples
#' Chr<-c("chrX","chrX" ,"chrX","chrX" )
#' Target.Start<-c(153996577,153999034,154004462,154002877)
#' Target.End<-c(153996707 ,153999154,154004599,154002980)
#' Strand<-rep("+",4)
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19
#' Convert_seq(Chr,Target.Start,Target.End,Strand,chromosome)
#'
Convert_seq<-function(Chr,Target.Start,Target.End,Strand,chromosome,NO.=NA){
  gr.exome <- GenomicRanges::GRanges(seqnames= Chr,
                                     ranges = IRanges::IRanges(start = Target.Start,
                                                      end = Target.End),
                                     strand = Strand)
  extract_seq <- BSgenome::getSeq(chromosome, gr.exome)
  target.seq<-mapply(as.character,extract_seq )
  x<-data.frame(NO.,seqnames=Chr,
                range=paste0(Target.Start,"-",Target.End),strand=Strand,
                Target.Seq=target.seq)
  if(is.na(NO.)==TRUE){y<-x[,2:5]} else{y<-x}
  return(y)
}
