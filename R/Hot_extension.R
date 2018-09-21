#' Hot_extension
#'
#'For the target sequecne dataframe, if you do not need to shear, you can directly use it to
#'extend the target sequence for designing the primer. Finally, amplified fragment will return
#'in the target sequecne dataframe, called Extend.fragment. Meanwhile, it additionally show the complete
#'sequence after amplifying, called Extend.seq
#' @param HOTSPOT is dataframe converted by Convert_seq funciton
#' @param n total lenght of target sequence after amplifying
#' @param chromosome is the chromosome information for amplified fragment
#' e.g chromosome <- BSgenome.Hsapiens.UCSC.hg19
#'
#' @return data.fame
#' @export
#'
#' @examples
#' Chr<-c("chrX","chrX" ,"chrX","chrX" )
#' Target.Start<-c(153996577,153999034,154004462,154002877)
#' Target.End<-c(153996707 ,153999154,154004599,154002980)
#' Strand<-rep("+",4)
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19
#' frame<-Convert_seq(Chr,Target.Start,Target.End,Strand,chromosome)
#' class<-Class_frame(frame,125)
#' sapply(class,nrow)
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' entend.hot<-Hot_extension(class$HOTSPOT,185,chromosome)
#' View(entend.hot)
Hot_extension<-function(HOTSPOT,n,chromosome){ #x=frame.eh$HOTSPOT}, n is length of extension
  ex_range<-function(x,n){
    y<-vector()
    for(i in 1:length(x)){
      y[i]<-as.numeric(unlist(strsplit(as.character(x[i]),"-")))[n]}
    return(y)
  } # get range
  gr.hotspot <- GenomicRanges::GRanges(seqnames= HOTSPOT$seqnames,
                        ranges = IRanges::IRanges(start = ex_range(HOTSPOT$range,1),
                                         end = ex_range(HOTSPOT$range,2)),
                        strand = HOTSPOT$strand)
  fix_resize<-function(gr.hotspot,n){
    TF<-vector()
    for(i in 1:length(gr.hotspot)){
      TF[i]<-ifelse(as.character(gr.hotspot@strand)[i]=="+","start","end")}
    y<-GenomicRanges::resize(gr.hotspot, n,fix=TF)
    return(y)
  } #extend the range
  ex.gr.hot<-fix_resize(gr.hotspot,n)
  extension_seq <- BSgenome::getSeq(chromosome, ex.gr.hot)
  extend<- GenomicRanges::GRanges(seqnames= HOTSPOT$seqnames,
                   ranges = IRanges::IRanges(start =  BSgenome::end(gr.hotspot),
                                    end =  BSgenome::end(ex.gr.hot)),strand = HOTSPOT$strand)
  Extend.fragment<-BSgenome::getSeq(chromosome, extend)
  HOTSPOT$Extend.fragment<-as.character(Extend.fragment)
  HOTSPOT$Extend.seq<-as.character(extension_seq)
  return(HOTSPOT)}
