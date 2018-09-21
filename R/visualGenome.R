#' visualGenome
#'
#' @param Chr Chromosome information，e.g "Chr1"
#' @param Target.Start Starting location of target sequence
#' @param Target.End ending location of target sequence
#' @param Extend.length toal length of extended target sequence
#' @param ID Name or NO. of target sequence
#' @param x a pair primer 3 avliable from Hot_output or Exome_outpur function
#' @param plot plot the genomic location, defaulr is to ‘TRUE’
#'
#' @return GRange of target and primer pair information
#' @export
#'
#' @examples
#' Chr<-"chr5"
#' Target.Start<-176939369
#' Target.End<-176939371
#' Extend.length<-185
#' ID<-"leukemia_36"
#' plist<- List_PrimerSet(hot.out)
#' hot.pair<-List_PairSet(plist)
#' rem_Pairs.inf<-Filter_Inf(hot.pair,plist,distance=20,overlap=10)
#' x<-Hot_output(rem_Pairs.inf,Write.out=FALSE)[3,]
#' visualGenome(Chr,Target.Start,Target.End,Extend.length,ID,x)
#' #only return local information
#' visualGenome(Chr,Target.Start,Target.End,Extend.length,ID,x,plot=FALSE)
visualGenome<-function(Chr,Target.Start,Target.End,Extend.length,ID,x,plot=TRUE){
  gr <- GenomicRanges::GRanges(Chr,
                               ranges = IRanges::IRanges(
                                 start = Target.Start,
                                 end = ((Target.End+Extend.length)-1)),ID=ID)
  gsp1.s<-(as.numeric(x$GSP1.End)+Target.Start-1)
  gsp2.s<-(as.numeric(x$GSP2.End)+Target.Start-1)
  range <- GenomicRanges::GRanges(rep(Chr,3),
                               ranges = IRanges::IRanges(
                                 start = c(Target.Start,gsp1.s,gsp2.s),
                                 end = c(Target.End,(as.numeric(x$GSP1.Length)+gsp1.s-1),c(gsp2.s+ as.numeric(x$GSP2.Length)-1))),ID=c(ID,"GSP1","GSP2"))

  ensGeneTrack <- TnT::FeatureTrack(range, tooltip = as.data.frame(pr),
                                    names = c(ID,"GSP1","GSP2"),
                                    color = TnT::mapcol(c(ID,"GSP1","GSP2"), palette.fun = grDevices::rainbow))
  if(plot==TRUE){
    TnT::TnTGenome(ensGeneTrack, view.range = gr)
  }else{return(range)}
}
