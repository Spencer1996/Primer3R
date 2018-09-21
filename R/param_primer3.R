#' param_primer3
#'
#'Due to there are many parameters required by Primer3, you can specifically supplement
#'the parameters you expected
#'
#' @param x is the Priemr3 input data
#' @param add supplementary parameters for Primer3
#'
#' @return
#' @export
#'
#' @examples
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19
#' entend.hot<-Hot_extension(hot,185,chromosome)
#' p.hot<-Primer3_Hotset(entend.hot,return=2)
#' View(p.hot)
#' add<-c("PRIMER_INTERNAL_MAX_TM=63","PRIMER_MAX_END_GC=5")
#' param_primer3(p.hot,add)

param_primer3<-function(x,add){
  List_primer3Output<-function(x){
    x<-as.data.frame(x)
    primer.list<-list()
    whic<-which(x[,1]=="=")
    for(i in 1:length(whic)){
      if(length(whic)<=1){primer.list[[i]]<-as.data.frame(x)} else{
        primer.list[[i]]<-na.omit(as.data.frame(x[ifelse(i==1,1,whic[i-1]+1):whic[i],]))}
    }
    Na<- x[grep("SEQUENCE_ID=",x[,]),]
    Nam<-sapply(strsplit(as.character(Na),"="),function(x){x[[2]]})
    names(primer.list)<-Nam
    return(primer.list)
  }
add_primer3<-function(x,add){
  y<-rbind(as.matrix(x[-nrow(x),1]),as.matrix(add))
  z<-rbind(y,as.matrix("="))
  return(z)
}
x<-List_primer3Output(x)
g<-list()
for(i in 1:length(x)){
  g[[i]]<-add_primer3(x[[i]],add)
}
gg<-rlist::list.rbind(g)
return(gg)}

