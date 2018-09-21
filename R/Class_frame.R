#' Class_frame
#'
#'Based on dataframe converted by Convert_seq funciton, this function can
#'divide the target seqence into Hotspot and Exome two list by thier sequence length
#' @param dataframe.seq is the dataframe converted by Convert_seq funciton
#' @param n is the lenght of sequence, if the target sequecne is less than
#' n, it will belong to Hotspot, otherwise is Exome
#'
#' @return list
#' @export
#'
#' @examples
#' Chr<-c("chrX","chrX" ,"chrX","chrX" )
#' Target.Start<-c(153996577,153999034,154004462,154002877)
#' Target.End<-c(153996707 ,153999154,154004599,154002980)
#' Strand<-rep("+",4)
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' frame<-Convert_seq(Chr,Target.Start,Target.End,Strand,chromosome)
#' class<-Class_frame(frame,125)
#' sapply(class,nrow)
Class_frame<-function(dataframe.seq,n){
  data.f<-list()
  con.list<-function(x){
    v_l<-function(x){
      lis<-list()
      for(i in 1:length(x)){lis[[i]]<-x[i]}
      return(lis)
    }
    if(class(x)=="list"){io=x}else{io=v_l(x)}
    return(lapply(io,as.character))}
  ex_width<-function(x){
    y=s=b=vector()
    for(i in 1:length(x)){
      s[i]<-as.numeric(unlist(strsplit(as.character(x[i]),"-")))[1]
      b[i]<-as.numeric(unlist(strsplit(as.character(x[i]),"-")))[2]
      y[i]<-1+b[i]-s[i]
    }
    return(y)
  }
  io<-con.list(dataframe.seq$range)
  data.f$HOTSPOT<-as.data.frame(dataframe.seq[which(ex_width(io)<=n) ,])
  data.f$EXOME<-as.data.frame(dataframe.seq[which(ex_width(io)>n),])
  return(data.f)
}
