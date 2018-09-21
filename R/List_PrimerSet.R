#' List_PrimerSet
#'
#' @param x is the output result of Primer3
#'
#' @return list
#' @export
#'
#' @examples
#' list<-List_PrimerSet(hot.out)
#' sapply(list,names)
List_PrimerSet<-function(x){
  List_primer3Output<-function(x){
    x<-as.data.frame(x)
    primer.list<-list()
    whic<-which(x[,1]=="=")
    for(i in 1:length(whic)){
      if(length(whic)<=1){primer.list[[i]]<-as.data.frame(x)} else{
        primer.list[[i]]<-as.data.frame(x[ifelse(i==1,1,whic[i-1]+1):whic[i],])%>% na.omit()}
    }
    Na<- x[grep("SEQUENCE_ID=",x[,]),]
    Nam<-sapply(strsplit(as.character(Na),"="),function(x){x[[2]]})
    names(primer.list)<-Nam
    return(primer.list)
  }
  List_PrimerSet<-function(PRIMER.list){
    primer3_list<-function(primer3){
      number<-function(primer3){
        primer3=as.data.frame(primer3)
        firs<-grep("PRIMER_RIGHT_0",primer3[,1])
        if(sum(firs)==0){firs=grep("PRIMER_LEFT_0",primer3[,1])}else{firs=firs}
        end<-grep("=",primer3[,1])
        n<-(end[length(end)]-firs[1]) /length(firs)
        return(n)
      }
      plist<-list()
      if(is.na(number(primer3))==TRUE) {NA}else{
        for (i in 1: number(primer3)){
          plist[[i]]<-primer3[grep(paste0("PRIMER_RIGHT_",i-1),primer3[,1]),1]
        }
        names(plist)<-paste0("Primer",c(1:number(primer3)))}
      return(plist)} #primer list
    plist<-list()
    for (i in 1:length(PRIMER.list)){
      plist[[i]]<-primer3_list(PRIMER.list[[i]])}
    names(plist)<-names(PRIMER.list)
    plist[which(sapply(plist,length)==0)]<-"NO primer"
    return(plist)
  }
  y<-List_PrimerSet(List_primer3Output(x))
  return(y)}
