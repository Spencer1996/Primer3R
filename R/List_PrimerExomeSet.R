#' List_PrimerExomeSet
#'
#'This function is used to integrate sheared fragements and corresponding primers into list formation
#'
#' @param split.exome is the list converted from Exome Primer3 result by 'Spliy_exome' function
#'Thereinto, GSP1 is the gene specifical primer1, which is target primer and far from extension
#' sequence；the other closer primer is GSP2
#' @return list
#' @export
#'
#' @examples
#' sp.ex<-Split_exome(exo.out,frame=exo)
#' x<-List_PrimerExomeSet(sp.ex)
#' sapply(x,function(x){c(names(x[1]),names(x[[1]]))})
#' only for one target
#' y<-List_PrimerSet(sp.ex[[1]])
#' sapply(y,names)
List_PrimerExomeSet<-function(split.exome){
  List_primer3Output<-function(x){
    x<-as.data.frame(x)
    primer.list<-list()
    whic<-which(x[,1]=="=")
    for(i in 1:length(whic)){
      if(length(whic)<=1){primer.list[[i]]<-as.data.frame(x)} else{
        primer.list[[i]]<- na.omit(as.data.frame(x[ifelse(i==1,1,whic[i-1]+1):whic[i],]))}
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
eplist<-list()
for (i in 1:length(split.exome)) {
  eplist[[i]]<-List_PrimerSet(List_primer3Output(split.exome[[i]]))}
names(eplist)<-names(split.exome)
return(eplist)}
