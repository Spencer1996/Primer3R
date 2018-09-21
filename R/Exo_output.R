#' Exo_output
#'
#'This function can write the filtered result exome sequecne into '.CSV' formation. Here, it only chooses the optimal
#'priemr pairs for each target sequence. If there is not primer, it only returns NA result; if  no
#' primer pairs fit the condiditons, it returns "None of them fit" at the first column
#'
#' @param fil.exo is the final filtered result of exome sequecne, namely non-sheared sequence
#' @param exo.ex is the exome targe sequence dataframe after extension processed by
#' 'Exo_extension' funciton
#' @param reserve can reserve the filtered row, default is FALSE .
#' @param Write.out  it can write the result into '.csv',default is FALSE.
#' @param NAME  name of output file with '.CSV' formation.
#'
#' @return
#' @export
#'
#' @examples
#' sp.ex<-Split_exome(exo.out,frame=exo)
#' exomepair.set<-List_PrimerExomeSet(sp.ex) %>% List_PairExomeSet()
#' fil.exo<-Filter_InfExome(exomepair.set,eplist,distance=20,overlap=10)
#' entend.exome<-Exo_extension(exo,280,200,chromosome)
#' Exo_output(fil.exo, entend.exome$Extend.seq)
#'
#' Output:
#' Exo_output(fil.exo, entend.exome$Extend.seq,Write.out=TRUE,NAME="exoResult")
Exo_output<-function(fil.exo,exo.ex,reserve=FALSE,Write.out=FALSE,NAME){
  Ex_first<-function(x){
    lis<-list()
    na<-data.frame(GSP1=NA,GSP2=NA,GSP1.Seq=NA,GSP1.Start=NA,
                   GSP1.End=NA, GSP1.Length=NA,GSP1.TM=NA,GSP2.Seq=NA,GSP2.Start=NA,
                   GSP2.End=NA,GSP2.Length=NA,GSP2.TM=NA)
    cna<-data.frame(GSP1="None of them fit",GSP2=NA,GSP1.Seq=NA,GSP1.Start=NA,
                    GSP1.End=NA, GSP1.Length=NA,GSP1.TM=NA,GSP2.Seq=NA,GSP2.Start=NA,
                    GSP2.End=NA,GSP2.Length=NA,GSP2.TM=NA)
    for(i in 1:length(x)){
      lis[[i]]<-if(class(x[[i]])!="data.frame"){cna}else{
        if(ncol(x[[i]])==1){na}else{
          as.data.frame(x[[i]][1,])}}}
    y<-rlist::list.rbind(lis)
    return(y)
  }
  use.fil<-function(Result,cut){Result[which(sapply(Result,length)==sapply(cut,length))]}
  fil.result<-use.fil(fil.exo,exo.ex)
  x<-lapply(fil.result,Ex_first)
  y<- as.data.frame(rlist::list.rbind(x))
      NO<-data.frame(NO.=rownames(y))
  res<-cbind(NO,y)
  rownames(res)<-NULL
  if(reserve==TRUE){out<-res}else{out<-res[which(res$GSP2 != "NA"),]}
  if(Write.out==TRUE){write.csv(out,paste0(NAME,".csv"),row.names = FALSE)} else{return(out)}
  }
