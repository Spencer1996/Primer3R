#' Hot_output
#'
#'This function can write the filtered result into '.CSV' formation. Here, it only chooses the optimal
#'priemr pairs for each target sequence. If there is not primer, it only returns NA result; if  no
#' primer pairs fit the condiditons, it returns "None of them fit" at the first column
#'
#'
#' @param fil.hot is the final filtered result of hotspot sequecne, namely non-sheared sequence.
#' @param reserve can reserve the filtered row, default is FALSE .
#' @param Write.out  it can write the result into '.csv',default is FALSE.
#' @param NAME  name of output file with '.CSV' formation.
#' @return
#' @export
#'
#' @examples
#' plist<-List_PrimerSet(hot.out)
#' hot.pair<-List_PairSet(plist)
#' rem_Pairs.inf<-Filter_Inf(hot.pair,plist,distance=20,overlap=10)
#' Hot_output(rem_Pairs.inf)
#'
#' Output:
#' Hot_output(rem_Pairs.inf,Write.out=TRUE,NAME="hotResult")
Hot_output<-function(fil.hot,reserve=FALSE,Write.out=FALSE,NAME){
  result.com<-list()
  na<-data.frame(GSP1=NA,GSP2=NA,GSP1.Seq=NA,GSP1.Start=NA,
                 GSP1.End=NA, GSP1.Length=NA,GSP1.TM=NA,GSP2.Seq=NA,GSP2.Start=NA,
                 GSP2.End=NA,GSP2.Length=NA,GSP2.TM=NA)

  nofit<-data.frame(GSP1="None of them fit",GSP2=NA,GSP1.Seq=NA,GSP1.Start=NA,
                    GSP1.End=NA, GSP1.Length=NA,GSP1.TM=NA,GSP2.Seq=NA,GSP2.Start=NA,
                    GSP2.End=NA,GSP2.Length=NA,GSP2.TM=NA)
  for (i in 1:length(fil.hot)){
    result.com[[i]]<-if(class(fil.hot[[i]])=="character"){nofit}else{
      if(ncol(fil.hot[[i]])==1){na}else{
        as.data.frame(fil.hot[[i]][1,])}}}
  result.frame<-rlist::list.rbind(result.com)
    NO.<-data.frame(NO.=names(fil.hot))
  y<-cbind(NO.,result.frame)
  rownames(y)<-NULL
  if(reserve==TRUE){out<-y}else{out<-y[which(y$GSP2 != "NA"),]}
  if(Write.out==TRUE){write_csv(out,paste0(NAME,".csv"))} else{return(out)}
  }
