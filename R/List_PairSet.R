#' List_PairSet
#'
#'List_PairSet can match primer pair, the high quality primers are ranked at top
#'
#' @param plist  is the list formation of  Primers under each target sequence.Thereinto,
#' GSP1 is the gene specifical primer1, which is target primer and far from extension
#' sequence；the other closer primer is GSP2.
#' @param dir   is the direction of extension of target sequence, if the strand of target sequence
#' is "+", the direction of extension of target sequence will be "left"; The direction of
#' "-" strand is "right". sort order of GSP1 and GSP2 depends on dir' parameter. Here, default is "left".
#' @return list
#' @export
#'
#' @examples
#' plist<-List_PrimerSet(hot.out)
#' List_PairSet(plist)
#' List_PairSet(plist,dir="right")
List_PairSet<-function(plist,dir="left"){
  ListTodata.frame<-function(plist){
    pair<-function(plist,sta){
      m<-matrix(0,length(plist),length(plist))
      rownames(m)=colnames(m)=sta
      for(i in 1:ncol(m)){
        m[,i]<-ifelse(as.numeric(colnames(m)[i])>as.numeric(rownames(m)),1,0)}
      colnames(m)=rownames(m)<-names(plist) #rename after runnning the 0/1 matrix
      g=igraph::graph.adjacency(t(m),mode="directed",weighted=T)
      pa<-igraph::get.edgelist(g)
      colnames(pa)<-c("GSP1","GSP2") #B large than A, is target
      return(pa)
    }
    sta<-function(plist){
      start_end<-function(plist){
        len<-lapply(plist,function(x){x[[3]]})
        x<-strsplit(as.character(unlist(len)),"=")
        y<-strsplit(unlist(lapply(x,function(x){x[2]})),",")
        z<-lapply(y,function(x){as.numeric(x[1])})
        wid<-unlist(lapply(y,function(x){as.numeric(x[2])}))
        names(z)<-names(len)
        for (i in 1:length(z)){
          z[[i]][2]<-z[[i]]+wid[i]-1
        }
        return(z)
      }
      sta<-unlist(lapply(start_end(plist), function(x){x[1]}))
      return(sta)
    } #start vector
    en<-function(plist){
      start_end<-function(plist){
        len<-lapply(plist,function(x){x[[3]]})
        x<-strsplit(as.character(unlist(len)),"=")
        y<-strsplit(unlist(lapply(x,function(x){x[2]})),",")
        z<-lapply(y,function(x){as.numeric(x[1])})
        wid<-unlist(lapply(y,function(x){as.numeric(x[2])}))
        names(z)<-names(len)
        for (i in 1:length(z)){
          z[[i]][2]<-z[[i]]+wid[i]-1
        }
        return(z)
      }
      en<-unlist(lapply(start_end(plist),function(x){x[2]}))
      return(en)
    } #end vector
    if(class(plist)!="list"){final=NA}else{
      ring<-IRanges::IRanges(sta(plist),en(plist))
      names(ring)<-names(plist) #range
      pri_pair<-as.data.frame(pair(plist,sta=sta(plist))) #pri-pair
      w=di=vector()
      if(nrow(pri_pair)==0){final=NA}else{
        for(i in 1:nrow(pri_pair)){
          di[i]<-if(sum(BSgenome::width(Biostrings::intersect(ring[pri_pair$GSP1[i]],ring[pri_pair$GSP2[i]])))==0) {"diff"} else{"overlaps"}
        }
        for (i in 1: nrow(pri_pair)){
          w[i]<-if(sum(BSgenome::width(Biostrings::intersect(ring[pri_pair$GSP1[i]],ring[pri_pair$GSP2[i]])))==0)
          {a<-sort(c(BSgenome::end(ring[pri_pair$GSP1[i]]),
                     BSgenome::end(ring[pri_pair$GSP2[i]]),
                     BSgenome::start(ring[pri_pair$GSP2[i]]),
                     BSgenome::start(ring[pri_pair$GSP1[i]])))
          a[3]-a[2]

          }
          else {BSgenome::width(Biostrings::intersect(ring[pri_pair$GSP1[i]],
                                ring[pri_pair$GSP2[i]]))}
        }
        evlua<-data.frame(width=w,type=di)
        final<-cbind(pri_pair,evlua)}}
    return(final)
  }
  hot.pair<-list()
  null<-function(x){if(class(x)=="NULL"){y=NA}else{y=x}
    return(y)}
  if(length(plist)==0){hot.pair=NA}else{
    for (i in 1:length(plist)){
      hot.pair[[i]]<- null(ListTodata.frame(plist[[i]]))}
    names(hot.pair)<-names(plist)}
  exchange<-function(x){y=x
  y[,1]<-x[,2]
  y[,2]<-x[,1]
  return(y)}
  if(dir=="left"){pair=hot.pair}else{
    pair<-lapply(hot.pair, exchange)}
  return(pair)
}
