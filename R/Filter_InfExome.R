#' Filter_InfExome
#'
#'Filter_InfExome is specifically applied in process of Exome sequence.
#'it not only filter the primer pairs, but also extract the basical information for each primer
#'
#' @param exomepair.set is the list formation of  Primers under each Exome sheared fragment
#' @param eplist list that integrates sheared fragements and corresponding primers
#' @param distance threshold for distance length between primer pairs, the distacne between primer pairs reqiure
#' to be less than the threshold
#' @param overlap threshold for ovelaps length between primer pairs,the overlap between primer pairs reqiure
#' to be less than the threshold
#'
#' @return list
#' @export
#'
#' @examples
#' sp.ex<-Split_exome(exo.out,frame=exo)
#' eplist<-List_PrimerExomeSet(sp.ex)
#' exomepair.set<-List_PairExomeSet(eplist)
#' Filter_InfExome(exomepair.set,eplist,distance=20,overlap=10)
Filter_InfExome<-function(exomepair.set,eplist,distance,overlap,inf=TRUE){
  Filter_Inf<-function(hot.pair,plist,distance,overlap,inf=TRUE){
    Filterset<-function(final,difference,overlap){
      Filter<-function(a,difference,overlap){
        if(a$type=="diff") {ifelse(a$width<difference,a$width,NA)
        } else {ifelse(a$width<overlap,a$width,NA)}
      }
      if(class(final)!="data.frame"){final.1=NA}else{
        for (i in 1:nrow(final)){
          final[i,3]<-Filter(final[i,],difference,overlap)
        }
        final.1=na.omit(final)}
      return(final.1)
    } #filer the data
    filer.hot.pair<-list()
    if(length(hot.pair)==0){result=NA} else{
      for (i in 1:length(hot.pair)){
        filer.hot.pair[[i]]<-Filterset(hot.pair[[i]],distance,overlap)}
      names(filer.hot.pair)<-names(hot.pair)
      result<-filer.hot.pair}
    result[which(sapply(result,nrow)==0)]<-"None of them fit"
    #get information
    if(inf==FALSE){INF=result}else{
      Pairs.inf<-function(fil_hot.pair,plist){
        pair_inf<-function(fil,primer3){
          fil=as.data.frame(fil)
          primer3=as.data.frame(primer3)
          primer_detail<-function(x,y){#x=x<-fil[1,], y=primer3,data.frame
            if(class(x)!="data.frame"){pr<-NA}else{
              x<-x[,1:2]
              frame.list<-list()
              tm=nu=g=wr=seq=start=end=vector()
              for (i in 1:length(x)){
                nu[i]<-as.numeric(unlist(strsplit(as.character(x[[i]]),"Primer"))[2])
                g[i]<-gsub("0",as.character(nu[i]-1),"PRIMER_RIGHT_0_SEQUENCE=")
                wr[i]<-grep(g[i],y[,1])
                seq[i]<-unlist(strsplit(as.character(y[wr[i],1]),"="))[2]
                tm[i]<-unlist(strsplit(as.character(y[(wr[i]+2),1]),"="))[2]
                #pp[i]<-y[(wr[i]+1),]
                start_width<-function(pp,sw){ #sw,s=1,w=2
                  unlist(strsplit(unlist(strsplit(as.character(pp),"="))[2],","))[sw]} #sw,s=1,w=2
                start[i]=start_width(y[(wr[i]+1),1],1)
                end[i]=start_width(y[(wr[i]+1),1],2)
                frame.list[[i]]<-data.frame(Seq=seq[i],Start=(as.numeric(start[i])+1),
                                            End=(as.numeric(start[i])+2-as.numeric(end[i])),Length=end[i],TM=tm[i])}

              names(frame.list)<-paste0("",colnames(x))
              pr<-cbind(x,rlist::list.cbind(frame.list))}
            return(pr)
          } #x=x<-fil[1,], y=primer3,data.frame
          primer_LIST<-list()
          if(nrow(fil)==1){y="None of them fit"}else{
            for(i in 1:nrow(fil)){
              primer_LIST[[i]]<-primer_detail(fil[i,],primer3)

            }
            y<-rlist::list.rbind(primer_LIST)}
          return(y)}
        unlist.primer<-function(x){
          names(x)=NULL
          y=rlist::list.rbind(lapply(x,as.data.frame))
          return(y)}
        PRIMER.list<-lapply(plist,unlist.primer)
        hot.inf<-list()
        for (i in 1:length(fil_hot.pair)){
          hot.inf[[i]]<-pair_inf(fil_hot.pair[[i]],PRIMER.list[names(fil_hot.pair[i])])
        }
        names(hot.inf)<-names(fil_hot.pair)
        return(hot.inf)}
      INF<-Pairs.inf(result,plist)
    }
    return(INF)}
Result=list()
for (i in 1:length(exomepair.set)) {#length(eplist)
  Result[[i]]= Filter_Inf(hot.pair=exomepair.set[[i]],plist=eplist[[i]],distance=distance,overlap=overlap,inf=inf)
}
names(Result)<-names(eplist)
return(Result)
}
