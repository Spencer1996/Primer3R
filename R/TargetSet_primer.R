#' TargetSet_primer
#'
#' @param x List set of Primers under each target sequence converted by 'List_PrimerSet' function
#'
#'
#' @return dataframe
#' @export
#'
#' @examples
#' plist<-List_PrimerSet(hot.out)
#' TargetSet_primer(plist)
#'
#' exome
#' sp.ex<-Split_exome(exo.out,frame=exo)
#' eplist<-List_PrimerExomeSet(sp.ex)
#' lapply(eplist,TargetSet_primer) %>% list.rbind()
TargetSet_primer<-function(x){
  Target_primer<-function(x){
    primers<-function(x){ #x<-hp$leukemia_16$Primer1 #factor
      start_width<-function(pp,sw){ #sw,s=1,w=2
        unlist(strsplit(unlist(strsplit(as.character(pp),"="))[2],","))[sw]}
      x<-as.vector(x)
      direction<-unlist(strsplit(x[1],"_"))[2]
      primer.seq<-unlist(strsplit(x[2],"="))[2]
      start<-(as.numeric(start_width(x[3],1)))+2-(as.numeric(start_width(x[3],2)))
      end<-(as.numeric(start_width(x[3],1))+1)
      width<-as.numeric(start_width(x[3],2))
      tm<-as.numeric(unlist(strsplit(x[4],"="))[2])
      GC<-unlist(strsplit(x[5],"="))[2]
      p.frame<-data.frame(Direction=direction,Primer.seq=primer.seq,Start=start,End=end,Width=width,
                          Temp=tm,GC.precent=GC)
      return(p.frame)
    }
    y<- rlist::list.rbind(lapply(x, primers))
    return(y)
  }
y<- rlist::list.rbind(lapply(x,Target_primer))
return(y)}
