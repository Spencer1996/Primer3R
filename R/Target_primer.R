#' Target_primer
#'
#'This fucntion is used to get information of each primer under the target sequecne
#'
#' @param x List of Primers under each target sequence converted by 'List_PrimerSet' function
#'
#' @return data.frame
#' @export
#'
#' @examples
#' plist<-List_PrimerSet(hot.out)
#'Target_primer(plist[[1]])
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
  y<- rlist::list.rbind(lapply(x, primers) )
  return(y)
}
