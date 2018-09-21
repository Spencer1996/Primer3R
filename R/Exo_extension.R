#' Exo_extension
#'
#'For the target sequecne dataframe, if you  need to shear, you can directly use Exo_extension function,
#'which not only shears the target sequence,but also amplifies the fragments for designing the primer. Finally,
#'sequence of sheared fragment and amplified fragment will be returned to the original dataframe, respectively
#'called Fragment.seq and Extend.fragment. Meanwhile, it additionally show the location and length information of sheared
#'fragments. Of course，each of complete amplifyed fragment sequence will also showed at dataframe, called Extend.seq.
#' @param EXOME is dataframe converted by Convert_seq funciton and need to be sheared
#' @param n total lenght of target sequence after amplifying
#' @param cut.length Clip length of target sequecne, the sheared fragment is reuqired to be less than "cut.length"
#' #e.g. for 345bp length of target sequence, it can be sheared into 1-178 and 179-354 two fragments
#' @param chromosome is the chromosome information for amplified fragment
#' e.g chromosome <- BSgenome.Hsapiens.UCSC.hg19
#'
#' @return
#' @export
#'
#' @examples
#' Chr<-c("chrX","chrX" ,"chrX","chrX" )
#' Target.Start<-c(153996577,153999034,154004462,154002877)
#' Target.End<-c(153996707 ,153999154,154004599,154002980)
#' Strand<-rep("+",4)
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19：：BSgenome.Hsapiens.UCSC.hg19
#' frame<-Convert_seq(Chr,Target.Start,Target .End,Strand,chromosome)
#' class<-Class_frame(frame,125)
#' sapply(class,nrow)
#' chromosome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' entend.exome<-Exo_extension(class$EXOME,280,200,chromosome)
#' View(entend.exome)
Exo_extension<-function(EXOME,n,cut.length,chromosome){ #x=frame.eh$EXOME, n is length of extension;length is cut length
  #range for real location, x=frame.eh$EXOME[i,]
  Exome.split<-function(EXOME,length){
    Irange<-function(e,length){
      spp.minus<-function(spp){
        if(nrow(spp)==1){spp=spp}else{
          for (i in 2: nrow(spp)){
            spp$Bin_starts[i]=(spp$Bin_starts[i]+1)
          }}
        return(spp)
      }
      sp.range<-function(x,length){ #x=frame.eh$EXOME[1,]
        ex_range<-function(x,n){
          y<-vector()
          for(i in 1:length(x)){
            y[i]<-as.numeric(unlist(strsplit(as.character(x[i]),"-")))[n]}
          return(y)
        }
        ex_width<-function(x){
          y=s=b=vector()
          for(i in 1:length(x)){
            s[i]<-as.numeric(unlist(strsplit(as.character(x[i]),"-")))[1]
            b[i]<-as.numeric(unlist(strsplit(as.character(x[i]),"-")))[2]
            y[i]<-1+b[i]-s[i]
          }
          return(y)
        }
        e<-GenomicRanges::GRanges(seqnames=x$seqnames,
                   ranges = IRanges::IRanges(start = ex_range(x$range,1),
                                    end = ex_range(x$range,2)),
                   strand = x$strand)
        a<-BSgenome::start(e)
        b<-BSgenome::end(e)
        L=ex_width(x$range)
        N = ceiling(L/length)
        Coordinates = round(seq(from=a,to=b,length.out=N+1))
        sd<-data.frame(Bin_starts = Coordinates[1:N],
                       Bin_ends = Coordinates[2:(N+1)])
        return(sd)
      }
      spp<-spp.minus(sp.range(e,length))
      c<-IRanges::IRanges(spp$Bin_starts,spp$Bin_ends)
      GenomicRanges::mcols(c)$names<-paste0("Exome_c",1:nrow(spp))
      return(c)}
    #
    x<-IRanges::IRangesList()
    for (i in 1:nrow(EXOME)){ #
      x[[i]]<-Irange(EXOME[i,],length)}
    return(x)
  }   #EXOME=frame.eh$EXOME;length=100
  exome.cut<-Exome.split(EXOME,cut.length)
  rowc<-sapply(exome.cut,length)
  gr.exome <-GenomicRanges::GRanges(seqnames= rep(EXOME$seqnames,rowc),
                      ranges = IRanges::IRanges(start = BSgenome::start(unlist(exome.cut)),
                                       end =BSgenome::end(unlist(exome.cut))),
                      strand = rep(EXOME$strand,rowc))
  fix_resize<-function(gr.hotspot,n){
    TF<-vector()
    for(i in 1:length(gr.hotspot)){
      TF[i]<-ifelse(as.character(gr.hotspot@strand)[i]=="+","start","end")}
    y<-GenomicRanges::resize(gr.hotspot, n,fix=TF)
    return(y)
  } #shift the range
  ex.gr.ex<-fix_resize(gr.exome,n)
  extension_seq <- BSgenome::getSeq(chromosome, ex.gr.ex)
  cut_sequence<-BSgenome::getSeq(chromosome, gr.exome)
  ex.n<-split(1:length(extension_seq), rep(1:length(rowc), rowc))

  extend<- GenomicRanges::GRanges(seqnames= rep(EXOME$seqnames,rowc),
                   ranges =IRanges:: IRanges(start = BSgenome::end(gr.exome),
                                    end = BSgenome::end(ex.gr.ex)), strand = rep(as.character(EXOME$strand),rowc))
  Extend.Seq<-BSgenome::getSeq(chromosome, extend)

  extension=exten_start=exten_end=cut_width=cut_seq=ExSeq=list()
  for(i in 1: nrow(EXOME)){
    extension[[i]]<-extension_seq[ex.n[[i]]]
    exten_start[[i]]<-BSgenome::start(ex.gr.ex[ex.n[[i]]])
    exten_end[[i]]<-BSgenome::end(ex.gr.ex[ex.n[[i]]])
    cut_width[[i]]<-BSgenome::width(gr.exome[ex.n[[i]]])
    cut_seq[[i]]<-cut_sequence[ex.n[[i]]]
    ExSeq[[i]]<-Extend.Seq[ex.n[[i]]]
  }
  EXOME$Extend.seq<-lapply(extension,as.character)
  EXOME$Extend.start<- exten_start
  EXOME$Extend.end<-exten_end
  EXOME$Fragment.width<-cut_width
  EXOME$Fragment.seq<-lapply(cut_seq,as.character)
  EXOME$Extend.Fragment<-lapply(ExSeq,as.character)
  return(EXOME)}
