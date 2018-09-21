#' Split_exome
#'
#' @param exo.out is the output exome result of Primer3
#' @param frame is the exome targe sequence dataframe
#'
#' @return list
#' @export
#'
#' @examples
#' library(xxx)
#' data(exo.out)
#' data(exo)
#' sp.ex<-Split_exome(exo.out,frame=exo)
Split_exome<-function(exo.out,frame){
  n<-nrow(frame)
  g<-list()
  list.ex0<-list()
  for (i in 1:n){
    g[[i]]=grep(paste0("SEQUENCE_ID=",frame$NO.[i],"."),exo.out[,1])}
  for (i in 1:n) {
    list.ex0[[i]]<-exo.out[g[[i]][1]:ifelse(i==n,nrow(exo.out),(g[[i+1]][1]-1)),]
  }

  names(list.ex0)<-frame$NO.
  return(list.ex0)
}
