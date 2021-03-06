% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Filter_InfExome.R
\name{Filter_InfExome}
\alias{Filter_InfExome}
\title{Filter_InfExome}
\usage{
Filter_InfExome(exomepair.set, eplist, distance, overlap, inf = TRUE)
}
\arguments{
\item{exomepair.set}{is the list formation of  Primers under each Exome sheared fragment}

\item{eplist}{list that integrates sheared fragements and corresponding primers}

\item{distance}{threshold for distance length between primer pairs, the distacne between primer pairs reqiure
to be less than the threshold}

\item{overlap}{threshold for ovelaps length between primer pairs,the overlap between primer pairs reqiure
to be less than the threshold}
}
\value{
list
}
\description{
Filter_InfExome is specifically applied in process of Exome sequence.
it not only filter the primer pairs, but also extract the basical information for each primer
}
\examples{
sp.ex<-Split_exome(exo.out,frame=exo)
eplist<-List_PrimerExomeSet(sp.ex)
exomepair.set<-List_PairExomeSet(eplist)
Filter_InfExome(exomepair.set,eplist,distance=20,overlap=10)
}
