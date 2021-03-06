% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Exo_output.R
\name{Exo_output}
\alias{Exo_output}
\title{Exo_output}
\usage{
Exo_output(fil.exo, exo.ex, reserve = FALSE, Write.out = FALSE, NAME)
}
\arguments{
\item{fil.exo}{is the final filtered result of exome sequecne, namely non-sheared sequence}

\item{exo.ex}{is the exome targe sequence dataframe after extension processed by
'Exo_extension' funciton}

\item{reserve}{can reserve the filtered row, default is FALSE .}

\item{Write.out}{it can write the result into '.csv',default is FALSE.}

\item{NAME}{name of output file with '.CSV' formation.}
}
\description{
This function can write the filtered result exome sequecne into '.CSV' formation. Here, it only chooses the optimal
priemr pairs for each target sequence. If there is not primer, it only returns NA result; if  no
primer pairs fit the condiditons, it returns "None of them fit" at the first column
}
\examples{
sp.ex<-Split_exome(exo.out,frame=exo)
exomepair.set<-List_PrimerExomeSet(sp.ex) \%>\% List_PairExomeSet()
fil.exo<-Filter_InfExome(exomepair.set,eplist,distance=20,overlap=10)
entend.exome<-Exo_extension(exo,280,200,chromosome)
Exo_output(fil.exo, entend.exome$Extend.seq)

Output:
Exo_output(fil.exo, entend.exome$Extend.seq,Write.out=TRUE,NAME="exoResult")
}
