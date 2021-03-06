% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hot_output.R
\name{Hot_output}
\alias{Hot_output}
\title{Hot_output}
\usage{
Hot_output(fil.hot, reserve = FALSE, Write.out = FALSE, NAME)
}
\arguments{
\item{fil.hot}{is the final filtered result of hotspot sequecne, namely non-sheared sequence.}

\item{reserve}{can reserve the filtered row, default is FALSE .}

\item{Write.out}{it can write the result into '.csv',default is FALSE.}

\item{NAME}{name of output file with '.CSV' formation.}
}
\description{
This function can write the filtered result into '.CSV' formation. Here, it only chooses the optimal
priemr pairs for each target sequence. If there is not primer, it only returns NA result; if  no
primer pairs fit the condiditons, it returns "None of them fit" at the first column
}
\examples{
plist<-List_PrimerSet(hot.out)
hot.pair<-List_PairSet(plist)
rem_Pairs.inf<-Filter_Inf(hot.pair,plist,distance=20,overlap=10)
Hot_output(rem_Pairs.inf)

Output:
Hot_output(rem_Pairs.inf,Write.out=TRUE,NAME="hotResult")
}
