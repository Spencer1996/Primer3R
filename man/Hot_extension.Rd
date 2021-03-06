% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hot_extension.R
\name{Hot_extension}
\alias{Hot_extension}
\title{Hot_extension}
\usage{
Hot_extension(HOTSPOT, n, chromosome)
}
\arguments{
\item{HOTSPOT}{is dataframe converted by Convert_seq funciton}

\item{n}{total lenght of target sequence after amplifying}

\item{chromosome}{is the chromosome information for amplified fragment
e.g chromosome <- BSgenome.Hsapiens.UCSC.hg19}
}
\value{
data.fame
}
\description{
For the target sequecne dataframe, if you do not need to shear, you can directly use it to
extend the target sequence for designing the primer. Finally, amplified fragment will return
in the target sequecne dataframe, called Extend.fragment. Meanwhile, it additionally show the complete
sequence after amplifying, called Extend.seq
}
\examples{
Chr<-c("chrX","chrX" ,"chrX","chrX" )
Target.Start<-c(153996577,153999034,154004462,154002877)
Target.End<-c(153996707 ,153999154,154004599,154002980)
Strand<-rep("+",4)
chromosome <- BSgenome.Hsapiens.UCSC.hg19
frame<-Convert_seq(Chr,Target.Start,Target.End,Strand,chromosome)
class<-Class_frame(frame,125)
sapply(class,nrow)
chromosome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
entend.hot<-Hot_extension(class$HOTSPOT,185,chromosome)
View(entend.hot)
}
