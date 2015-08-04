\name{qualityAssessment}
\alias{qualityAssessment}
\title{Quality Assessment for MicroRNA Expression Estimates}
\description{
  This function produces several figures to examine quality scores.
}
\usage{
qualityAssessment(object1, object2=NULL, cloglog1=FALSE, cloglog2=FALSE, 
                  na.rm=FALSE, plotType=c("scatterplot","boxplot"), 
                  label1=NULL, label2=NULL) 
}
\arguments{
  \item{object1}{a list containing two elements: ct (the expression
    estiamtes) and qc (quality scores)}
  \item{object2}{an optional second list of the same format as
    object1, used to compare two methods.}
  \item{cloglog1}{if TRUE, the -log(-log(object1$qc)) is plotted. This is useful
    to visualize certain quality scores such as R-squared.}
  \item{cloglog2}{if TRUE, the -log(-log(object2$qc)) is plotted. This is useful
      to visualize certain quality scores such as R-squared.}
  \item{na.rm}{if TRUE, quality scores corresponding to NA expression estimates are       removed from analysis. This is only used if plotType=="boxplot".}
  \item{plotType}{type of plot to return. Options are a scatterplot or a
    boxplot.}
  \item{label1}{optional label corresponding to object 1 to be used in plotting.}
  \item{label2}{optional label corresponding to object 2 to be used in plotting.}
}
\value{
  The function generates a plot of the desired type.
}
\author{Matthew N. McCall}
\examples{
  data(lifetech)
  qualityAssessment(lifetech, plotType="boxplot")

  data(qpcRdefault)
  qualityAssessment(object1=lifetech,object2=qpcRdefault,cloglog2=TRUE)
}
\keyword{manip}

