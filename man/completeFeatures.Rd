\name{completeFeatures}
\alias{completeFeatures}
\title{Determine the Number of Completely Observed Features}
\description{
  This function determines the number of features that are good quality and non-NA across all samples using a given quality threshold.
}
\usage{
completeFeatures(object1, qcThreshold1, object2=NULL, qcThreshold2=NULL,
label1=NULL, label2=NULL)
}
\arguments{
  \item{object1}{a list containing two elements: ct (the expression
    estiamtes) and qc (quality scores)}
  \item{qcThreshold1}{a numeric threshold corresponding to object1$qc below which values are considered low quality.}
  \item{object2}{an optional second list of the same format as
    object1, used to compare two methods.}
  \item{qcThreshold2}{a numeric threshold corresponding to object2$qc
    below which values are considered low quality.}
  \item{label1}{optional label corresponding to object 1 to be used in plotting.}
  \item{label2}{optional label corresponding to object 2 to be used in plotting.}
}
\value{
  The function generates a table of the number of complete, partial, and absent features. 
}
\author{Matthew N. McCall}
\examples{
  data(lifetech)
  completeFeatures(object1=lifetech,qcThreshold1=1.25)
  data(qpcRdefault)
  completeFeatures(object1=lifetech,qcThreshold1=1.25,
           object2=qpcRdefault,qcThreshold2=0.99)
}
\keyword{manip}

