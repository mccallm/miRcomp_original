\name{accuracy}
\alias{accuracy}
\title{Assess the accuracy of the expression estimates}
\description{
  This function assesses the signal detect slope (a measure of accuracy)
  for each feature.
}
\usage{
accuracy(object1, qcThreshold1, object2=NULL, qcThreshold2=NULL, 
commonFeatures=TRUE, bins=3, label1=NULL, label2=NULL)
}
\arguments{
  \item{object1}{a list containing two elements: ct (the expression
    estimates) and qc (quality scores)}
  \item{qcThreshold1}{a numeric threshold corresponding to object1$qc below which values are considered low quality.}
  \item{object2}{an optional second list of the same format as
    object1, used to compare two methods.}
  \item{qcThreshold2}{a numeric threshold corresponding to object2$qc
    below which values are considered low quality.}
  \item{commonFeatures}{if TRUE and object2 is non-NULL, only high quality non-NA       features in common between both objects are used.}  
  \item{bins}{the number of bins to divide the data into.}
  \item{label1}{optional label corresponding to object 1 to be used in plotting.}
  \item{label2}{optional label corresponding to object 2 to be used in plotting.}
}
\value{
  A plot of signal detect slopes stratified by difference in pure sample
  expression is produced, and a summary (median and MAD) of the signal
  detect slopes in each bin is returned.
}
\author{Matthew N. McCall}
\examples{
  data(lifetech)
  accuracy(object1=lifetech,qcThreshold1=1.25)
  data(qpcRdefault)
  accuracy(object1=lifetech,qcThreshold1=1.25,
           object2=qpcRdefault,qcThreshold2=0.99)
}
\keyword{manip}

