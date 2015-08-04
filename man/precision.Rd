\name{precision}
\alias{precision}
\title{Assess the precision of the expression estimates}
\description{
  This function assesses the within-replicate precision for each feature.
}
\usage{
precision(object1, qcThreshold1, object2=NULL, qcThreshold2=NULL,
          commonFeatures=TRUE, statistic=c("sd","cv"),
          scale=c("none","log","log10"), bins=3,
          label1=NULL, label2=NULL) 
}
\arguments{
  \item{object1}{a list containing two elements: ct (the expression
    estiamtes) and qc (quality scores)}
  \item{qcThreshold1}{a numeric threshold corresponding to object1$qc
    below which values are considered low quality.} 
  \item{object2}{an optional second list of the same format as
    object1, used to compare two methods.}
  \item{qcThreshold2}{a numeric threshold corresponding to object2$qc
    below which values are considered low quality.}
  \item{commonFeatures}{if TRUE and object2 is non-NULL, only high quality non-NA       features in common between both objects are used.}  
  \item{statistic}{whether to compute the standard deviation (sd) or
    coefficient of variation (cv).} 
  \item{scale}{optional scaling of the values. This can help with
    visualizing the distributions.}
  \item{bins}{the number of bins to divide the data into.}
  \item{label1}{optional label corresponding to object 1 to be used in plotting.}
  \item{label2}{optional label corresponding to object 2 to be used in plotting.}
}
\value{
  A boxplot of either the standard deviation or coefficient of variation
  stratified by expression is produced. The values plotted in each box
  of the boxplot are returned.
}
\author{Matthew N. McCall}
\examples{
  data(lifetech)
  tmp1 <- precision(object1=lifetech,qcThreshold1=1.25)
  data(qpcRdefault)
  tmp2 <- precision(object1=lifetech,qcThreshold1=1.25,
           object2=qpcRdefault,qcThreshold2=0.99)
}
\keyword{manip}

