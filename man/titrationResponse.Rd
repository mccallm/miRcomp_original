\name{titrationResponse}
\alias{titrationResponse}
\title{Assess monotone signal across titration}
\description{
  This function determines 
}
\usage{
titrationResponse(object1, qcThreshold1, object2=NULL, qcThreshold2=NULL, 
                  commonFeatures=TRUE, label1=NULL, label2=NULL)
}
\arguments{
  \item{object1}{a list containing two elements: ct (the expression
    estiamtes) and qc (quality scores)}
  \item{qcThreshold1}{a numeric threshold corresponding to object1$qc below which values are considered low quality.}
  \item{object2}{an optional second list of the same format as
    object1, used to compare two methods.}
  \item{qcThreshold2}{a numeric threshold corresponding to object2$qc
    below which values are considered low quality.}
  \item{commonFeatures}{if TRUE and object2 is non-NULL, only high quality non-NA       features in common between both objects are used.}  
  \item{label1}{optional label corresponding to object 1 to be used in plotting.}
  \item{label2}{optional label corresponding to object 2 to be used in plotting.}
}
\value{
A table listing the number of features showing a monotone
titration response and a figure showing the proportion
of features showing a monotone titration response vs the difference in
pure sample expression.
}
\author{Matthew N. McCall}
\examples{
  data(lifetech)
  titrationResponse(object1=lifetech,qcThreshold1=1.25)
  data(qpcRdefault)
  titrationResponse(object1=lifetech,qcThreshold1=1.25,
           object2=qpcRdefault,qcThreshold2=0.99) 
}
\keyword{manip}

