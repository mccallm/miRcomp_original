\name{limitOfDetection}
\alias{limitOfDetection}
\title{Assess the limit of detection of a given method}
\description{
  This function assesses the limit of detection in one of two ways: (1)
  the distribution of expression estimates stratified by the proportion
  of poor quality values within replicates, (2) the average vs expected
  expression for the two most diluted sample types.
}
\usage{
limitOfDetection(object, qcThreshold, plotType=c("boxplot","scatterplot","MAplot"))
}
\arguments{
  \item{object}{a list containing two elements: ct (the expression
    estiamtes) and qc (quality scores)}
  \item{qcThreshold}{a numeric threshold corresponding to object1$qc
    below which values are considered low quality.} 
  \item{plotType}{the desired output type -- boxplot is option (1); 
    scatterplot is option (2), an MA-plot is option (3).}
}
\value{
  This function assesses the limit of detection in several ways. If
  plotType is boxplot, then boxplots of expression estimates stratified
  by the proportion of poor quality values within replicates is
  displayed. The function also outputs a list with the values plotted in
  each box of the boxplot. If plotType is scatterplot, then the average within
  replicate expression vs expected expression (based on pure sample
  expression) is displayed for the 0.1/0.1 dilution and 0.01/0.01
  dilution. If plotType is MAplot, then the difference in expression (average within
  replicate expression - expected expression) is displayed for the
  0.1/0.1 dilution and 0.01/0.01 dilution. For both plotTypes,
  scatterplot and MAplot, the function outputs a matrix containing
  estimates of the limit of detection for four different
  tolerances. Specifically, the two columns
  correspond to the two dilutions (0.1/0.1 and 0.01/0.01) and rows
  correspond to the median difference between the observed and expected
  values. The values in the matrix are the expected expression values
  such that the median absolute difference of all larger expected
  expression values is approximately equal to the given tolerance. 
}
\author{Matthew N. McCall}
\examples{
  data(lifetech)
  tmp <- limitOfDetection(object=lifetech,qcThreshold=1.25)
  data(qpcRdefault)
  limitOfDetection(object=qpcRdefault,qcThreshold=0.99,plotType="scatter")
}
\keyword{manip}

