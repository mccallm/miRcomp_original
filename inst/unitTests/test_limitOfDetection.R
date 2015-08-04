test_limitOfDetection <- function(){
    data(lifetech)
    tst <- limitOfDetection(object=lifetech, qcThreshold=1.25, 
                            plotType="scatterplot")
    checkEquals(as.vector(tst), c(27.6, 28.9, 29.1, 
                                  26.8, 28.4, 29.0, 
                                  26.3, 28.4, 29.2))
}
