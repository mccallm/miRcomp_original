test_completeFeatures <- function(){
    data(lifetech)
    tst <- completeFeatures(object1=lifetech, qcThreshold1=1.25)
    checkEquals(as.vector(tst), c(165,375,214))
}
