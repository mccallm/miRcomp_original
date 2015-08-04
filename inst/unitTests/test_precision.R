test_precision <- function(){
    data(lifetech)
    tst <- precision(object1=lifetech, qcThreshold1=1.25)
    checkEqualsNumeric(as.vector(tst[[1]][1:3]), 
                       c(0.1683349, 0.1147788, 0.1715636), 
                       tolerance=1.0e-4)
}
