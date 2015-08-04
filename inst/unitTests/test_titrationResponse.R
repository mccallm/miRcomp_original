test_titrationResponse <- function(){
    data(lifetech)
    tst <- titrationResponse(object1=lifetech, qcThreshold1=1.25)
    checkEquals(as.vector(tst), c(98, 109, 164, 43))
}
