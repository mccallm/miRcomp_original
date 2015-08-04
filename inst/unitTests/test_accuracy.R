test_accuracy <- function(){
    data(lifetech)
    tst <- accuracy(object1=lifetech, qcThreshold1=1.25)
    checkEquals(as.vector(tst[2:3,]), c("0.85", "0.28", "0.9",
                                        "0.18", "0.9", "0.14"))
}
