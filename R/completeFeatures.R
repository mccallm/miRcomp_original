completeFeatures <- function(object1, qcThreshold1,
                             object2=NULL, qcThreshold2=NULL,
                             label1=NULL, label2=NULL){

    object1 <- checkObject(object1)
    
    object1$ct[object1$qc<qcThreshold1] <- NA
    tst1 <- rowSums(is.na(object1$ct))

    ## xs <- as.numeric(gsub("KW","",gsub(":.+","",colnames(object1$ct))))
    ## sna <- t(apply(is.na(object1$ct),1,by,xs,sum))

    if(is.null(object2)){
        tab <- matrix(nrow=3, ncol=1)
        tab[1,1] <- sum(tst1==0)
        tab[2,1] <- sum(tst1>0 & tst1<ncol(object1$ct))
        tab[3,1] <- sum(tst1==ncol(object1$ct))
        rownames(tab) <- c("Complete miRNAs (all good quality & non-NA)",
                           "Partial miRNAs (some good quality & non-NA)",
                           "Absent miRNAs (no good quality & non-NA)")
        colnames(tab) <- label1
    } else{
        object2 <- checkObject(object2)
        object2$ct[object2$qc<qcThreshold2] <- NA
        tst2 <- rowSums(is.na(object2$ct))
        tab <- matrix(nrow=3, ncol=3)
        tab[1,1] <- sum(tst1==0 & tst2==0)
        tab[2,1] <- sum((tst1>0 & tst1<ncol(object1$ct)) & tst2==0)
        tab[3,1] <- sum(tst1==ncol(object1$ct) & tst2==0)
        tab[1,2] <- sum(tst1==0 & (tst2>0 & tst2<ncol(object2$ct)))
        tab[2,2] <- sum((tst1>0 & tst1<ncol(object1$ct)) &
                            (tst2>0 & tst2<ncol(object2$ct)))
        tab[3,2] <- sum(tst1==ncol(object1$ct) & 
                          (tst2>0 & tst2<ncol(object2$ct)))
        tab[1,3] <- sum(tst1==0 & tst2==ncol(object2$ct))
        tab[2,3] <- sum((tst1>0 & tst1<ncol(object1$ct)) & 
                          tst2==ncol(object2$ct))
        tab[3,3] <- sum(tst1==ncol(object1$ct) & tst2==ncol(object2$ct))
        
        if(is.null(label1)) label1 <- "Method 1"
        if(is.null(label2)) label2 <- "Method 2"
        rownames(tab) <- paste(label1, c("Complete", "Partial", "Absent"),
                               sep=":")
        colnames(tab) <- paste(label2, c("Complete", "Partial", "Absent"),
                               sep=":")
    }
    return(tab)
}
