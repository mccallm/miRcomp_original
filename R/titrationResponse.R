## titration response
## look for a monotone relationship across the mixtures
titrationResponse <- function(object1, qcThreshold1, 
                              object2=NULL, qcThreshold2=NULL, 
                              commonFeatures=TRUE,
                              label1=NULL, label2=NULL){

    object1 <- checkObject(object1)
    xs <- as.numeric(gsub("KW","",gsub(":.+","",colnames(object1$ct))))
    ## mixture samples
    pA <- c(1,1,1,1,0,0.2,0.4,0.8)[xs[xs<9]]
    pB <- c(0,0.2,0.4,0.8,1,1,1,1)[xs[xs<9]]

    cts1 <- object1$ct[,xs<9]
    qc1 <- object1$qc[,xs<9]
    
    if(is.null(label1)) label1 <- "Method 1"
    if(is.null(label2)) label2 <- "Method 2"
    
    if(is.null(object2)){
        ## filter miRNAs with any NAs or poor quality values
        i.rm <- which(apply(qc1<qcThreshold1,1,any) | apply(is.na(cts1),1,any))
        cts1 <- cts1[-i.rm,]
    
        ## compute within replicate means for the mixed samples
        iMixB <- which(pA==1 & pB>0)
        iMixA <- which(pB==1 & pA>0)
        repsA <- as.numeric(xs[iMixA])
        repsB <- as.numeric(xs[iMixB])
        ctsMixB <- t(apply(cts1[,iMixB],1,by,repsB,mean))
        ctsMixA <- t(apply(cts1[,iMixA],1,by,repsA,mean))

        ## monotone in A
        monoMixA <- apply(ctsMixA,1,function(x) (x[1]>x[2])&(x[2]>x[3]))
        ## monotone in B
        monoMixB <- apply(ctsMixB,1,function(x) (x[1]>x[2])&(x[2]>x[3]))
        
        ## output table
        tab <- matrix(nrow=2,ncol=2)
        rownames(tab) <- c("Mono","Non-Mono")
        colnames(tab) <- c("A","B")
        tab[1,1] <- sum(monoMixA)
        tab[1,2] <- sum(monoMixB)
        tab[2,1] <- sum(!monoMixA)
        tab[2,2] <- sum(!monoMixB)

        ## figure
        
        ## look at monotonicity by expression
        iPureA <- which(pA==1 & pB==0)
        iPureB <- which(pB==1 & pA==0)

        ctspureA <- rowMeans(cts1[,iPureA])
        ctspureB <- rowMeans(cts1[,iPureB])

        ## strat by difference in expression in pure samples
        d <- c(ctspureA-ctspureB, ctspureB-ctspureA)
        monoqc <- c(monoMixB,monoMixA)

        ix <- sort(d,index.return=T)$ix
        d <- d[ix]
        monoqc <- monoqc[ix]

        dbins <- cut(d, breaks=quantile(d,probs=seq(0,1,0.05)))
        monoqcBins <- unlist(lapply(split(monoqc,dbins),mean))
        xx <- quantile(d,probs=seq(0,1,0.05))
        xx <- (xx[2:length(xx)]+xx[1:(length(xx)-1)])/2
        plot(x=xx,y=monoqcBins,pch=19,ylim=c(0,1),
             xlab="Difference in Pure Sample Expression",
             ylab="Proportion Monotone Increasing")
            
        fit <- glm(monoqc~d,family=binomial(link=logit))
        xx <- seq(min(xx),max(xx),by=0.01)
        lines(x=xx,y=predict(fit,newdata=data.frame(d=xx),type="response"),
              lwd=2, lty=3)
    } else{
        object2 <- checkObject(object2)
        cts2 <- object2$ct[,xs<9]
        qc2 <- object2$qc[,xs<9]
    
        ## filter miRNAs with any NAs or poor quality values
        i.rm1 <- which(apply(qc1<qcThreshold1,1,any) | apply(is.na(cts1),1,any))
        i.rm2 <- which(apply(qc2<qcThreshold2,1,any) | apply(is.na(cts2),1,any))
        if(commonFeatures){
          i.rm <- union(i.rm1,i.rm2)
          cts1 <- cts1[-i.rm,]
          cts2 <- cts2[-i.rm,]
        } else{
          cts1 <- cts1[-i.rm1,]
          cts2 <- cts2[-i.rm2,]
        }
    
        ## compute within replicate means for the mixed samples
        iMixB <- which(pA==1 & pB>0)
        iMixA <- which(pB==1 & pA>0)
        repsA <- as.numeric(xs[iMixA])
        repsB <- as.numeric(xs[iMixB])
        cts1MixB <- t(apply(cts1[,iMixB],1,by,repsB,mean))
        cts1MixA <- t(apply(cts1[,iMixA],1,by,repsA,mean))
        cts2MixB <- t(apply(cts2[,iMixB],1,by,repsB,mean))
        cts2MixA <- t(apply(cts2[,iMixA],1,by,repsA,mean))

        ## monotone in A
        monoMixA1 <- apply(cts1MixA,1,function(x) (x[1]>x[2])&(x[2]>x[3]))
        monoMixA2 <- apply(cts2MixA,1,function(x) (x[1]>x[2])&(x[2]>x[3]))
        ## monotone in B
        monoMixB1 <- apply(cts1MixB,1,function(x) (x[1]>x[2])&(x[2]>x[3]))
        monoMixB2 <- apply(cts2MixB,1,function(x) (x[1]>x[2])&(x[2]>x[3]))

        ## output table
        tab <- matrix(nrow=2,ncol=4)
        rownames(tab) <- c("Mono","Non-Mono")
        colnames(tab) <- paste(rep(c(label1,label2),2),
                               c("A","A","B","B"), sep=":")
        tab[1,1] <- sum(monoMixA1)
        tab[1,2] <- sum(monoMixA2)
        tab[1,3] <- sum(monoMixB1)
        tab[1,4] <- sum(monoMixB2)
        tab[2,1] <- sum(!monoMixA1)
        tab[2,2] <- sum(!monoMixA2)
        tab[2,3] <- sum(!monoMixB1)
        tab[2,4] <- sum(!monoMixB2)

        ## figure
        
        ## look at monotonicity by expression
        iPureA <- which(pA==1 & pB==0)
        iPureB <- which(pB==1 & pA==0)

        cts1pureA <- rowMeans(cts1[,iPureA])
        cts1pureB <- rowMeans(cts1[,iPureB])
        cts2pureA <- rowMeans(cts2[,iPureA])
        cts2pureB <- rowMeans(cts2[,iPureB])

        ## strat by difference in expression in pure samples
        d1 <- c(cts1pureA-cts1pureB, cts1pureB-cts1pureA)
        monoqc1 <- c(monoMixB1,monoMixA1)
        d2 <- c(cts2pureA-cts2pureB, cts2pureB-cts2pureA)
        monoqc2 <- c(monoMixB2,monoMixA2)

        ix1 <- sort(d1,index.return=T)$ix
        d1 <- d1[ix1]
        monoqc1 <- monoqc1[ix1]
        ix2 <- sort(d2,index.return=T)$ix
        d2 <- d2[ix2]
        monoqc2 <- monoqc2[ix2]

        dbins1 <- cut(d1, breaks=quantile(c(d1,d2),probs=seq(0,1,0.05)))
        dbins2 <- cut(d2, breaks=quantile(c(d1,d2),probs=seq(0,1,0.05)))

        monoqcBins1 <- unlist(lapply(split(monoqc1,dbins1),mean))
        monoqcBins2 <- unlist(lapply(split(monoqc2,dbins2),mean))

        xx <- quantile(c(d1,d2),probs=seq(0,1,0.05))
        xx <- (xx[2:length(xx)]+xx[1:(length(xx)-1)])/2
        plot(x=xx,y=monoqcBins1,pch=19,ylim=c(0,1),
             xlab="Difference in Pure Sample Expression",
             ylab="Proportion Monotone Increasing")
        points(x=xx,y=monoqcBins2,pch=22)
            
        fit1 <- glm(monoqc1~d1,family=binomial(link=logit))
        fit2 <- glm(monoqc2~d2,family=binomial(link=logit))
        xx <- seq(min(xx),max(xx),by=0.01)
        lines(x=xx,y=predict(fit1,newdata=data.frame(d1=xx),
                     type="response"), lwd=2, lty=1)
        lines(x=xx,y=predict(fit2,newdata=data.frame(d2=xx),
                     type="response"), lwd=2, lty=3)
        legend("topleft", c(label1,label2), pch=c(19,22), lty=c(1,3))
    }        
    return(tab)
}
