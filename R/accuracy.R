## look at accuracy of the signal detect curves
## using the nominal change in expression
accuracy <- function(object1, qcThreshold1,
                     object2=NULL, qcThreshold2=NULL, 
                     commonFeatures=TRUE,
                     bins=3, label1=NULL, label2=NULL){

    object1 <- checkObject(object1)
    
    xs <- as.numeric(gsub("KW","",gsub(":.+","",colnames(object1$ct))))
    cts1 <- object1$ct[,xs<9]
    qc1 <- object1$qc[,xs<9]
    if(!is.null(object2)){
        object2 <- checkObject(object2)
        cts2 <- object2$ct[,xs<9]
        qc2 <- object2$qc[,xs<9]
    }
    xs <- xs[xs<9]

    ## mixture samples
    pA <- c(1,1,1,1,0,0.2,0.4,0.8)[xs]
    pB <- c(0,0.2,0.4,0.8,1,1,1,1)[xs]

    if(is.null(object2)){
        ## filter miRNAs with any NAs or poor quality values
        i.rm <- which(apply(qc1<qcThreshold1,1,any) | apply(is.na(cts1),1,any))
        cts1 <- cts1[-i.rm,]

        ## expression in pure samples
        ctspureA <- rowMeans(cts1[,pA==1 & pB==0])
        ctspureB <- rowMeans(cts1[,pA==0 & pB==1])
        
        coefmixA <- matrix(nrow=nrow(cts1),ncol=4)
        colnames(coefmixA) <- c("slope","SE","tstat","pval")
        ctsmixA <- cts1[,xs%in%c(6,7,8)]
        pAmixA <- pA[xs%in%c(6,7,8)]
        for(k in 1:nrow(ctsmixA)){
            X <- -log2(((2^-ctspureB[k])+(pAmixA*(2^-ctspureA[k]))))
            if(var(X)<0.01){
              coefmixA[k,] <- rep(NA,4)
            } else{
              coefmixA[k,] <- coef(summary(lm(ctsmixA[k,]~X)))[2,]
            }
        }    
        coefmixB <- matrix(nrow=nrow(cts1),ncol=4)
        colnames(coefmixB) <- c("slope","SE","tstat","pval")
        ctsmixB <- cts1[,xs%in%c(2,3,4)]
        pBmixB <- pB[xs%in%c(2,3,4)]
        for(k in 1:nrow(ctsmixB)){
            X <- -log2(((2^-ctspureA[k])+(pBmixB*(2^-ctspureB[k]))))
            if(var(X)<0.01){
              coefmixB[k,] <- rep(NA,4)
            } else{
              coefmixB[k,] <- coef(summary(lm(ctsmixB[k,]~X)))[2,]
            }
        }    
        coef <- rbind(coefmixA,coefmixB)
        slope <- coef[,"slope"]

        ## stratify by pure sample expression difference
        d <- c(ctspureB-ctspureA, ctspureA-ctspureB)

        ## consider only the titration series where the titrating sample
        ## has higher expression than the sample held constant
        ind <- which(d>0)

        grp <- cut(d[ind],breaks=quantile(d[ind],
                         probs=seq(0,1,length.out=bins+1))
                   -c(0.1,rep(0,bins-1),-0.1))
        ysBin <- split(slope[ind],grp)
        cols <- ifelse(coef[ind,"pval"]<0.05,"black","grey")
        plot(x=jitter(rep(1:bins,lapply(ysBin,length)),factor=0.75),
             y=unlist(ysBin), col=cols, pch=21,
             ylab="Signal Detect Slope",xlab="",xaxt="n")
        if(bins>1){
          if(bins==2){
            axis(side=1,labels=c("Low","High"),at=c(1,bins), tick=FALSE, 
                 line=1)
          } else{
            axis(side=1,labels=c("Low","Medium","High"),
                 at=c(1,(bins+1)/2,bins), tick=FALSE, line=1)
          }
          mtext("Difference in Expression (Titrating - Constant)",
                side=1,line=3)
        }
        legend("topleft",pch=21,col=c("black","grey"),
               c("pval<0.05","pval>0.05"))
        abline(h=1,lwd=2,lty=3)

        tmp <- round(sapply(ysBin,function(x) c(median(x),mad(x))),digits=2)
        out1 <- rbind(colnames(tmp),tmp)
        rownames(out1) <- c("Bin","Median","MAD")
        colnames(out1) <- paste0("bin",1:bins)
        return(out1)
        
    } else{

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

        ## compute pure sample expression
        ctspureA1 <- rowMeans(cts1[,pA==1 & pB==0])
        ctspureB1 <- rowMeans(cts1[,pA==0 & pB==1])
        ctspureA2 <- rowMeans(cts2[,pA==1 & pB==0])
        ctspureB2 <- rowMeans(cts2[,pA==0 & pB==1])
        
        coefmixA1 <- matrix(nrow=nrow(cts1),ncol=4)
        coefmixA2 <- matrix(nrow=nrow(cts2),ncol=4)
        colnames(coefmixA1) <- colnames(coefmixA2) <- c("slope","SE",
                                                        "tstat","pval")
        ctsmixA1 <- cts1[,xs%in%c(6,7,8)]
        ctsmixA2 <- cts2[,xs%in%c(6,7,8)]
        pAmixA <- pA[xs%in%c(6,7,8)]
        for(k in 1:nrow(ctsmixA1)){
            X1 <- -log2(((2^-ctspureB1[k])+(pAmixA*(2^-ctspureA1[k]))))
            if(var(X1)<0.01){
              coefmixA1[k,] <- rep(NA,4)
            } else{
              coefmixA1[k,] <- coef(summary(lm(ctsmixA1[k,]~X1)))[2,]
            }
        }
        for(k in 1:nrow(ctsmixA2)){
            X2 <- -log2(((2^-ctspureB2[k])+(pAmixA*(2^-ctspureA2[k]))))
            if(var(X2)<0.01){
              coefmixA2[k,] <- rep(NA,4)
            } else{
              coefmixA2[k,] <- coef(summary(lm(ctsmixA2[k,]~X2)))[2,]
            }
        }    
        
        coefmixB1 <- matrix(nrow=nrow(cts1),ncol=4)
        coefmixB2 <- matrix(nrow=nrow(cts2),ncol=4)
        colnames(coefmixB1) <- colnames(coefmixB2) <- c("slope","SE",
                                                        "tstat","pval")
        ctsmixB1 <- cts1[,xs%in%c(2,3,4)]
        ctsmixB2 <- cts2[,xs%in%c(2,3,4)]
        pBmixB <- pB[xs%in%c(2,3,4)]
        for(k in 1:nrow(ctsmixB1)){
          X1 <- -log2(((2^-ctspureA1[k])+(pBmixB*(2^-ctspureB1[k]))))
          if(var(X1)<0.01){
            coefmixB1[k,] <- rep(NA,4)
          } else{
            coefmixB1[k,] <- coef(summary(lm(ctsmixB1[k,]~X1)))[2,]
          }
        }       
        for(k in 1:nrow(ctsmixB2)){
          X2 <- -log2(((2^-ctspureA2[k])+(pBmixB*(2^-ctspureB2[k]))))
          if(var(X2)<0.01){
            coefmixB2[k,] <- rep(NA,4)
          } else{
            coefmixB2[k,] <- coef(summary(lm(ctsmixB2[k,]~X2)))[2,]
          }
        }    
        coef1 <- rbind(coefmixA1,coefmixB1)
        coef2 <- rbind(coefmixA2,coefmixB2)
        slope1 <- coef1[,"slope"]
        slope2 <- coef2[,"slope"]

        ## stratify by pure sample expression difference
        d1 <- c(ctspureB1-ctspureA1, ctspureA1-ctspureB1)
        d2 <- c(ctspureB2-ctspureA2, ctspureA2-ctspureB2)

        ## consider only the titration series where the titrating sample
        ## has higher expression than the sample held constant
        ind1 <- which(d1>0)
        ind2 <- which(d2>0)
        
        grp1 <- cut(d1[ind1],breaks=quantile(d1[ind1],
                         probs=seq(0,1,length.out=bins+1))-c(0.1,rep(0,bins-1),-0.1))
        ysBin1 <- split(slope1[ind1],grp1)
        cols1 <- ifelse(coef1[ind1,"pval"]<0.05,"black","grey")

        grp2 <- cut(d2[ind2],breaks=quantile(d2[ind2],
                         probs=seq(0,1,length.out=bins+1))-c(0.1,rep(0,bins-1),-0.1))
        ysBin2 <- split(slope2[ind2],grp2)
        cols2 <- ifelse(coef2[ind2,"pval"]<0.05,"black","grey")
        
        plot(x=c(jitter(rep(seq(1,(bins*2)-1,by=2),lapply(ysBin1,length)),factor=0.5),
               jitter(rep(seq(2,(bins*2),by=2),lapply(ysBin2,length)),factor=0.5)),
             y=c(unlist(ysBin1),unlist(ysBin2)), col=c(cols1,cols2), pch=21,
             ylab="Signal Detect Slope",xlab="",xaxt="n")
        if(is.null(label1)) label1 <- "M1"
        if(is.null(label2)) label2 <- "M2"
        axis(side=1, labels=rep(c(label1,label2),bins), at=1:(2*bins))
        if(bins>1){
          if(bins==2){
            axis(side=1,labels=c("Low","High"),at=c(1.5,2*bins-0.5), 
                 tick=FALSE, line=1)
          } else{
            axis(side=1,labels=c("Low","Medium","High"),
                 at=c(1.5,bins+0.5,2*bins-0.5), tick=FALSE, line=1)
          }
          mtext("Difference in Expression (Titrating - Constant)",
                side=1,line=3)
        }
        legend("topleft",pch=21,col=c("black","grey"),
               c("pval<0.05","pval>0.05"))
        abline(h=1,lwd=2,lty=3)
    
        tmp1 <- round(sapply(ysBin1,function(x) c(median(x),mad(x))),digits=2)
        out1 <- rbind(colnames(tmp1),tmp1)
        rownames(out1) <- c("Bin","Median","MAD")
        colnames(out1) <- paste0("bin",1:bins)
        
        tmp2 <- round(sapply(ysBin2,function(x) c(median(x),mad(x))),digits=2)
        out2 <- rbind(colnames(tmp2),tmp2)
        rownames(out2) <- c("Bin","Median","MAD")
        colnames(out2) <- paste0("bin",1:bins)
        
        out <- list("Method1"=out1, "Method2"=out2)
        if(!is.null(label1)) names(out)[1] <- label1
        if(!is.null(label2)) names(out)[2] <- label2
        return(out)
    }
}
