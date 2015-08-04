## figures to assess the limit of detection
limitOfDetection <- function(object, qcThreshold,
                             plotType=c("boxplot","scatterplot","MAplot")){
    plotType <- match.arg(plotType)

    object <- checkObject(object)
    object$ct[object$qc<qcThreshold] <- NA
    xs <- as.numeric(gsub("KW","",gsub(":.+","",colnames(object$ct))))

    if(plotType=="boxplot"){
        uconds <- unique(xs)
        mus <- pnd <- matrix(nrow=nrow(object$ct),ncol=length(uconds))
        for(k in 1:length(uconds)){
            mus[,k] <- apply(object$ct[,xs==uconds[k]],1,mean,na.rm=TRUE)
            pnd[,k] <- rowMeans(is.na(object$ct[,xs==uconds[k]]))
        }
        mus <- as.vector(mus)
        pnd <- as.vector(pnd)
        i.rm <- which(is.na(mus))
        mus <- mus[-i.rm]
        pnd <- pnd[-i.rm]
        
        boxes <- list("0.00"=mus[pnd==0],
                     "0.25"=mus[pnd==0.25],
                     "0.5"=mus[pnd==0.5],
                     "0.75"=mus[pnd==0.75],
                     "1.00"=rep(NA,length(i.rm)))
        boxplot(boxes, ylab="Average Observed Ct Value")
        axis(side=1,tick=F,at=1:5,
             labels=paste0("(",sapply(boxes,length),")"),line=1.15)
        mtext("Proportion Poor Quality\n(number of data points)",
              side=1,line=5)

        return(boxes)
    }

    if(plotType%in%c("scatterplot","MAplot")){
        ## use pure samples to compute expected expression in
        ## 0.1/0.1 and 0.01/0.01 mixtures
        pA <- c(1,1,1,1,0,0.2,0.4,0.8,0.1,0.01)[xs]
        pB <- c(0,0.2,0.4,0.8,1,1,1,1,0.1,0.01)[xs]
        
        ctspureA <- object$ct[,pA==1 & pB==0]
        qcpureA <- object$qc[,pA==1 & pB==0]
        ctspureB <- object$ct[,pA==0 & pB==1]
        qcpureB <- object$qc[,pA==0 & pB==1]

        sna <- apply(is.na(cbind(ctspureA,ctspureB)),1,sum)
        ind <- which(sna==0)

        ctspureA <- ctspureA[ind,]
        ctspureB <- ctspureB[ind,]

        ## examine quality and expression in the 0.1/0.1 samples 
        cts1 <- object$ct[ind,pA==0.1 & pB==0.1]
        qc1 <- object$qc[ind,pA==0.1 & pB==0.1]

        ## examine quality and expression in the 0.01/0.01 samples 
        cts01 <- object$ct[ind,pA==0.01 & pB==0.01]
        qc01 <- object$qc[ind,pA==0.01 & pB==0.01]

        ## any NAs
        iNA1 <- apply(is.na(cts1),1,sum)
        iNA01 <- apply(is.na(cts01),1,sum)
        ## all NAs
        iAllNA1 <- which(iNA1==4)
        iAllNA01 <- which(iNA01==4)
        
        ## plot average expression in pure samples
        
        Exp1 <- -log2(((2^-rowMeans(ctspureA))+(2^-rowMeans(ctspureB)))/10)
        ys1 <- rowMeans(cts1,na.rm=T)
        M1 <- ys1-Exp1

        Exp01 <- -log2(((2^-rowMeans(ctspureA))+(2^-rowMeans(ctspureB)))/100)
        ys01 <- rowMeans(cts01,na.rm=T)
        M01 <- ys01-Exp01
        
        ## for the 0.1/0.1 vs 0.01/0.01 comparison
        ## remove ones that are all NAs in the 0.1/0.1 dilution
        Exp101 <- ys1[-iAllNA1]+log2(10)
        M101 <- ys01[-iAllNA1]-Exp101
        
        ## calculate output table
        es1 <- seq(min(Exp1),max(Exp1),by=0.01)
        ps1 <- sapply(es1,function(k) median(abs(M1)[Exp1>k],na.rm=TRUE))
        es01 <- seq(min(Exp01),max(Exp01),by=0.01)
        ps01 <- sapply(es01,function(k) median(abs(M01)[Exp01>k],na.rm=TRUE))
        es101 <- seq(min(Exp101),max(Exp101),by=0.01)
        ps101 <- sapply(es101,function(k) median(abs(M101)[Exp101>k],na.rm=TRUE))
        
        out <- matrix(ncol=3,nrow=3)
        out[,1] <- es1[c(which.min(abs(ps1-0.5)),
                         which.min(abs(ps1-0.75)),
                         which.min(abs(ps1-1)))]
        out[,2] <- es01[c(which.min(abs(ps01-0.5)),
                          which.min(abs(ps01-0.75)),
                          which.min(abs(ps01-1)))]
        out[,3] <- es01[c(which.min(abs(ps101-0.5)),
                          which.min(abs(ps101-0.75)),
                          which.min(abs(ps101-1)))]
        colnames(out) <- c("0.1/0.1 vs pure",
                           "0.01/0.01 vs pure",
                           "0.01/0.01 vs 0.1/0.1")
        rownames(out) <- c("0.50","0.75","1.00")
        
        ## make plots
        if(plotType=="scatterplot"){
          ## 0.1/0.1 vs Expected (pure)
            toptick <- ceiling(max(c(Exp1,ys1),na.rm=T)/5)*5
            lowtick <- floor(min(c(Exp1,ys1),na.rm=T)/5)*5
            if(length(iAllNA1)>0) ys1[iAllNA1] <- toptick+5
            plot(x=Exp1,y=ys1,pch=c(19,17,15,18,4)[iNA1+1],
                 col=ifelse(apply(qc1>qcThreshold,1,all),"black","red"),
                 main="0.1/0.1 Dilution",
                 xlab="Expected Expression\n(based on pure samples)",
                 ylab="Average Observed Expression", yaxt="n",
                 ylim=c(lowtick,toptick+5),xlim=c(lowtick,toptick))
            axis(side=2,at=seq(lowtick,toptick+5,by=5),
                 labels=c(as.character(seq(lowtick,toptick,by=5)),"NA"),las=2)
            abline(a=0,b=1,lwd=2,lty=3,col="blue")
            legend("topleft",pch=c(19,17,15,18,4),title="Proporiton Poor Quality",
                   c("0/4","1/4","2/4","3/4","4/4"),ncol=2,cex=0.75,
                   col=c("black",rep("red",4)))

            ## 0.01/0.01 vs Expected (pure)
            toptick <- ceiling(max(c(Exp01,ys01),na.rm=T)/5)*5
            lowtick <- floor(min(c(Exp01,ys01),na.rm=T)/5)*5
            if(length(iAllNA01)>0) ys01[iAllNA01] <- toptick+5
            plot(x=Exp01,y=ys01,pch=c(19,17,15,18,4)[iNA01+1],
                 col=ifelse(apply(qc01>qcThreshold,1,all),"black","red"),
                 main="0.01/0.01 Dilution",
                 xlab="Expected Expression\n(based on pure samples)",
                 ylab="Average Observed Expression",yaxt="n",
                 ylim=c(lowtick,toptick+5),xlim=c(lowtick,toptick))
            axis(side=2,at=seq(lowtick,toptick+5,by=5),
                 labels=c(as.character(seq(lowtick,toptick,by=5)),"NA"),las=2)
            abline(a=0,b=1,lwd=2,lty=3,col="blue")
            legend("topleft",pch=c(19,17,15,18,4),title="Proporiton Poor Quality",
                   c("0/4","1/4","2/4","3/4","4/4"),ncol=2,cex=0.75,
                   col=c("black",rep("red",4)))
            
            ## 0.01/0.01 vs Expected (0.1/0.1)
            toptick <- ceiling(max(c(Exp101,ys01),na.rm=T)/5)*5
            lowtick <- floor(min(c(Exp101,ys01),na.rm=T)/5)*5
            if(length(iAllNA01)>0) ys01[iAllNA01] <- toptick+5
            plot(x=Exp101,y=ys01[-iAllNA1],pch=c(19,17,15,18,4)[iNA01+1][-iAllNA1],
                 col=ifelse(apply(qc01>qcThreshold,1,all)[-iAllNA1],"black","red"),
                 main="0.01/0.01 Dilution",
                 xlab="Expected Expression\n(based on 0.1/0.1)",
                 ylab="Average Observed Expression",yaxt="n",
                 ylim=c(lowtick,toptick+5),xlim=c(lowtick,toptick))
            axis(side=2,at=seq(lowtick,toptick+5,by=5),
                 labels=c(as.character(seq(lowtick,toptick,by=5)),"NA"),las=2)
            abline(a=0,b=1,lwd=2,lty=3,col="blue")
            legend("topleft",pch=c(19,17,15,18,4),title="Proporiton Poor Quality",
                   c("0/4","1/4","2/4","3/4","4/4"),ncol=2,cex=0.75,
                   col=c("black",rep("red",4)))
        } else{
          ## 0.1/0.1 vs Expected (pure)
            toptick <- ceiling(max(M1,na.rm=T)/5)*5
            lowtick <- floor(min(M1,na.rm=T)/5)*5
            if(length(iAllNA1)>0) M1[iAllNA1] <- toptick+5
            plot(x=Exp1,y=M1,pch=c(19,17,15,18,4)[iNA1+1],
                 col=ifelse(apply(qc1>qcThreshold,1,all),"black","red"),
                 main="0.1/0.1 Dilution",
                 xlab="Expected Expression\n(based on pure samples)",
                 ylab="Expression Difference (Observed - Expected)", yaxt="n",
                 ylim=c(lowtick,toptick+5))
            axis(side=2,at=seq(lowtick,toptick+5,by=5),
                 labels=c(as.character(seq(lowtick,toptick,by=5)),"NA"),las=2)
            abline(a=0,b=1,lwd=2,lty=3,col="blue")
            legend("topleft",pch=c(19,17,15,18,4),title="Proporiton Poor Quality",
                   c("0/4","1/4","2/4","3/4","4/4"),ncol=2,cex=0.75,
                   col=c("black",rep("red",4)))

            ## 0.01/0.01 vs Expected (pure)
            toptick <- ceiling(max(M01,na.rm=T)/5)*5
            lowtick <- floor(min(M01,na.rm=T)/5)*5
            if(length(iAllNA01)>0) M01[iAllNA01] <- toptick+5
            plot(x=Exp01,y=M01,pch=c(19,17,15,18,4)[iNA01+1],
                 col=ifelse(apply(qc01>qcThreshold,1,all),"black","red"),
                 main="0.01/0.01 Dilution",
                 xlab="Expected Expression\n(based on pure samples)",
                 ylab="Expression Difference (Observed - Expected)", yaxt="n",
                 ylim=c(lowtick,toptick+5))
            axis(side=2,at=seq(lowtick,toptick+5,by=5),
                 labels=c(as.character(seq(lowtick,toptick,by=5)),"NA"),las=2)
            abline(a=0,b=1,lwd=2,lty=3,col="blue")
            legend("topleft",pch=c(19,17,15,18,4),title="Proporiton Poor Quality",
                   c("0/4","1/4","2/4","3/4","4/4"),ncol=2,cex=0.75,
                   col=c("black",rep("red",4)))
            
            ## 0.01/0.01 vs Expected (0.1/0.1)
            toptick <- ceiling(max(M101,na.rm=T)/5)*5
            lowtick <- floor(min(M101,na.rm=T)/5)*5
            if(length(iAllNA01[-iAllNA1])>0) M101[iNA01[-iAllNA1]==4] <- toptick+5
            plot(x=Exp101,y=M101,pch=c(19,17,15,18,4)[iNA01+1][-iAllNA1],
                 col=ifelse(apply(qc01>qcThreshold,1,all)[-iAllNA1],"black","red"),
                 main="0.01/0.01 Dilution",
                 xlab="Expected Expression\n(based on 0.1/0.1)",
                 ylab="Expression Difference (Observed - Expected)", yaxt="n",
                 ylim=c(lowtick,toptick+5))
            axis(side=2,at=seq(lowtick,toptick+5,by=5),
                 labels=c(as.character(seq(lowtick,toptick,by=5)),"NA"),las=2)
            abline(a=0,b=1,lwd=2,lty=3,col="blue")
            legend("topleft",pch=c(19,17,15,18,4),title="Proporiton Poor Quality",
                   c("0/4","1/4","2/4","3/4","4/4"),ncol=2,cex=0.75,
                   col=c("black",rep("red",4)))
        }
        
        return(round(out,digits=1))
    }
}
       

    
