## assess precision 
precision <- function(object1, qcThreshold1,
                      object2=NULL, qcThreshold2=NULL,
                      commonFeatures=TRUE,
                      statistic=c("sd","cv"),
                      scale=c("none","log","log10"),
                      bins=3, label1=NULL, label2=NULL){

    scale <- match.arg(scale)
    statistic <- match.arg(statistic)

    object1 <- checkObject(object1)
    cts1 <- object1$ct
    qc1 <- object1$qc
    xs <- as.numeric(gsub("KW","",gsub(":.+","",colnames(cts1))))

    if(is.null(object2)){
        ## filter miRNAs with any NAs or poor quality values
        i.rm <- which(apply(qc1<qcThreshold1,1,any) | apply(is.na(cts1),1,any))
        cts1 <- cts1[-i.rm,]
    
        ## calculate within condition means and variances
        mu <- as.vector(apply(cts1,1,by,xs,mean))
        s <- as.vector(apply(cts1,1,by,xs,sd))

        if(statistic=="cv") stat <- s/mu else stat <- s

        mubins <- cut(mu, breaks=quantile(mu, probs=seq(0,1,length.out=bins+1))
                      -c(0.1,rep(0,bins-1),-0.1))
        statbins <- split(stat, mubins)
        ylabel <- ifelse(statistic=="cv","Coefficient of Variation",
                         "Standard Deviation")
        if(scale=="log"){
          statbins <- lapply(statbins,log)
          ylabel <- paste(ylabel,"(log)")
        }
        if(scale=="log10"){
          statbins <- lapply(statbins,log10)
          ylabel <- paste(ylabel,"(log10)")
        }
        boxplot(statbins[bins:1], ylab=ylabel)
        if(bins>1){
          if(bins==2){
            axis(side=1,labels=c("Low","High"),at=c(1,bins), tick=FALSE, 
                 line=2)
          } else{
            axis(side=1,labels=c("Low","Medium","High"),
                 at=c(1,(bins+1)/2,bins), tick=FALSE, line=2)
          }
        }
    } else{
        object2 <- checkObject(object2)
        cts2 <- object2$ct
        qc2 <- object2$qc

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
    
        ## calculate within condition means and variances
        mu1 <- as.vector(apply(cts1,1,by,xs,mean))
        s1 <- as.vector(apply(cts1,1,by,xs,sd))
        mu2 <- as.vector(apply(cts2,1,by,xs,mean))
        s2 <- as.vector(apply(cts2,1,by,xs,sd))

        if(statistic=="cv") stat1 <- s1/mu1 else stat1 <- s1
        if(statistic=="cv") stat2 <- s2/mu2 else stat2 <- s2
        
        mubins1 <- cut(mu1, breaks=quantile(mu1,
                              probs=seq(0,1,length.out=bins+1))
                      -c(0.1,rep(0,bins-1),-0.1))
        mubins2 <- cut(mu2, breaks=quantile(mu2,
                              probs=seq(0,1,length.out=bins+1))
                      -c(0.1,rep(0,bins-1),-0.1))
        statbins1 <- split(stat1, mubins1)
        statbins2 <- split(stat2, mubins2)
        ylabel <- ifelse(statistic=="cv","Coefficient of Variation",
                         "Standard Deviation")
        if(scale=="log"){
            statbins1 <- lapply(statbins1,log)
            statbins2 <- lapply(statbins2,log)
            ylabel <- paste(ylabel,"(log)")
        }
        if(scale=="log10"){
            statbins1 <- lapply(statbins1,log10)
            statbins2 <- lapply(statbins2,log10)
            ylabel <- paste(ylabel,"(log10)")
        }
        statbins <- list()
        for(k in bins:1){
            i2 <- 2*k
            i1 <- i2-1
            statbins[[i1]] <- statbins1[[k]]
            names(statbins)[i1] <- names(statbins1)[k]
            statbins[[i2]] <- statbins2[[k]]
            names(statbins)[i2] <- names(statbins2)[k]
        }
        tmp <- names(statbins)
        if(is.null(label1)) label1 <- "M1"
        if(is.null(label2)) label2 <- "M2"
        names(statbins) <- rep(c(label1,label2),bins)
        boxplot(statbins, ylab=ylabel)
        if(bins>1){
          if(bins==2){
            axis(side=1, labels=c("Low","High"), at=c(1.5,2*bins-0.5), 
                 tick=FALSE, line=2)
          } else{
            axis(side=1,labels=c("Low","Medium","High"),
                 at=c(1.5,bins+0.5,2*bins-0.5), tick=FALSE, line=2)
          }
        }
        names(statbins) <- paste(names(statbins),tmp,sep=":")
    }
    return(statbins)
}
