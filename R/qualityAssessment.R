qualityAssessment <- function(object1, object2=NULL,
                              cloglog1=FALSE, cloglog2=FALSE,
                              na.rm=FALSE,
                              plotType=c("scatterplot", "boxplot"),
                              label1=NULL, label2=NULL){
    plotType <- match.arg(plotType)

    object1 <- checkObject(object1)
    if(cloglog1) object1$qc <- -log(-log(object1$qc))
    if(is.null(object2)){
        if(is.null(label1)) label1 <- "Quality Score"
        if(plotType=="scatterplot"){
            toptick <- ceiling(max(object1$ct, na.rm=T)/5)*5
            lowtick <- floor(min(object1$ct, na.rm=T)/5)*5
            naValue <- toptick+5
            object1$ct[is.na(object1$ct)] <- naValue
            ind <- which(object1$ct==naValue |
                             object1$qc==min(object1$qc, na.rm=T))
            smoothScatter(x=object1$ct[-ind], y=object1$qc[-ind], pch=19,
                          cex=0.25, xlab="Expression Estimate", ylab=label1, 
                          main="", xaxt="n", yaxt="n",
                          xlim=c(lowtick,naValue),
                          ylim=range(object1$qc, na.rm=T))
            points(x=object1$ct[ind], y=object1$qc[ind], pch=19, cex=0.25)
            axis(side=1, at=seq(lowtick, naValue, by=5),
                 labels=c(as.character(seq(lowtick, toptick, by=5)), "NA"))
            if(cloglog1){
                axis(side=2, at=-log(-log(c(0.01,0.5,0.9,0.99,0.999,0.9999))),
                     labels=c("0.01","0.50","0.90","0.99","0.999","0.9999"),
                     las=1)
            } else axis(side=2)
        }
        if(plotType=="boxplot"){
            if(na.rm) object1$qc[is.na(object1$ct)] <- NA
            boxplot(object1$qc, las=2, pch=20, ylab=label1, xlab="", 
                    yaxt="n")
            if(cloglog1){
                axis(side=2, at=-log(-log(c(0.01,0.5,0.9,0.99,0.999,0.9999))),
                     labels=c("0.01","0.50","0.90","0.99","0.999","0.9999"),
                     las=1)
            } else axis(side=2)
        }
    } else{
        object2 <- checkObject(object2)
        if(cloglog2) object2$qc <- -log(-log(object2$qc))
        if(is.null(label1)) label1 <- "Object 1 Quality Score "
        if(is.null(label2)) label2 <- "Object 2 Quality Score"
        if(plotType=="scatterplot"){
            smoothScatter(y=object2$qc, x=object1$qc, xaxt="n", yaxt="n",
                          xlab=label1, ylab=label2)
            if(cloglog1){
                axis(side=1, at=-log(-log(c(0.01,0.5,0.9,0.99,0.999,0.9999))),
                     labels=c("0.01","0.50","0.90","0.99","0.999","0.9999"),
                     las=1)
            } else axis(side=1)
            if(cloglog2){
                axis(side=2, at=-log(-log(c(0.01,0.5,0.9,0.99,0.999,0.9999))),
                     labels=c("0.01","0.50","0.90","0.99","0.999","0.9999"),
                     las=1)
            } else axis(side=2)
        }
        if(plotType=="boxplot"){
            if(na.rm){
                object1$qc[is.na(object1$ct)] <- NA
                object2$qc[is.na(object2$ct)] <- NA
            }
            boxplot(object1$qc, las=2, pch=20, yaxt="n", ylab=label1,
                    xlab="")
            if(cloglog1){
                axis(side=2, at=-log(-log(c(0.01,0.5,0.9,0.99,0.999,0.9999))),
                     labels=c("0.01","0.50","0.90","0.99","0.999","0.9999"),
                     las=1)
            } else axis(side=2, las=1)
            axis(side=1, at=1:ncol(object1$qc), labels=colnames(object1$qc),
                 las=2)
            boxplot(object2$qc, las=2, pch=20, yaxt="n", ylab=label2,
                    xlab="")
            if(cloglog2){
                axis(side=2, at=-log(-log(c(0.01,0.5,0.9,0.99,0.999,0.9999))),
                     labels=c("0.01","0.50","0.90","0.99","0.999","0.9999"),
                     las=1)
            } else axis(side=2, las=1)
            axis(side=1, at=1:ncol(object2$qc), labels=colnames(object2$qc),
                 las=2)
        }
    }
}
