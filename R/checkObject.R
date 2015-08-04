checkObject <- function(object){

    ## load lifetech data object to compare
    data("lifetech", package="miRcomp", envir=environment())
    lifetech <- get("lifetech", envir=environment())

    ## check object format
    if(length(object)!=2){
        stop("Object must be a list with two components: ct and qc.")
    }
    names(object) <- tolower(names(object))
    if(!identical(names(object),c("ct","qc"))){
        warning("Object names are not recognizable as ct and qc.")
        warning(paste("Setting",names(object)[1],"as the ct matrix."))
        warning(paste("Setting",names(object)[2],"as the qc matrix."))
        names(object) <- c("ct","qc")
    }

    ## check that ct and qc matrices have matching rows / cols
    if(!identical(rownames(object$ct),rownames(object$qc))){
        stop("Rownames of ct and qc objects must be identical.")
    }
    if(!identical(colnames(object$ct),colnames(object$qc))){
        stop("Colnames of ct and qc objects must be identical.")
    }

    ## check that ct and qc matrics have correct cols and order
    if(!identical(colnames(object$ct),colnames(lifetech$ct))){
        if(!all(colnames(object$ct)%in%colnames(lifetech$ct))){
            stop("Object does not contain all 40 samples.")
        }
        message("Attempting to reorder data columns to match required order.")
        map <- match(colnames(lifetech$ct),colnames(object$ct))
        object$ct <- object$ct[,map]
        object$qc <- object$qc[,map]
        if(!identical(colnames(object$ct),colnames(lifetech$ct))){
            stop("Not able to reorder data columns to match required order.")
        }
    }

    return(object)
}
        
    
