
################################################################
#
# LPE functions with modifications
#
################################################################
fixbound.predict.smooth.spline <- function (object, x, deriv = 0)
{

    if (missing(x)) {
        if (deriv == 0) {
            return(object[c("x", "y")])
        }
        else {
            x <- object$x
        }
    }
    if (is.null(object)) {
        stop("not a valid smooth.spline object")
    }
    else {
        out <- predict(object, x, deriv)
        maxpredY <- object$y[object$x == max(object$x, na.rm=TRUE)]
        out$y[out$x > max(object$x, na.rm=TRUE)] <- maxpredY
        minpredY <- object$y[object$x == min(object$x, na.rm=TRUE)]
        out$y[out$x < min(object$x, na.rm=TRUE)] <- minpredY
        invisible(out)

    }

}


permut <- function (a)
{
    aa <- matrix(NA, length(a) - 1, length(a))
    for (i in 1:(length(a) - 1)) {
        aa[i, ] <- a[c((i + 1):length(a), 1:i)]
    }
    return(aa)
}


am.tran <- function (y) 
{
    A <- c()
    M <- c()
    n <- ncol(y)

    if (n < 2) {     
        stop("There are no replicated arrays!")
    } else
    {
       cc <- permut(1:n)
       for (i in 1:(n - 1)) {
           A <- c(A, c((y + y[, cc[i, ]])/2), recursive = TRUE)
           M <- c(M, c(y - y[, cc[i, ]]), recursive = TRUE)
       }
    }

    AM <- cbind(A,M)
    AM <- na.omit(AM)

    return(AM)
}


base.error.Olig <- function (y, q = 0.01, r = rep(1,ncol(y)), type="total") 
{

    nrep <- unique(r)  
    n <- length(nrep)
    nchip <- length(r)

    if(type=="exp") {                            #exp variation

        AM <-matrix(NA,1,2)
        for (i in 1:n) {
             j <- which(r==nrep[i])
             if(length(j) >1) AM <- rbind(AM,am.tran(y[,j]))
        }
        if(nrow(AM) <= 1) stop("No experimental replication")
        AM <- AM[-1,] 

    } else                                         #total variation
    {
        AM <-matrix(NA,1,2)
        if (n > 1){
            for (i in 1:(nchip-1)) {
                 for (j in (i+1):nchip) {
                      k<-c(i,j)
                      if(r[i]!=r[j]) {
                         #print(k)
                         AM <- rbind(AM,am.tran(y[,k]))
                      }
                 }
            }
            AM <- AM[-1,] 

        } else
        {
            AM <- am.tran(y)
        }
    }


    A <- AM[, 1]
    M <- AM[, 2]
    quantile.A <- quantile(A, probs = seq(0, 1, q), na.rm = TRUE)
    quan.n <- length(quantile.A) - 1
    var.M <- rep(NA, length = quan.n)
    medianAs <- rep(NA, length = quan.n)
 
    if (sum(A == min(A)) > (q * length(A))) {
        tmpA <- A[!(A == min(A))]
        quantile.A <- c(min(A), quantile(tmpA, probs = seq(q, 
            1, q), na.rm = TRUE))
    }
 
   for (i in 1:quan.n) {
        bin <- which(A > quantile.A[i] & A <= quantile.A[i+1])
        n.i <- length(!is.na(M[bin]))
        mult.factor <- 0.5 * ((n.i - 0.5)/(n.i - 1))
        var.M[i] <- mult.factor * var(M[bin], na.rm = TRUE)
        medianAs[i] <- median(A[bin], na.rm = TRUE)
    }
    return(cbind(A = medianAs, var.M = var.M, quantile.A = quantile.A[-1]))
}



boot.base.error.Olig  <- function (y, quantile.A, r = rep(1,ncol(y)))  
{


    nrep <- unique(r)  
    AM <-matrix(NA,1,2)
    for (i in 1:length(nrep)) {
         j <- which(r==nrep[i])
         if(length(j) >1) AM <- rbind(AM,am.tran(y[,j]))
    }
    if(nrow(AM) <= 1) stop("No experimental replication")
    AM <- AM[-1,]  

    A <- AM[, 1]
    M <- AM[, 2]
    quan.n <- length(quantile.A)
    var.M <- rep(NA, length = quan.n)
    #medianAs <- rep(NA, length = quan.n)


    bin <- which(A <= quantile.A[1])
    n.i <- length(!is.na(M[bin]))
    if(n.i > 1) {
       mult.factor <- 0.5 * ((n.i - 0.5)/(n.i - 1))
       var.M[1] <- mult.factor * var(M[bin], na.rm = TRUE)
       #medianAs[1] <- median(A[bin], na.rm = TRUE)
    }


    for (i in 2:quan.n) {
        bin <- which(A > quantile.A[i-1] & A <= quantile.A[i])
        n.i <- length(!is.na(M[bin]))
        if(n.i > 1) {
           mult.factor <- 0.5 * ((n.i - 0.5)/(n.i - 1))
           var.M[i] <- mult.factor * var(M[bin], na.rm = TRUE)
           #medianAs[i] <- median(A[bin], na.rm = TRUE)
        }
    }

    return(var.M)

}
