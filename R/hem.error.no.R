############################################################
### To estimate total error variance by LPE (parametirc) ###
############################################################
par.no.error.Olig <- function(x, q=0.01,  df=10, baseline.var="BLPE", p.remove=0)
{ 

      median.x <- x
      var.total <- x

      if(p.remove >0) {
         basevar.x <- base.error.Olig.quanOnly(x, q = q, r = (1:ncol(x)), type="total")
         quan.A <- basevar.x[,3]
         quan.n <- length(quan.A)
         x <- remove.sig.genes(x, q=q, quantile.A=quan.A, p.remove=p.remove)
      } 
      if(baseline.var=="BLPE") basevar.total <- base.error.Olig(x, q = q, r=(1:ncol(x)), type="total")
      if(baseline.var=="PSE") basevar.total <- base.PSE.Olig(x, q = q)
      if(baseline.var=="ASE") basevar.total <- base.ASE.Olig(x, q = q)

      sf.x <- smooth.spline(basevar.total[,1], basevar.total[,2], df = df)
      max.A <- max(basevar.total[,1])
      min.A <- min(basevar.total[,1])

      for(j in 1:ncol(median.x)){

          k <- which(median.x[,j] > max.A); median.x[k,j] <- max.A 
          k <- which(median.x[,j] < min.A); median.x[k,j] <- min.A
 
          var.total[,j] <- fixbound.predict.smooth.spline(sf.x, median.x[,j])$y

          #Correction on low intensities
          i.max <- max(which(var.total[,j]==max(var.total[,j])))      
          x.max <- median.x[i.max,j]
          var.total[which(x[,j] < x.max),j] <- max(var.total[,j])
      }

      return(list(m=x, var.total=var.total))
}


################################################################
### To estimate total error variance by BLPE (nonparametric) ###
################################################################
nonpar.no.error.Olig <- function(x, q=0.01, B=100, baseline.var="BLPE", p.remove=0,
                                 print.message.on.screen=TRUE)
{ 

       r=(1:ncol(x)) # 1,2,3,... in one-layer HEM.
       basevar.x <- base.error.Olig.quanOnly(x, q = q, r = r, type="total")
       quan.A <- basevar.x[,3]
       quan.n <- length(quan.A)

       r = rep(1,ncol(x)) # 1,1,1,.... total error is like exp error
       basevar <- matrix(NA,nrow=(B+1), ncol=quan.n) 


       for (i in 1:B) {
          bootstrap.row <- sample(nrow(x), replace=TRUE)
          data.generated <- x[bootstrap.row,]
          if(p.remove >0) data.generated <- remove.sig.genes(data.generated, q=q, quantile.A=quan.A, p.remove=p.remove) 
          if(baseline.var=="BLPE") result <- boot.base.error.Olig(data.generated, quantile.A=quan.A, r=r)
          if(baseline.var=="PSE")  result <- boot.base.PSE.Olig(data.generated, quantile.A=quan.A)
          if(baseline.var=="ASE")  result <- boot.base.ASE.Olig(data.generated, quantile.A=quan.A)
          basevar[i,] <- t(result)
          if(print.message.on.screen) cat(".")
          #cat("ITERATION  = ", i, " FINISHED \n")
       }

       #Correction on low intensities
       x <- basevar.x[,3]
       v <- basevar.x[,2]
       i.max <- max(which(v==max(v)))      
       x.max <- x[i.max]
       i <- which(x < x.max)
       if(length(i) >=1) basevar[,i] <- basevar[,i.max] 

       basevar[(B+1),] <- quan.A 
       return(basevar)
}

################################################################
# Remove significant genes
################################################################
remove.sig.genes <- function(x, q=0.01, quantile.A, p.remove=0.5) {

    x.rank <- apply(x, 2, rank)
    max.rank <- apply(x.rank, 1, max) 
    min.rank <- apply(x.rank, 1, min) 
    diff.rank <- abs(max.rank-min.rank)
    keep.genes <- c()

    A <- apply(x, 1, mean)
    quan.n <- length(quantile.A)

    bin <- which(A <= quantile.A[1])
    n.i <- length(!is.na(diff.rank[bin]))
    if(n.i > 1) {
        I <- which(diff.rank <= quantile(diff.rank[bin], probs=(1-p.remove)))
        keep.genes <- c(keep.genes, intersect(bin,I))
    }


    for (i in 2:quan.n) {
        bin <- which(A > quantile.A[i-1] & A <= quantile.A[i])
        n.i <- length(!is.na(diff.rank[bin]))
        if(n.i > 1) {
           I <- which(diff.rank <= quantile(diff.rank[bin], probs=(1-p.remove)))
           keep.genes <- c(keep.genes, intersect(bin,I))
        }
     }

     x.sub <- x[keep.genes,]
     return(x.sub)

}

################################################################
# Estimate baseline error by PSE
################################################################
boot.base.PSE.Olig  <- function (y, quantile.A)  
{

    AM <- am.tran.half(y)
    A <- AM[,1]
    M <- abs(AM[,2])
    quan.n <- length(quantile.A)
    var.M <- rep(NA, length = quan.n)


    bin <- which(A <= quantile.A[1])
    n.i <- length(!is.na(M[bin]))
    if(n.i > 1) {
        s0 <- 1.5 * median(M[bin], na.rm = TRUE)
        I <- which(M <= 2.5*s0)
        var.M[1] <- 1.5 * median(M[intersect(bin,I)])
        var.M[1] <- var.M[1]*var.M[1]/2.0 #fixed 1-17-05    
}



    for (i in 2:quan.n) {
        bin <- which(A > quantile.A[i-1] & A <= quantile.A[i])
        n.i <- length(!is.na(M[bin]))
        if(n.i > 1) {
           s0 <- 1.5 * median(M[bin], na.rm = TRUE)
           I <- which(M <= 2.5*s0)
           var.M[i] <- 1.5 * median(M[intersect(bin,I)])
           var.M[i] <- var.M[i]*var.M[i]/2.0 #fixed 1-17-05
        }
    }

    return(var.M)

}

################################################################
# Estimate baseline error by ASE
################################################################
boot.base.ASE.Olig  <- function (y, quantile.A)  
{

    AM <- am.tran.half(y)
    A <- AM[,1]
    M <- abs(AM[,2])
    quan.n <- length(quantile.A)
    var.M <- rep(NA, length = quan.n)


    bin <- which(A <= quantile.A[1])
    n.i <- length(!is.na(M[bin]))
    if(n.i > 1) {
        s0 <- 1.5 * median(M[bin], na.rm = TRUE)
        I <- which(M <= 2.5*s0)
        MM <- M[intersect(bin,I)]
        var.M[1] <- sum(MM * MM) / length(MM)
        var.M[1] <- var.M[1]/2.0 #fixed 1-17-05
    }


    for (i in 2:quan.n) {
        bin <- which(A > quantile.A[i-1] & A <= quantile.A[i])
        n.i <- length(!is.na(M[bin]))
        if(n.i > 1) {
           s0 <- 1.5 * median(M[bin], na.rm = TRUE)
           I <- which(M <= 2.5*s0)
           MM <- M[intersect(bin,I)]
           var.M[i] <- sum(MM * MM) / length(MM)
           var.M[i] <- var.M[i]/2.0 #fixed 1-17-05
        }
    }

    return(var.M)

}


################################################################
# Obtain AM values 
################################################################
am.tran.half <- function (y, max.chip=4) #default??? 
{
    A <- c()
    M <- c()
    n <- ncol(y)

    if (n < 2) stop("There are no replicated arrays!")
    if (n > max.chip) {
       j <- sample(1:n, size=max.chip)
       y <- y[,j]  
       n <- ncol(y)
    }

       
    for (i in 1:(n-1)) {
         for (j in (i+1):n) {
              A <- c(A, c((y[,i] + y[,j])/2))
              M <- c(M, c( y[,i] - y[,j]))
         }
    }

    AM <- cbind(A,M)
    AM <- na.omit(AM)
    return(AM)
}

################################################################
# Obtain quantiles
################################################################
base.error.Olig.quanOnly <- function (y, q = 0.01, r = rep(1,ncol(y)), type="total") 
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
        if (n== 1) AM <- am.tran(y)
        if (n > 1){
            AM <-matrix(NA,1,2)
            for (i in 1:(nchip-1)) {
                 for (j in (i+1):nchip) {
                      k <- c(i,j)
                      if(r[i]!=r[j]) {
                         #print(k)
                         AM <- rbind(AM,am.tran(y[,k]))
                      }
                 }
            }
            AM <- AM[-1,] 
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




################################################################
# Estimate baseline error by PSE
################################################################
base.PSE.Olig <- function (y, q = 0.01) 
{


    AM <- am.tran.half(y)
    A <- AM[,1]
    M <- abs(AM[,2])

    quantile.A <- quantile(A, probs = seq(0, 1, q), na.rm = TRUE)
    quan.n <- length(quantile.A) - 1
    var.M <- rep(NA, length = quan.n)
    medianAs <- rep(NA, length = quan.n)

    for (i in 1:quan.n) {
           bin <- which(A > quantile.A[i] & A <= quantile.A[i+1])
           s0 <- 1.5 * median(M[bin], na.rm = TRUE)
           I <- which(M <= 2.5*s0)
           var.M[i] <- 1.5 * median(M[intersect(bin,I)])
           var.M[i] <- var.M[i]*var.M[i]/2.0 #fixed 1-17-05
           medianAs[i] <- median(A[intersect(bin,I)], na.rm = TRUE)
    }

    return(cbind(A = medianAs, var.M = var.M, quantile.A = quantile.A[-1]))
}

################################################################
# Estimate baseline error by ASE
################################################################
base.ASE.Olig <- function (y, q = 0.01) 
{


    AM <- am.tran.half(y)
    A <- AM[,1]
    M <- abs(AM[,2])

    quantile.A <- quantile(A, probs = seq(0, 1, q), na.rm = TRUE)
    quan.n <- length(quantile.A) - 1
    var.M <- rep(NA, length = quan.n)
    medianAs <- rep(NA, length = quan.n)


    for (i in 1:quan.n) {
           bin <- which(A > quantile.A[i] & A <= quantile.A[i+1])
           s0 <- 1.5 * median(M[bin], na.rm = TRUE)
           II <- which(M <= 2.5*s0)
           MM <- M[intersect(bin,II)]
           AA <- A[intersect(bin,II)]
           var.M[i] <- sum(MM * MM) / length(MM)
           var.M[i] <- var.M[i]/2.0 #fixed 1-17-05
           medianAs[i] <- median(AA, na.rm = TRUE)
    }
    return(cbind(A = medianAs, var.M = var.M, quantile.A = quantile.A[-1]))
}

