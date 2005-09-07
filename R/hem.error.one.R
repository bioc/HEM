############################################################
### To estimate total error variance by LPE (parametirc) ###
############################################################
par.error.Olig <- function(x, q=0.01,  df=10)
{ 

      r=(1:ncol(x)) # 1,2,3,... in one-layer HEM
      median.x <- apply(x, 1, median, na.rm=TRUE)
      basevar.total <- base.error.Olig(x, q = q, r=r, type="total")
      sf.x <- smooth.spline(basevar.total[,1], basevar.total[,2], df = df)

      var.total <- fixbound.predict.smooth.spline(sf.x, median.x)$y

      #Correction on low intensities
      i.max <- max(which(var.total==max(var.total)))      
      x.max <- median.x[i.max]
      var.total[which(median.x < x.max)] <- max(var.total)

      return(list(m=median.x, var.total=var.total))
}

###############################################################
### To estimate total error variance by LPE (nonparametric) ###
###############################################################
nonpar.error.Olig <- function(x, q=0.01, B=100, print.message.on.screen=TRUE)
{ 

       r=(1:ncol(x)) # 1,2,3,... in one-layer HEM.
       basevar.x <- base.error.Olig(x, q = q, r = r, type="total")
       quan.A <- basevar.x[,3]
       quan.n <- length(quan.A)

       r = rep(1,ncol(x)) # 1,1,1,.... total error is like exp error
       basevar <- matrix(NA,nrow=(B+1), ncol=quan.n) 


      for (i in 1:B) {

          bootstrap.row <- sample(nrow(x), replace=TRUE)
          data.generated <- x[bootstrap.row,]
          result <- boot.base.error.Olig(data.generated,quantile.A=quan.A,r=r)
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

