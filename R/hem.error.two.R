######################################################
### To estimate biological error variance by LPE ###
######################################################
par.rep.error.Olig <- function(x, q=0.01, r=rep(1,ncol(x)), df=10)
{ 

      median.x <- apply(x, 1, median,na.rm=TRUE)


      #Experimental  variation
      basevar.exp <- base.error.Olig(x, q = q, r=r, type="exp")
      sf.x <- smooth.spline(basevar.exp[,1],basevar.exp[,2], df = df)
      var.exp <- fixbound.predict.smooth.spline(sf.x, median.x)$y

      #Total variation
      basevar.total <- base.error.Olig(x, q = q, r=r, type="total")
      sf.x <- smooth.spline(basevar.total[,1], basevar.total[,2], df = df)
      var.total <- fixbound.predict.smooth.spline(sf.x, median.x)$y

      #Correction on low intensities
      i.max <- max(which(var.total==max(var.total)))      
      x.max <- median.x[i.max]
      var.total[which(median.x < x.max)] <- max(var.total)

      i.max <- max(which(var.exp==max(var.exp)))      
      x.max <- median.x[i.max]
      var.exp[which(median.x < x.max)] <- max(var.exp)

      x <- basevar.exp[,1]
      v <- basevar.exp[,2]
      i.max <- which(v==max(v))      
      x.max <- x[i.max]
      v[which(x < x.max)] <- max(v)
      basevar.exp[,2] <- v


      #Biological variation
      var.bio <- var.total - var.exp
      
      basevar <- t(basevar.exp[,2:3])

      return(list(m=median.x, var.exp=var.exp, var.bio=var.bio, var.total=var.total, var.e=basevar))
      ###return(var.bio=var.bio, var.exp=basevar)   

}


######################################################
### To estimate experimental error variance by LPE ###
######################################################
nonpar.rep.error.Olig <- function(x, q=0.01, r=rep(1,ncol(x)), B=100, print.message.on.screen=TRUE)
{ 

       basevar.x <- base.error.Olig(x, q = q, r = r, type="exp")
       quan.A <- basevar.x[,3]
       quan.n <- length(quan.A)

       basevar <- matrix(NA,nrow=(B+1), ncol=quan.n) 
       for (i in 1:B) {

          bootstrap.row <- sample(nrow(x), replace=TRUE)
          data.generated <- x[bootstrap.row,]
          result <- boot.base.error.Olig(data.generated, quantile.A=quan.A, r=r)
          basevar[i,] <- t(result)
          if(print.message.on.screen) cat(".")
          #cat("ITERATION  = ", i, " FINISHED \n")

      }


      #Correction on low intensities
      x <- basevar.x[,3]
      v <- basevar.x[,2]
      i.max <- max(which(v==max(v)))+1      
      x.max <- x[i.max]
      basevar[,which(x < x.max)] <- basevar[,i.max] 


      #returning
      basevar[(B+1),] <- quan.A 
      return(basevar)
}

