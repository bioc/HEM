
########################################################################################### 
#
# Computes HEM-F scores under null by re-sampling
#
########################################################################################### 

hem.null.one <- function (dat, n.layer, design, burn.ins = 1000, n.samples = 3000,
                          method.var.t = "gam", var.t = NULL, var.tot = NULL, 
                          var.g = 1, var.c = 1, var.r = 1, alpha.t = 3, beta.t = 0.2,
                          q=0.01, B=25, n.iter=5, 
                          print.message.on.screen=TRUE)
{

      
        method.vart <- ifelse(method.var.t=="gam", 1, ifelse(method.var.t=="peb", 2, 3))

        ###size
        n.gene <- nrow(dat)
        n.chip <- ncol(dat)
        cond <- design[,1]
        u.cond <- unique(cond)
        n.cond <- length(u.cond)
        c.size <- table(cond)

        ###Obatin within-gene-condition means
        mean.dat <- matrix(NA, n.gene, n.cond)
        for(j in 1:n.cond) {
            jj <- which(u.cond[j]==cond)
            mean.dat[,j] <- apply(dat[,jj], 1, mean)
        }


        ###Sort and obtain sigma
        j <-1
        jj <- which(u.cond[j]==cond)
        ii <- sort.list(mean.dat[,j])
        mean.dat[,j] <- mean.dat[ii,j]
        dat[,jj] <- dat[ii,jj]      

        sigma <- rep(sqrt(beta.t/alpha.t), nrow(dat))
        if(method.vart >1) {
           nn <- length(unique(design[jj,2]))
           sigma <- sqrt(var.tot[ii,1]/nn)#NOTE
        }

        for(j in 2:n.cond) {
            jj <- which(u.cond[j]==cond)
            ii <- sort.list(mean.dat[,j])
            mean.dat[,j] <- mean.dat[ii,j]
            dat[,jj] <- dat[ii,jj]
        }


        ###generate null data and apply HEM to the null data
        F.null <- matrix(NA, n.gene, n.iter)
        null.data <- dat
        for(it in 1:n.iter) {
            if(print.message.on.screen) cat(".")
            #null.data <- matrix(sample(dat), nrow=n.gene, ncol=n.chip) #full permutation

            for (i in 1:n.gene) { 
                   mu1 <- rnorm((n.cond-1), mean.dat[i,1], sigma[i])
                   kk <- 0
                   while(min(mu1) < min(mean.dat[,1]) | max(mu1) > max(mean.dat[,1]))
                   {
                         mu1 <- rnorm((n.cond-1), mean.dat[i,1], sigma[i])
                         if(kk > 100) break
                         kk <- kk+1
                   }

                   for (j in 2:n.cond) {
                        jj <- which(u.cond[j]==cond)
                        mm <- min(abs(mean.dat[,j]-mu1[j-1]))
                        i2 <- min(which(abs(mean.dat[,j]-mu1[j-1])==mm))                         
                        null.data[i,jj] <- dat[i2,jj]
                   }
            }

  

            # write.table(null.data, file = "null.csv", append = FALSE, quote = FALSE, sep = ",",
            # row.names = FALSE,col.names = FALSE)

            ##apply HEM to null data
            eb.out <- NULL
            if(method.vart > 1) {
               eb.out <- hem.eb.prior(null.data,  n.layer=n.layer, design=design,  
                                   method.var.t=method.var.t, q=q, B=B,
                                   print.message.on.screen=FALSE)              
            }


            F.null[,it] <- hem(null.data, n.layer=n.layer, design=design, 
                               burn.ins = burn.ins, n.samples = n.samples,
                               method.var.t = method.var.t, var.t = eb.out$var.t,
                               var.g = var.g, var.c = var.c, var.r = var.r,
                               alpha.t = alpha.t, beta.t = beta.t,
                               print.message.on.screen=FALSE)$F
        }
 
        return(F.null)

}




