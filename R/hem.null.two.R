
########################################################################################### 
#
# Computes HEM-F scores under null by re-sampling
#
########################################################################################### 


hem.null.two <- function (dat, n.layer, design, burn.ins = 1000, n.samples = 3000,
                method.var.e = "gam", method.var.b = "gam", 
                var.e = NULL, var.b = NULL, var.exp = NULL, var.bio = NULL, var.tot= NULL,
                var.g = 1, var.c = 1, var.r = 1, 
                alpha.e = 3, beta.e = 0.1, alpha.b = 3, beta.b = 0.1,
                q=0.01, B=25,  n.iter=5,
                print.message.on.screen=TRUE)
{


        method.vare <- ifelse(method.var.e=="gam", 1, ifelse(method.var.e=="peb", 2, 3))
        method.varb <- ifelse(method.var.b=="gam", 1, 2)

        ###Size
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
        sigma <- rep(sqrt(beta.e/alpha.e + beta.b/alpha.b), nrow(dat))

        if(method.vare >1) {
           mm <- length(unique(design[jj,2]))
           jj1<-jj[1]  
           nn <- length(which(design[jj,2]==design[jj1,2]))
           sigma <- sqrt(var.bio[ii,1]/mm+var.exp[ii,1]/nn)
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

             #null.full.data <- matrix(sample(dat), nrow=n.gene, ncol=n.chip) #full permutation
             #write.table(null.full.data, file = "null.full.csv", append = FALSE, quote = FALSE,
             # sep = ",", row.names = FALSE, col.names = FALSE)

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
  

             #write.table(null.data, file = "null.csv", append = FALSE, quote = FALSE, sep = ",",
             #row.names = FALSE, col.names = FALSE)


            ###apply HEM to null data
            lpe.out <- NULL
            if(method.vare >1 & method.varb >1) {
               lpe.out <- hem.eb.prior(null.data,  n.layer=n.layer, design=design, 
                                   method.var.e=method.var.e, method.var.b=method.var.b, 
                                   q=q, B=B, print.message.on.screen=FALSE)
            }

            F.null[,it] <- hem(null.data, n.layer=n.layer, design=design, 
                               burn.ins = burn.ins, n.samples = n.samples,
                               method.var.e = method.var.e, method.var.b = method.var.b,
                               var.e = lpe.out$var.e, var.b = lpe.out$var.b, 
                               var.g = var.g, var.c = var.c, var.r = var.r,
                               alpha.e = alpha.e, beta.e = beta.e, 
                               alpha.b = alpha.b, beta.b = beta.b, 
                               print.message.on.screen=FALSE)$F
       }

       return(F.null)

}

















