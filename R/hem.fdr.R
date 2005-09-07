
########################################################################################### 
#
# Computes FDR 
# Updated on 2005.07.28
#
########################################################################################### 
hem.fdr <- function(dat, n.layer, design, rep=TRUE,
           hem.out, eb.out=NULL, n.iter=5, q.trim=0.9,
           target.fdr=c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.20, 0.30, 0.40, 0.50), 
           n.digits=10, print.message.on.screen=TRUE)
{

        if(print.message.on.screen) cat("Computing FDR...Please wait.\n")

        ###design matrix
        ncol.design <- ncol(design)
        if(ncol.design > 3 | ncol.design < 2 ) stop("Error: incorrect size of design") 
        for (j in ncol.design) {
             u <- unique(design[,j])
             for (i in 1:length(u)){
                  w <- which(design[,j]==u[i])
                  design[w,j] <- i 
             }
        }
        #design <- as.numeric(design)
        if(ncol.design==3) i <- order(design[,1],design[,2],design[,3])
        else  i <- order(design[,1],design[,2])
        dat <- dat[,i]
        design <- design[i,]


        #####################################################################
        ###compute H-scores under null with iterations

        if(n.layer==2) {
               H.null <- hem.null.two(dat,
                         n.layer      = n.layer, 
                         design       = design,
                         burn.ins     = hem.out$burn.ins, 
                         n.samples    = hem.out$n.samples,
                         method.var.e = hem.out$method.var.e, 
                         method.var.b = hem.out$method.var.b, 
                         var.g        = hem.out$priors[1], 
                         var.c        = hem.out$priors[2], 
                         var.r        = hem.out$priors[3], 
                         alpha.e      = hem.out$priors[4], 
                         beta.e       = hem.out$priors[5],
                         alpha.b      = hem.out$priors[6], 
                         beta.b       = hem.out$priors[7],
                         var.e        = eb.out$var.e, 
                         var.b        = eb.out$var.b, 
                         var.exp      = eb.out$var.exp, ###NOTE
                         var.bio      = eb.out$var.bio, ###NOTE
                         var.tot      = eb.out$var.tot, ###NOTE
                         q            = eb.out$q, 
                         B            = eb.out$B,
                         n.iter       = n.iter,  
                         print.message.on.screen=print.message.on.screen)
        }


        if(n.layer==1 & rep) {
               H.null <- hem.null.one(dat,
                         n.layer      = n.layer, 
                         design       = design,
                         burn.ins     = hem.out$burn.ins, 
                         n.samples    = hem.out$n.samples,
                         method.var.t = hem.out$method.var.t, 
                         var.g        = hem.out$priors[1], 
                         var.c        = hem.out$priors[2], 
                         var.r        = hem.out$priros[3], 
                         alpha.t      = hem.out$priors[4], 
                         beta.t       = hem.out$priors[5],
                         q            = eb.out$q, 
                         B            = eb.out$B,                        
                         var.t        = eb.out$var.t,
                         var.tot      = eb.out$var.tot, ###NOTE 
                         n.iter       = n.iter,  
                         print.message.on.screen=print.message.on.screen)
        }

        if(n.layer==1 & !rep) {
               H.null <- hem.null.no(dat,
                         n.layer      = n.layer, 
                         design       = design,
                         burn.ins     = hem.out$burn.ins, 
                         n.samples    = hem.out$n.samples,
                         method.var.t = hem.out$method.var.t, 
                         var.g        = hem.out$priors[1], 
                         var.c        = hem.out$priors[2], 
                         var.r        = hem.out$priros[3], 
                         alpha.t      = hem.out$priors[4], 
                         beta.t       = hem.out$priors[5],
                         baseline.var = eb.out$baseline.var, 
                         p.remove     = eb.out$p.remove,                        
                         q            = eb.out$q, 
                         B            = eb.out$B,                        
                         var.t        = eb.out$var.t,
                         var.tot      = eb.out$var.tot, ###NOTE 
                         n.iter       = n.iter,  
                         print.message.on.screen=print.message.on.screen)
        }


        #####################################################################
        ###compute FDR
        len.H <- length(hem.out$H[,1]) 
        H.null <- as.vector(H.null)
        min.fdr <- 1.0/(len.H*n.iter)

        m.c <- quantile(H.null, probs=q.trim)
        pi0 <- length(which(hem.out$H[,1] <= m.c))/(len.H*q.trim)
        pi0 <- min(1, pi0)

        i <- sort.list(hem.out$H[,1], decreasing=TRUE) 
        id <- 1:length(hem.out$H[,1]); id <- id[i]
        H <- hem.out$H[i,] 

        H.all <- c(H, H.null)
        i.all <- c(rep(0, length(hem.out$H[,1])), rep(1, length(H.null)))
        i <- sort.list(H.all, decreasing=TRUE)
        all <- data.frame(H.all[i], i.all[i])

        V <- 1.0*cumsum(all[,2])/n.iter; V <- V[which(all[,2]==0)]
        R <- 1:length(V)
        FDR <- V / R * pi0; FDR[FDR <= 0] <- min.fdr; FDR[FDR > 1] <- 1
        fdr <- round(data.frame(H, FDR),n.digits)

        ###find critical values at target FDRs
        n.target <- length(target.fdr)
        H.critical <- rep(NA, n.target)
        n.sig.genes <- rep(NA, n.target)
        min.FDR <- min(FDR)

        for(i in 1:n.target) {
            if(target.fdr[i] >= min.FDR){
               k <- max(which(fdr$FDR <= target.fdr[i]))
               H.critical[i] <- round(H[k],n.digits)
               n.sig.genes[i] <- length(which(hem.out$H[,1] >= H.critical[i]))
            }
        }

        target.FDR <- target.fdr
        targets <- data.frame(target.FDR, H.critical, n.sig.genes)

        #Put FDRs in the data order
        i <- sort.list(id)
        fdr <- fdr[i,]
        rownames(fdr) <- rownames(hem.out$H) #H CHO, 05-27-05
         

        if(print.message.on.screen) cat("Done.\n")
        return(list(fdr=fdr, pi0=pi0, q.trim=q.trim, n.iter=n.iter, H.null=H.null, targets=targets))

}

######END############################################################################
