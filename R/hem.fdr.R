
########################################################################################### 
#
# Computes FDR 
#
########################################################################################### 
hem.fdr <- function(dat, tr=" ", n.layer, design, 
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


        if(tr==" ") {
            min.value <- -9
            if(any(is.na(dat))) dat[is.na(dat)] <- min.value
            if(!is.matrix(dat)) dat <- as.matrix(dat)
            dat[which(dat <= min.value, arr.ind=TRUE)] <- min.value
        }

        if(tr!=" "){
            min.value <- 0.001
            if(any(is.na(dat))) dat[is.na(dat)] <- min.value
            if(!is.matrix(dat)) dat <- as.matrix(dat)
            dat[which(dat <= min.value, arr.ind=TRUE)] <- min.value

            if(tr=="log2" ) dat=log2(dat)
            if(tr=="log10") dat=log10(dat)
            if(tr=="loge" ) dat=log(dat)
        }

        ###compute F under null with iterations
        if(n.layer==2) {
               F.null <- hem.null.two(dat,
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


        if(n.layer==1) {
               F.null <- hem.null.one(dat,
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



        ###compute FDR
        len.F <- length(hem.out$F) 
        F.null <- as.vector(F.null)
        min.fdr <- 1.0/(len.F*n.iter)

        m.c <- quantile(F.null, probs=q.trim)
        pi0 <- length(which(hem.out$F <= m.c))/(len.F*q.trim)
        pi0 <- min(1, pi0)

        F <- sort(hem.out$F, decreasing=TRUE) 
        F.all <- c(F, F.null)
        i.all <- c(rep(0, length(hem.out$F)), rep(1, length(F.null)))
        i <- sort.list(F.all, decreasing=TRUE)
        all <- data.frame(F.all[i], i.all[i])

        V <- 1.0*cumsum(all[,2])/n.iter; V <- V[which(all[,2]==0)]
        R <- 1:length(V)
        FDR <- V/R*pi0; FDR[FDR<=0] <- min.fdr; FDR[FDR> 1] <- 1
        fdr <- round(data.frame(F, FDR),n.digits)



        ###find critical values at target FDRs
        n.target <- length(target.fdr)
        F.critical <- rep(NA, n.target)
        n.sig.genes <- rep(NA, n.target)
        min.FDR <- min(FDR)

        for(i in 1:n.target) {
            if(target.fdr[i] >= min.FDR){
               k <- max(which(fdr$FDR <= target.fdr[i]))
               F.critical[i] <- round(F[k],n.digits)
               n.sig.genes[i] <- length(which(hem.out$F >= F.critical[i]))
            }
        }

        target.FDR <- target.fdr
        targets <- data.frame(target.FDR, F.critical, n.sig.genes)
         
        if(print.message.on.screen) cat("Done.\n")

        return(list(fdr=fdr, pi0=pi0, q.trim=q.trim, n.iter=n.iter, F.null=F.null, targets=targets))

}
