
##########################################################################
#
#        Heterogeneous Error Model (HEM) for Identification of 
#        Differentially Expressed Genes Under Multiple Conditions
#
#          developed by HyungJun Cho, PhD, and Jae K. Lee, PhD,
#              (hcho@virginia.edu; jaeklee@virginia.edu)           
#              Division of Biostatistics and Epidemiology
#              University of Virginia School of Medicine
#
##########################################################################

.First.lib <- function(lib, pkg) { 
   if(.Platform$OS.type=="windows" && interactive() && .Platform$GUI=="Rgui") 
      {addVigs2WinMenu("HEM") }
}


#####main function
hem <- function(dat, probe.ID=NULL, n.layer, design, burn.ins=1000, n.samples=3000,
                method.var.e="gam", method.var.b="gam", method.var.t="gam",           
                var.e=NULL, var.b=NULL, var.t=NULL, 
                var.g=1, var.c=1, var.r=1,
                alpha.e=3, beta.e=.1, alpha.b=3, beta.b=.1, alpha.t=3, beta.t=.2,
                n.digits=10, print.message.on.screen=TRUE)
{

        if(n.layer < 1 | n.layer >2) print("ERROR: n.layer = 1 or 2")       
        if(length(probe.ID) ==0) probe.ID <- 1:nrow(dat)                

        #Organize design
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

        #parameters 
        nMCMC  <- c(burn.ins,n.samples)
        n.Runs <- 1 + burn.ins + n.samples
        n.gene <- nrow(dat)
        n.chip <- ncol(dat)
        cond   <- design[,1]
        u.cond <- unique(design[,1])
        n.cond <- length(u.cond)
        c.size <- table(design[,1])
        ndat   <- c(n.gene, n.chip, n.cond)

     ########################################################################################
     #Two-layer HEM
     if(n.layer==2) 
     {
          par <-  c(var.g, var.c, var.r, alpha.e, beta.e, alpha.b, beta.b)
          grp <- rep(NA, n.cond) 
          nrp <- rep(NA,1)
          for(i in 1:n.cond) {
              w <- which(u.cond[i]==cond)
              grp[i] <-  length(unique(design[w,2]))
              nrp <- c(nrp, table(design[w,2]))
          }
          nrp <- nrp[-1]
          t.cond <- sum(grp)


          #Priors for variances
          method.vare <- ifelse(method.var.e=="gam", 1, ifelse(method.var.e=="peb", 2, 3))
          method.varb <- ifelse(method.var.b=="gam", 1, 2)

          opt <- c(method.vare, method.varb,0) 
          vare <- 1
          Bsize <-c(1,1)
          if(method.vare==2) {vare <- var.e}
          if(method.vare==3) {
              Brep <- round(nrow(var.e)/n.cond)-1
              quan <- ncol(var.e)
              Bsize <- c(Brep,quan)
              vare <- var.e
          }
 
          varb <- 1
          if(method.varb==2) varb <- var.b 

          #Fit HEM
          if(print.message.on.screen) cat("Modeling fitting...Please wait.\n")
          fit.hem <- .C("twolayerhem", 
                       as.double(t(dat)),  as.integer(opt),
                       as.integer(ndat), as.integer(grp), as.integer(nrp), 
                       as.integer(nMCMC), as.double(par),
                       as.integer(Bsize), as.double(t(vare)), as.double(t(varb)),
                       fstat = double(n.gene),
                       mexprest = matrix(double(1),n.cond, n.gene), 
                       mexpr    = matrix(double(1),t.cond, n.gene),
                       msigma2b = matrix(double(1),n.cond, n.gene),
                       msigma2e = matrix(double(1),t.cond, n.gene),
                       MCMCsamp = matrix(double(1), 7, n.Runs),
                       PACKAGE = "HEM"
                       )  

          m.mu    <- t(fit.hem$mexprest) 
          m.x     <- t(fit.hem$mexpr) 
          m.var.b <- t(fit.hem$msigma2b) 
          m.var.e <- t(fit.hem$msigma2e) 
          samples <- data.frame(t(fit.hem$MCMCsamp))
          names(samples) <- c("expr", "var.e", "mu", "gene", "cond", "inter", "var.b") 

          priors <- par
          H <- data.frame(matrix(fit.hem$fstat))
          rownames(H) <- probe.ID

          return(list(n.layer=n.layer, method.var.e=method.var.e, method.var.b=method.var.b, 
                 n.gene = n.gene, n.chip = n.chip, n.cond = n.cond, design=design, 
                 burn.ins = burn.ins, n.samples = n.samples, priors = priors, 
                 m.mu = round(m.mu,n.digits), 
                 m.x = round(m.x,n.digits), 
                 m.var.b = round(m.var.b,n.digits), 
                 m.var.e = round(m.var.e,n.digits), 
                 H = round(H,n.digits), 
                 samples = round(samples, n.digits))
            )
    }

    ########################################################################################
    #One-layer HEM
    if(n.layer==1) 
    {

          method.vart <- ifelse(method.var.t=="gam", 1, ifelse(method.var.t=="peb", 2, 3))

          par <-  c(var.g, var.c, var.r, alpha.t, beta.t)
          opt <- c(0, 0, method.vart) 
          grp <- c.size

          vart <- 1
          Bsize <-c(1,1)
          if(method.vart==2) {
             vart <- var.t
          }
          if(method.vart==3) {
              vart<- var.t
              Brep <- round(nrow(vart)/n.cond)-1
              quan <- ncol(vart)
              Bsize <- c(Brep,quan)
          }

         if(print.message.on.screen) cat("Modeling fitting...Please wait.\n")
         fit.hem <- .C("onelayerhem", 
                       as.double(t(dat)), as.integer(opt),
                       as.integer(ndat), as.integer(grp), 
                       as.integer(nMCMC), as.double(par), 
                       as.integer(Bsize), as.double(t(vart)),
                       fstat = double(n.gene),
                       mexprest = matrix(double(1),n.cond,n.gene), 
                       msigma2 = matrix(double(1),n.cond,n.gene),
                       MCMCsamp = matrix(double(1), 5, n.Runs),
                       PACKAGE = "HEM"
                       )  

          m.mu <- t(fit.hem$mexprest) 
          m.var.t <- t(fit.hem$msigma2) 
          priors <- par
          samples <- data.frame(t(fit.hem$MCMCsamp))
          names(samples) <- c("mu", "gene", "cond", "inter", "var.t") 

          H <- fit.hem$fstat
          H <- data.frame(matrix(fit.hem$fstat))
          rownames(H) <- probe.ID

          if(print.message.on.screen) cat(" Done.\n")
 
          return(list(n.layer=n.layer,  method.var.t=method.var.t, 
                 n.gene = n.gene, n.chip = n.chip, n.cond = n.cond, design=design, 
                 burn.ins = burn.ins, n.samples = n.samples, priors = priors, 
                 m.mu = round(m.mu,n.digits), 
                 m.var.t = round(m.var.t,n.digits), 
                 H=round(H,n.digits), 
                 samples = round(samples, n.digits)))
    }
}

######END############################################################################
