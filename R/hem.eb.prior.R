##################################################
#
#  Empirical Bayes (EB) Prior Specification
#  Updated on 2005.08.12
#
##################################################

hem.eb.prior <- function(dat,  n.layer, design, 
                method.var.e="neb", method.var.b="peb", method.var.t="neb", 
                rep=TRUE, baseline.var="LPE", p.remove=0, max.chip=4,
                q=0.01, B=25, n.digits=10, print.message.on.screen=TRUE)
{

     if(n.layer != 1 & n.layer != 2) stop("n.layer should be 1 or 2.")       
     if(print.message.on.screen) cat("EB...Please wait.\n")

     #Organize design
     ncol.design <- ncol(design)
     if(ncol.design > 3 | ncol.design < 2 ) stop(" incorrect size of design") 
     for (j in ncol(design)) {
        u <- unique(design[,j])
        for (i in 1:length(u)){
             w <- which(design[,j]==u[i])
             design[w,j] <- i 
        }
     }

     #design <- as.numeric(design)
     if(ncol(design)==3) i <- order(design[,1],design[,2],design[,3])
     else  i <- order(design[,1],design[,2])
     dat <- dat[,i]
     design <- design[i,]
     if(!is.data.frame(dat)) dat <- data.frame(dat)


     if(nrow(dat)*q <10) stop("q is too small")
     if(q > 0.5) stop("q is too large")

     #parameters 
     n.gene <- nrow(dat)
     n.chip <- ncol(dat)
     cond   <- design[,1]
     u.cond <- unique(cond)
     n.cond <- length(u.cond)
     c.size <- table(design[,1])
     n.bin <- length(seq(0,1,q)) - 1


     ########################################################################################
     #EB for two-layer HEM
     if(n.layer==2) 
     {

          method.vare <- ifelse(method.var.e=="peb", 2, ifelse(method.var.e=="neb", 3, 0))
          method.varb <- ifelse(method.var.b=="peb", 2, 0)

          if(method.vare < 2 | method.vare > 3) stop("method.var.e = peb or neb")       
          if(method.varb != 2) stop("method.var.b = peb")       

          #Parametric EB priors for experimental, biological, total error var
          med.expr <- matrix(NA,n.gene,n.cond)
          var.bio  <- matrix(NA,n.gene,n.cond)
          var.exp  <- matrix(NA,n.gene,n.cond)
          var.tot  <- matrix(NA,n.gene,n.cond)
          if(method.vare==2) var.e <- matrix(NA,1,n.bin)

          if(print.message.on.screen) cat("  conditions ")
          for(j in 1:n.cond){

              if(print.message.on.screen) cat(" ", j)
              w <- which(u.cond[j]==cond)
              x <- dat[,w] 
              r <- design[w,2]
              var.par <- par.rep.error.Olig(x, q=q, r=r)

              med.expr[,j] <- round(var.par$m, n.digits)
              var.bio[,j]  <- round(var.par$var.bio, n.digits)
              var.exp[,j]  <- round(var.par$var.exp, n.digits)
              var.tot[,j]  <- round(var.par$var.total, n.digits)
              if(method.vare==2) var.e <- rbind(var.e, round(var.par$var.e,n.digits))

          }
          if(method.vare==2) var.e <- var.e[-1,]
          var.bio[which(var.bio <=0)] <- min(var.bio[which(var.bio >0)])/2 ###NOTE
          if(print.message.on.screen) cat(" \n")


          if(method.vare==3) {
             #Nonparametric EB prior for experimental error var
             var.e <- matrix(NA,1,n.bin)
             if(print.message.on.screen) cat("  conditions ")
             for(j in 1:n.cond){
                 if(print.message.on.screen) cat(" ", j)
                 w <- which(u.cond[j]==design[,1])
                 x <- dat[,w] 
                 r <- design[w,2]
                 var.nonpar <- nonpar.rep.error.Olig(x, q=q, r=r, B=B, print.message.on.screen=FALSE) 
                 var.e <- rbind(var.e,round(var.nonpar,  n.digits))
             }
             var.e <- var.e[-1,]
          }

          if(print.message.on.screen) cat(" Done.\n")
          #return(var.b=round(var.bio,n.digits), var.e=round(var.exp,n.digits), q=q, B=B) 
          return(list(med.expr = med.expr, 
                 var.bio  = var.bio, var.exp  = var.exp, var.tot= var.tot,
                 var.b    = var.bio, var.e    = var.e,   q=q, B=B)) 
    }


    ########################################################################################
    #EB for one-layer HEM with replication
    if((n.layer==1) & rep) 
    {

          method.vart <- ifelse(method.var.t=="peb", 2, ifelse(method.var.t=="neb", 3, 0))
          if(method.vart < 2 | method.vart > 3) stop("method.var.t = peb or neb")
       
          #Parametric EB priors for total error var
          med.expr <- matrix(NA,n.gene,n.cond)
          var.tot  <- matrix(NA,n.gene,n.cond)


          if(print.message.on.screen) cat("  conditions ")
          for(j in 1:n.cond){
              if(print.message.on.screen) cat(" ", j)
              w <- which(u.cond[j]==cond)
              x <- dat[,w] 
              var.par <- par.error.Olig(x, q=q) 
              med.expr[,j] <- round(var.par$m, n.digits)
              var.tot[,j] <-  round(var.par$var.total, n.digits)                 
          }
          if(method.vart==2) var.t <- var.tot
          if(print.message.on.screen) cat(" \n")
 
          #Nonparametric EB prior for experimental error var
          if(method.vart==3) {
             var.t <- matrix(NA,1,n.bin) 
             if(print.message.on.screen) cat("  conditions ")
             for(j in 1:n.cond){
                 if(print.message.on.screen) cat(" ",j)
                 w <- which(u.cond[j]==design[,1])
                 x <- dat[,w] 
                 var.nonpar <- nonpar.error.Olig(x, q=q, B=B, print.message.on.screen=FALSE) 
                 var.t <- rbind(var.t,var.nonpar)
             }
             var.t <- round(var.t[-1,], n.digits)
          }

          if(print.message.on.screen) cat(" Done.\n")
          return(list(med.expr = med.expr, var.tot=var.tot, var.t=var.t, q=q, B=B)) 

    }


    ########################################################################################
    #EB for one-layer HEM without replication
    if((n.layer==1) & !rep) 
    {

          method.vart <- ifelse(method.var.t=="peb", 2, ifelse(method.var.t=="neb", 3, 0))
          if(method.vart < 2 | method.vart > 3) stop("method.var.t = peb or neb")

          #Parametric EB priors
          med.expr <- round(dat, n.digits)
          var.par <- par.no.error.Olig(dat, q=q, baseline.var=baseline.var, p.remove=p.remove) 
          var.tot <-  round(var.par$var.total, n.digits)                 
          if(method.vart==2) var.t <- var.tot
 
          #Nonparametric EB prior
          if(method.vart==3) {
             var.t <- matrix(NA,1,n.bin) 
             var.nonpar <- nonpar.no.error.Olig(dat, q=q, B=B, 
                           baseline.var=baseline.var, p.remove=p.remove, 
                           print.message.on.screen=FALSE) 
             for(j in 1:n.cond) var.t <- rbind(var.t,var.nonpar)
             var.t <- round(var.t[-1,], n.digits)
          }

          if(print.message.on.screen) cat(" Done.\n")
          return(list(med.expr = med.expr, var.tot=var.tot, var.t=var.t, q=q, B=B, 
                      rep=rep, baseline.var=baseline.var, p.remove=p.remove))  

    }


}

######END############################################################################


