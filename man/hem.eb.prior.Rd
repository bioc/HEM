\name{hem.eb.prior}
\alias{hem.eb.prior}
\title{Empirical Bayes (EB) Prior Specification}
\description{
Estimates experimental and biological variances by LPE and resampling
}

\usage{
hem.eb.prior(dat,  n.layer, design, 
             method.var.e="neb", method.var.b="peb", method.var.t="neb", 
             rep=TRUE, baseline.var="LPE", p.remove=0, max.chip=4,
             q=0.01, B=25, n.digits=10, print.message.on.screen=TRUE)
}

\arguments{
\item{dat}{data}
\item{n.layer}{number of layers}
\item{design}{design matrix}
\item{method.var.e}{prior specification method for experimental variance;  
"peb"=parametric EB prior specification, "neb"=nonparametric EB prior specification}
\item{method.var.b}{prior specification method for biological variance; 
"peb"=parametric EB prior specification}
\item{method.var.t}{prior specification method for total variance; 
"peb"=parametric EB prior specification, "neb"=nonparametric EB prior specification}        
\item{rep}{no replication if FALSE}
\item{baseline.var}{baseline variance estimation method: LPE for replicated data and BLPE, PSE, or ASE for unreplicated data}
\item{p.remove}{percent of removed rank-variance genes for BLPE}
\item{max.chip}{maximum number of chips to estimate errors}
\item{q}{quantile for paritioning genes based on expression levels}
\item{B}{number of iterations for resampling}
\item{n.digits}{number of digits}
\item{print.message.on.screen}{if TRUE, process status is shown on screen.}
}

\value{
\item{var.b}{prior estimate matrix for biological variances (n.layer=2)}
\item{var.e}{prior estimate matrix for experiemtnal variances (n.layer=2)}
\item{var.t}{prior estimate matrix for total variances (n.layer=1)}
}

\author{
HyungJun Cho and Jae K. Lee
}

\seealso{
\code{\link{hem}}, 
\code{\link{hem.fdr}}
}


\examples{

#Example 1: Two-layer HEM with EB prior specification

data(pbrain)

##construct a design matrix
cond <- c(1,1,1,1,1,1,2,2,2,2,2,2)
ind  <- c(1,1,2,2,3,3,1,1,2,2,3,3)
rep  <- c(1,2,1,2,1,2,1,2,1,2,1,2)
design <- data.frame(cond,ind,rep)

##normalization
pbrain.nor <- hem.preproc(pbrain[,2:13])

##take a subset for a testing purpose;
##use all genes for a practical purpose
pbrain.nor <- pbrain.nor[1:1000,]

##estimate hyperparameters of variances by LPE
pbrain.eb  <- hem.eb.prior(pbrain.nor, n.layer=2,  design=design,
                           method.var.e="neb", method.var.b="peb")     

#fit HEM with two layers of error
#using the small numbers of burn-ins and MCMC samples for a testing purpose;
#but increase the numbers for a practical purpose 
pbrain.hem <- hem(pbrain.nor, n.layer=2,  design=design,burn.ins=10, n.samples=30, 
              method.var.e="neb", method.var.b="peb", 
              var.e=pbrain.eb$var.e, var.b=pbrain.eb$var.b)


#Example 2: One-layer HEM with EB prior specification

data(mubcp)

##construct a design matrix
cond <- c(rep(1,6),rep(2,5),rep(3,5),rep(4,5),rep(5,5))
ind  <- c(1:6,rep((1:5),4))
design <- data.frame(cond,ind)

##normalization
mubcp.nor <- hem.preproc(mubcp)

##take a subset for a testing purpose;
##use all genes for a practical purpose
mubcp.nor <- mubcp.nor[1:1000,] 

##estimate hyperparameters of variances by LPE
mubcp.eb  <- hem.eb.prior(mubcp.nor, n.layer=1, design=design,
             method.var.t="neb")                                

#fit HEM with two layers of error
#using the small numbers of burn-ins and MCMC samples for a testing purpose;
#but increase the numbers for a practical purpose 
mubcp.hem <- hem(mubcp.nor, n.layer=1, design=design,  burn.ins=10, n.samples=30, 
             method.var.t="neb", var.t=mubcp.eb$var.t)


###NOTE: see help(hem.fdr) for examples with FDR evaluation
###NOTE: see help(hem) for examples without EB

}



\keyword{models}

