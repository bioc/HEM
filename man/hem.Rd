\name{hem}
\alias{hem}
\title{Heterogeneous Error Model for Identification of Differential Expressed Genes Under Multiple Conditions}
\description{
Fits an error model with heterogeneous experimental and biological variances.
}
\usage{
hem(dat, probe.ID=NULL, n.layer, design, burn.ins=1000, n.samples=3000,
    method.var.e="gam", method.var.b="gam", method.var.t="gam",           
    var.e=NULL, var.b=NULL, var.t=NULL, var.g=1, var.c=1, var.r=1,
    alpha.e=3, beta.e=.1, alpha.b=3, beta.b=.1, alpha.t=3, beta.t=.2,
    n.digits=10, print.message.on.screen=TRUE)
}
\arguments{
\item{dat}{data}
\item{probe.ID}{a vector of probe set IDs}
\item{n.layer}{number of layers; 1=one-layer EM, 2=two-layer EM}
\item{design}{design matrix}
\item{burn.ins}{number of burn-ins for MCMC}
\item{n.samples}{number of samples for MCMC} 
\item{method.var.e}{prior specification method for experimental variance; 
"gam"=Gamma(alpha,beta), 
"peb"=parametric EB prior specification, 
"neb"=nonparametric EB prior specification}
\item{method.var.b}{prior specification method for biological variance; 
"gam"=Gamma(alpha,beta), 
"peb"=parametric EB prior specification}
\item{method.var.t}{prior specification method for total variance; 
"gam"=Gamma(alpha,beta), 
"peb"=parametric EB prior specification, 
"neb"=nonparametric EB prior specification}        
\item{var.e}{prior estimate matrix for experimental variance} 
\item{var.b}{prior estimate matrix for biological variance} 
\item{var.t}{prior estimate matrix for total variance}
\item{var.g}{N(0, var.g); prior parameter for gene effect} 
\item{var.c}{N(0, var.c); prior parameter for condition effect} 
\item{var.r}{N(0, var.r); prior parameter for interaction effect of gene and condition}
\item{alpha.e, beta.e}{Gamma(alpha.e,alpha.e); prior parameters for inverse of experimental variance} 
\item{alpha.b, beta.b}{Gamma(alpha.b,alpha.b); prior parameters for inverse of biological variance}
\item{alpha.t, beta.t}{Gamma(alpha.b,alpha.b); prior parameters for inverse of total variance}
\item{n.digits}{number of digits}
\item{print.message.on.screen}{if TRUE, process status is shown on screen.}
}

\value{
\item{n.gene}{numer of genes}
\item{n.chip}{number of chips}
\item{n.cond}{number of conditions} 
\item{design}{design matrix} 
\item{burn.ins}{number of burn-ins for MCMC} 
\item{n.samples}{number of samples for MCMC} 
\item{priors}{prior parameters} 
\item{m.mu}{estimated mean expression intensity for each gene under each condition}
\item{m.x}{estimated unobserved expression intensity for each combination of genes, conditions, and individuals (n.layer=2)}
\item{m.var.b}{estimated biological variances (n.layer=2)}
\item{m.var.e}{estimated experiemental variances (n.layer=2)}
\item{m.var.t}{estimated total variances (n.layer=1)}
\item{H}{H-scores} 
}

\references{
Cho, H. and Lee, J.K. (2004) Bayesian Hierarchical Error Model for Analysis of Gene Expression Data, 
Bioinformatics, 20: 2016-2025.
}

\author{
 HyungJun Cho and Jae K. Lee
}

\seealso{
\code{\link{hem.eb.prior}}, 
\code{\link{hem.fdr}}
}


\examples{

#Example 1: Two-layer HEM

data(pbrain)

##construct a design matrix
cond <- c(1,1,1,1,1,1,2,2,2,2,2,2) #condition
ind  <- c(1,1,2,2,3,3,1,1,2,2,3,3) #biological replicate
rep  <- c(1,2,1,2,1,2,1,2,1,2,1,2) #experimental replicate
design <- data.frame(cond,ind,rep)

##normalization
pbrain.nor <- hem.preproc(pbrain[,2:13])

##fit HEM with two layers of error
##using the small numbers of burn-ins and MCMC samples for a testing purpose;
##but increase the numbers for a practical purpose 
#pbrain.hem <- hem(pbrain.nor, n.layer=2, design=design, 
#                  burn.ins=10, n.samples=30)

##print H-scores
#pbrain.hem$H 


#Example 2: One-layer HEM

data(mubcp)

##construct a design matrix
cond <- c(rep(1,6),rep(2,5),rep(3,5),rep(4,5),rep(5,5))
ind  <- c(1:6,rep((1:5),4))
design <- data.frame(cond,ind)

##construct a design matrix
mubcp.nor <- hem.preproc(mubcp)

#fit HEM with one layers of error
#using the small numbers of burn-ins and MCMC samples for a testing purpose;
#but increase the numbers for a practical purpose 
#mubcp.hem <- hem(mubcp.nor, n.layer=1,design=design, burn.ins=10, n.samples=30)

##print H-scores
#mubcp.hem$H


###NOTE: Use 'hem.fdr' for FDR evaluation
###NOTE: Use 'hem.eb.prior' for Empirical Bayes (EB) prior sepecification
###NOTE: Use EB-HEM ('hem' after 'hem.eb.prior') for small data sets

}

\keyword{models}


