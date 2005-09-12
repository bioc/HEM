\name{hem.fdr}
\alias{hem.fdr}
\title{FDR Evaluation}
\description{
Computes resampling-based False Discovery Rate (FDR) 
}
\usage{
hem.fdr(dat,  n.layer, design, rep=TRUE, hem.out, eb.out=NULL, n.iter=5, q.trim=0.9,
        target.fdr=c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.20, 0.30, 0.40, 0.50),
        n.digits=10, print.message.on.screen=TRUE)
}
\arguments{
\item{dat}{data}
\item{n.layer}{number of layers: 1=one-layer EM; 2=two-layer EM}
\item{design}{design matrix}
\item{rep}{no replication if FALSE}  
\item{hem.out}{output from hem function}
\item{eb.out}{output from hem.eb.prior function}
\item{n.iter}{number of iterations}
\item{q.trim}{quantile used for estimtaing the proportion of true negatives (pi0)}
\item{target.fdr}{Target FDRs}
\item{n.digits}{number of digits}
\item{print.message.on.screen}{if TRUE, process status is shown on screen.}
}

\value{
\item{fdr}{H-values and corresponding FDRs}
\item{pi0}{estimated proportion of true negatives}
\item{H.null}{H-scores from null data}
\item{targets}{given target FDRs, corrsponding critical values and numbers of significant genes are provided}
}

\author{
 HyungJun Cho and Jae K. Lee
}

\seealso{
\code{\link{hem.eb.prior}} 
\code{\link{hem}} 
}


\examples{

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
#pbrain.eb  <- hem.eb.prior(pbrain.nor, n.layer=2,  design=design,
#                           method.var.e="neb", method.var.b="peb")     

##fit HEM with two layers of error
##using the small numbers of burn-ins and MCMC samples for a testing purpose;
##but increase the numbers for a practical purpose 
#pbrain.hem <- hem(pbrain.nor, n.layer=2,  design=design,burn.ins=10, n.samples=30, 
#              method.var.e="neb", method.var.b="peb", 
#              var.e=pbrain.eb$var.e, var.b=pbrain.eb$var.b)

##Estimate FDR based on resampling 
#pbrain.fdr <- hem.fdr(pbrain.nor, n.layer=2,  design=design, 
#                  hem.out=pbrain.hem, eb.out=pbrain.eb)

}

\keyword{models}
