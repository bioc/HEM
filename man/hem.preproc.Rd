\name{hem.preproc}
\alias{hem.preproc}
\title{
Preprocessing
}
\description{
Performs IQR normalization, thesholding, and log2-transformation 
}
\usage{
function (x, cond, data.type = "MAS5")
}

\arguments{
\item{x}{data}
\item{cond}{condion vector}
}

\author{
HyungJun Cho and Jae K. Lee
}
\seealso{
\code{\link{hem}}, 
\code{\link{hem.eb.prior}}, 
\code{\link{hem.fdr}}, 
\code{\link{lpe}} 
}



\examples{

library(HEM)

data(pbrain)
cond <- c(1,1,1,1,1,1,2,2,2,2,2,2)
pbrain.nor <- hem.preproc(pbrain[,2:13], cond=cond)

}

\keyword{models}
