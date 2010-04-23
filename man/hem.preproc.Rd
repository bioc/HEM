\name{hem.preproc}
\alias{hem.preproc}
\title{
Preprocessing
}
\description{
Performs IQR normalization, thesholding, and log2-transformation 
}
\usage{
hem.preproc(x, data.type = "MAS5")
}

\arguments{
\item{x}{data}
\item{data.type}{data type: MAS5 or MAS4}
}

\author{
HyungJun Cho and Jae K. Lee
}
\seealso{
\code{\link{hem}}, 
\code{\link{hem.eb.prior}}, 
\code{\link{hem.fdr}}
}



\examples{

data(pbrain)
pbrain.nor <- hem.preproc(pbrain[,2:13])

}

\keyword{models}
