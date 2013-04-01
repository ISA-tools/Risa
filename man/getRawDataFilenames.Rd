\name{getRawDataFilenames}
\alias{getRawDataFilenames}
\title{
Retrieves the raw data filenames in the ISATab dataset
}
\description{
Retrieves the raw data filenames in the ISATab dataset
}
\usage{
getRawDataFilenames(x)
}
\arguments{
  \item{x}{
an ISATab object
}
}
%\details{
%}
%\value{
%}
%\references{
%}
\author{
Alejandra Gonzalez-Beltran
}
%\note{
%}
\seealso{
\code{ISATab-class}
}
\examples{
data.dir <- system.file("extdata", package="Risa")
isa <- readISAtab(file.path(data.dir, "ARMSTRONG-S-3"))
getRawDataFilenames(isa, full.path=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
