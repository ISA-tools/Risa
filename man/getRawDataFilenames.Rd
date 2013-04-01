\name{getRawDataFilenames}
\alias{getRawDataFilenames}
\title{
Retrieves the raw data filenames in the ISATab dataset
}
\description{
Retrieves the raw data filenames in the ISATab dataset
}
\usage{
getRawDataFilenames(.Object, full.path=TRUE)
}
\arguments{
  \item{.Object}{ an \code{ISATab-class} object }
\item{full.path}{a \code{boolean} value indicating if the result should include the full path of the raw data filenames or not. The default value is TRUE.
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
