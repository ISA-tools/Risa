\name{getMSAssayFilenames}
\alias{getMSAssayFilenames}
\title{
Retrieves a vector with the assay filenames corresponding to mass spectrometry.
}
\description{
Given an object from the ISAtab-class, this method retrieves a vector with the assay filenames corresponding to mass spectrometry (identified by having a column called 'Raw Spectral Data File'.
}
\usage{
getMSAssayFilenames(isa)
}
\arguments{
  \item{isa}{
an object from the ISAtab class
}
}
%\details{
%}
\value{
A character vector with the assay filenames corresponding to mass spectrometry.
}
%\references{
%}
\author{
Alejandra Gonzalez-Beltran
}
%\note{
%}
%\seealso{
%}
\examples{
  data.dir <- system.file("extdata", package="Risa")
  isazip <- system.file("extdata","faahKO-metadata.zip", package="Risa")
  faahkoISA <- readISAtab(zipfile = isazip, path = file.path(data.dir,"faahKOISA"), verbose =TRUE)  
  ms.assay.filenames <- getMSAssayFilenames(faahkoISA)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
