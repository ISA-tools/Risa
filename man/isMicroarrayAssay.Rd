\name{isMicroarrayAssay}
\alias{isMicroarrayAssay}
\title{
Determines if an assay is microarray-based or not.
}
\description{
Retrieves a boolean value indicating wether the assay filename given as parameter corresponds to a microarray-based assay or not.
}
\usage{
isMicroarrayAssay(isa, assay.filename)
}
\arguments{
  \item{isa}{
An object of the \code{ISAtab-class}.
}
  \item{assay.filename}{
A string indicating the name of an assay file. It must correspond to a microarray-based assay for the function to return a valid value.
}
}
%\details{
%}
\value{
Retrieves a boolean value.
}
%\references{
%}
\author{
Alejandra Gonzalez-Beltran
}
%\note{
%}
\seealso{
\code{getMIAMEMetadata}
}
\examples{
data.dir <- system.file("extdata", package="Risa")
isazip <- system.file("extdata","ARMSTRONG-S-3-metadata.zip", package="Risa")
isa <- readISAtab(zipfile = isazip, path = file.path(data.dir,"ARMSTRONG_ISA"), verbose =TRUE)
    
assay.filename <- isa@assay.filenames[[1]]
isMicroarrayAssay(isa, assay.filename)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
