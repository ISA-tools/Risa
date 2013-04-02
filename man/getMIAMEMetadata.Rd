\name{getMIAMEMetadata}
\alias{getMIAMEMetadata}
\title{
Retrieves the MIAME information about a microarray experiment
}
\description{
If the assay.filename uses the microarray technology, it retrieves an object of the class \code{MIAME} with the minimum information about the experiment. Otherwise, it outputs a message indicating that the method is not valid for assays not using the microarray technology.
}
\usage{
getMIAMEMetadata(isa, assay.filename)
}
%- maybe also 'usage' for other objects documented here.
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
It retrieves an object from the \code{MIAME} class or an error message.
}
%\references{
%}
\author{
Alejandra Gonzalez-Beltran
}
%\note{
%}

\seealso{
\code{isMicroarrayAssay}
}
\examples{
 armstrongISAzip <- system.file("extdata","ARMSTRONG-S-3-metadata.zip", package="Risa")
 armstrongISA <- readISAtab(zipfile = armstrongISAzip, path = file.path(system.file("extdata", package="Risa"),"armstrongISA"), verbose =TRUE)
 assay.filename <- armstrongISA@assay.filenames[[1]]
 miame <- getMIAMEMetadata(armstrongISA, assay.filename)
 str(miame)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
