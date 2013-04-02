\name{getAnnotatedDataFrameAssay}
\alias{getAnnotatedDataFrameAssay}
\title{
Retrieves an AnnotatedDataFrameAssay.
}
\description{
Retrieves an AnnotatedDataFrameAssay.
}
\usage{
getAnnotatedDataFrameAssay(isa, assay.filename)
}
\arguments{
  \item{isa}{
An object of the \code{ISAtab-class}.
}
  \item{assay.filename}{
%A \code{character} .
}
}
%\details{
%}
\value{
An annotated data frame.
}
%\references{
%}
\author{
Alejandra Gonzalez-Beltran
}
%\note{
%}
\seealso{
\code{ISATab-class}, \code{AnnotatedDataFrame}
}
\examples{
 faahkoISAzip <- system.file("extdata","faahKO-metadata.zip", package="Risa")
 faahkoISA <- readISAtab(zipfile = faahkoISAzip, path = file.path(system.file("extdata", package="Risa"),"faahkoISA"), verbose =TRUE)
 faahkoADF <- getAnnotatedDataFrameAssay(faahkoISA, faahkoISA@assay.filenames[[1]])
 str(faahkoADF)
 
 
 armstrongISAzip <- system.file("extdata","ARMSTRONG-S-3-metadata.zip", package="Risa")
 armstrongISA <- readISAtab(zipfile = armstrongISAzip, path = file.path(system.file("extdata", package="Risa"),"armstrongISA"), verbose =TRUE)
 armstrongADF <- getAnnotatedDataFrameAssay(armstrongISA, armstrongISA@assay.filenames[[1]])
 str(armstrongADF)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
