\name{getStudyFilename}
\alias{getStudyFilename}
\title{
Retrieves the study filename to which an assay.filename corresponds.
}
\description{
Retrieves the study filename to which an assay.filename corresponds.
}
\usage{
getStudyFilename(isa, assay.filename)
}
\arguments{
  \item{isa}{
an object from the \code{ISATab-class}.
}
  \item{assay.filename}{
a string with the filename of an assay
}
}
%\details{
%}
\value{
a string with the study filename.
}
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
faahkoISA <- readISAtab(find.package("faahKO"))
getStudyFilename(faahkoISA, faahkoISA@assay.filenames[[1]])
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
