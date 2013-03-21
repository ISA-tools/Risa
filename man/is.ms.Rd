\name{is.ms}
\alias{is.ms}
\title{
Indicates if an assay filename corresponds to a mass spectrometry assay.
}
\description{
Indicates if an assay filename corresponds to a mass spectrometry assay.
}
\usage{
is.ms(isa, assay.filename)
}
\arguments{
  \item{isa}{
an object from the ISAtab-class.
}
  \item{assay.filename}{
the filename of an assay file
}
}
%\details{
%}
\value{
It retrieves a boolean value indicating wether the assay is a mass spectrometry assay or not.
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
faahkoISA = readISAtab(find.package("faahKO"))
assay.filename <- faahkoISA["assay.filenames"][[1]]
is.ms(faahkoISA, assay.filename)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line