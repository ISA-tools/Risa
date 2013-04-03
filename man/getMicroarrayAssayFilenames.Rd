\name{getMicroarrayAssayFilenames}
\alias{getMicroarrayAssayFilenames}
\title{
Listing filenames of microarray-based asssay from an ISAtab dataset.
}
\description{
Retrieves the list of asssay filenames in an ISAtab dataset, which correspond to microarray-based assays.
}
\usage{
getMicroarrayAssayFilenames(isa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isa}{
An object of the \code{ISAtab-class}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a list of assay filenames
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
getMicroarrayAssayFilenames(armstrongISA)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
