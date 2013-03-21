\name{getRawDataFiles}
\alias{getRawDataFiles}
\title{
Retrieve raw data files for mass spectrometry assays.
}
\description{
The method \code{getRawDataFiles} retrieves raw data files for mass spectrometry assays, i.e. the files reported in the column 'Raw Spectral Data File' from the assay files.
}
\usage{
getRawDataFiles(isa, full.path=TRUE)
}
\arguments{
  \item{isa}{
  An object from the ISAtab class. 
}
  \item{full.path}{
  A boolean value to indicate if the data files in the result should be given with their full path (if TRUE, which is the default value) or just as they appear in the ISAtab dataset (if FALSE).
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
It retrieves a matrix, where each column contains the 'Raw Spectral Data Files' for each assay in the ISAtab archive. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Alejandra Gonzalez-Beltran
}
\note{
It will only work for mass spectrometry assays (containing the 'Raw Spectral Data File' column)
}
%\seealso{
%}
\examples{
faahkoISA <- readISAtab(find.package("faahKO"))
getRawDataFiles(faahkoISA)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
