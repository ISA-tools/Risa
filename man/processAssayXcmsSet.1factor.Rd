\name{processAssayXcmsSet.1factor}
\alias{processAssayXcmsSet.1factor}
\title{ 
Build an xcmsSet object given a mass spectrometry assay and considering the first factor defined 
}
\description{
\code{processAssayXcmsSet.1factor} retrieves an xcmsSet object given and ISA-tab object and one of its assay filenames.
}
\usage{
processAssayXcmsSet.1factor(isa, assay.filename, ...)
}
\arguments{
  \item{isa}{ an ISA object, as retrieved by the function \code{readISAtab}}
  \item{assay.filename}{ a boolean indicating to show messages for the different steps, if TRUE, or not to show them, if FALSE}
  \item{...}{ extra arguments that can be pass down to the xcmsSet function from the xcms package}
  }
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  the xcmsSet object built from the assay metadata, considering all the factors defined, and data files
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran  \email{isatools@googlegroups.com} 
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{readISAtab} 
}
\examples{
faahkoISA = readISAtab(find.package("faahKO"))
assay.filename <- faahkoISA["assay.filenames"][[1]]
xset = processAssayXcmsSet.1factor(faahkoISA, assay.filename)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
