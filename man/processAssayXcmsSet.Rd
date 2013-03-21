\name{processAssayXcmsSet}
\alias{processAssayXcmsSet}
\docType{data}
\title{ Build an xcmsSet object given a mass spectrometry assay and considering all factors defined}
\description{
  \code{processAssayXcmsSet} retrieves an xcmsSet object given and ISA-tab object and one of its assay filenames.
}
\usage{
processAssayXcmsSet(isa, assay.filename, ...)
}
\arguments{
  \item{isa}{ an ISA object, as retrieved by the function \code{readISAtab}}
  \item{assay.filename}{ a boolean indicating to show messages for the different steps, if TRUE, or not to show them, if FALSE}
  \item{...}{ extra arguments that can be pass down to the xcmsSet function from the xcms package}
  }
\value{
  the xcmsSet object built from the assay metadata, considering all the factors defined, and data files
}
\author{
  Steffen Neumann, Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran  \email{isatools@googlegroups.com} 
}
\seealso{ \code{readISAtab} }
\examples{
 faahkoISA = readISAtab(find.package("faahKO"))
 assay.filename <- faahkoISA["assay.filenames"][[1]]
 xset = processAssayXcmsSet(faahkoISA, assay.filename)
}