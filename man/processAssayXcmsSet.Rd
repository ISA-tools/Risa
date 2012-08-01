\name{processAssayXcmsSet}
\alias{processAssayXcmsSet}
\docType{data}
\title{ Build an xcmsSet object given a mass spectrometry assay }
\description{
  \code{processAssayXcmsSet} reads a ISA-tab file.
}
\usage{
processAssayXcmsSet(isa, assay.filename, ...)
}
\arguments{
  \item{isa}{ an ISA object, as retrieved by the functions \code{isatab2bioc} and \code{isatab2bioczip}}
  \item{assay.filename}{ a boolean indicating to show messages for the different steps, if TRUE, or not to show them, if FALSE}
  \item{...}{ extra arguments that can be pass down to the xcmsSet function from the xcms package}
  }
\value{
  the xcmsSet object built from the assay metadata and data files
}
\author{
  Steffen Neumann, Alejandra Gonzalez-Beltran / 
  Maintainer: Audrey Kauffmann, Gabriella Rustici, Philippe Rocca-Serra, <isatools@googlegroups.com>  
}
\seealso{ \code{isatab2bioc}, \code{isatab2bioczip} }
