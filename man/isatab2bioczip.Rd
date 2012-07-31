\name{isatab2bioczip}
\alias{isatab2bioczip}
\docType{data}
\title{ Build R objects from a zipped ISA-tab dataset}
\description{
  \code{isatab2bioczip} reads a ISA-tab zip file.
}
\usage{
isatab2bioczip(zip, path = getwd(),verbose=FALSE)
}
\arguments{
  \item{zip}{ a zip archive containing ISAtab files.}
  \item{path}{ the name of the directory in which the files
    from the zip archive will be extracted. The default is the working directory.}
  \item{verbose}{ a boolean indicating to show messages for the different steps, if TRUE, or not to show them, if FALSE}
  }
\value{
  The output is as in the function isatab2bioc.
}
\author{
  Alejandra Gonzalez-Beltran, Audrey Kauffmann / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com> 
}
\seealso{ \code{isatab2bioc} }
\examples{
     ## Example for microarray and flow cytometry experiment
     isazip = system.file("extdata","Yox1.zip", package="rISAtab")
     isaset = isatab2bioczip(isazip, path = "ISAexample")
}
