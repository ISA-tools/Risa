\name{risazip}
\alias{risazip}
\docType{data}
\title{ Build R objects from ISA-tab format files }
\description{
  \code{risazip} reads a ISA-tab zip file.
}
\usage{
risazip(zip, path = getwd())
}
\arguments{
  \item{zip}{ a zip archive containing ISAtab files.}
  \item{path}{ the name of the directory in which the files
    from the zip archive will be extracted. The
    default is the working directory.}
  }
\value{
  The output is a list composed, from a data.frame containing the metadata information originally contained in the s_ file of ISAtab, followed by an object per assay of the study which will be an \code{\link[affy:AffyBatch-class]{AffyBatch}} if the assay is a raw Affymetrix microarray dataset, an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} for a raw one colour microarray dataset or for a processed microarray dataset, a \code{\link[Biobase:class.NChannelSet]{NChannelSet}} from a raw two colour dataset, a \code{\link[flowCore:flowSet-class]{flowSet}} for a flow cytometry dataset and a \code{\link[xcms:xcmsSet-class]{xcmsSet}} for a mass spectrometry dataset.
}
\author{
  Audrey Kauffmann / Maintainer: Audrey Kauffmann, Gabriella Rustici, Philippe Rocca-Serra, <isatools@googlegroups.com>  
}
\seealso{ \code{risatab} }
\examples{
     ## Example for microarray and flow cytometry experiment
     isazip = system.file("extdata","Yox1.zip", package="rISAtab")
     isaset = risazip(isazip, path = "ISAexample")
}
