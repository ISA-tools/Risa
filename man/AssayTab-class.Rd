\name{AssayTab-class}
\Rdversion{1.1}
\docType{class}
\alias{AssayTab-class}
\alias{[,AssayTab-method}
\alias{[,AssayTab,ANY,ANY,ANY-method}
\alias{[<-,AssayTab-method}
\alias{[<-,AssayTab,ANY,ANY,ANY-method}
\alias{MSAssayTab-class}
\alias{MicroarrayAssayTab-class}
\alias{SeqAssayTab}
\alias{NMRAssayTab}
\title{Class \code{"AssayTab"}, a class for maintaing the information of AssayTab}
\description{
  This class maintains information from an AssayTab dataset, including an investigation file, one or more study files, and one or more assay files for each study file. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AssayTab", ...)}.
}
\section{Slots}{
  \describe{
   \item{\code{path}}{Object of class \code{"character"}, it represents the path of the ISA-Tab dataset to which this \code{AssayTab} belongs}
   \item{\code{study.filename}}{Object of class \code{"character"} with the filename of the study to which this assay belongs to.}
   \item{\code{study.identifier}}{Object of class \code{"character"} with the identifier of the study to which this assay belongs to.}
   \item{\code{assay.filename}}{Object of class \code{"character"} with the filename of the assay.}
   \item{\code{assay.file}}{Object of class \code{"data.frame"} with the contents of the assay file.}
   \item{\code{assay.technology.type}}{Object of class \code{"character"} with the technology type of this assay.}
   \item{\code{assay.measurement.type}}{Object of class \code{"character"} with the measurement type of this assay.}
   \item{\code{assay.names}}{Object of class \code{"data.frame"} with the assay names for this assay file.}
   \item{\code{data.filenames}}{Object of class \code{"data.frame"} with the data filenames.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "AssayTab")}: ... }
    \item{getAssayRawDataFilenames}{\code{signature(.Object = "AssayTab", full.path="logical")}: returns the raw data filenames for this assay}
   }
}
%\references{
%}
\author{
Alejandra Gonzalez-Beltran, \email{isatools@googlegroups.com}
}
%\note{
%No notes yet.
%}
\seealso{
\code{ISATab-class},\code{ISATab-methods}, 
}
\examples{
showClass("AssayTab")
showClass("MSAssayTab")
showClass("MicroarrayAssayTab")
showClass("SeqAssayTab")
showClass("NMRAssayTab")
}
\keyword{classes}