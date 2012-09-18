\name{ISAtab-class}
\Rdversion{1.1}
\docType{class}
\alias{ISAtab-class}
\alias{[,ISAtab-method}
\alias{[<-,ISAtab-method}
\title{Class \code{"ISAtab"}, a class for maintaing ISAtab information}
\description{
  This class maintains informatino from an ISAtab dataset, including an investigation file, one or more study files, and one or more assay files for each study file. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ISAtab", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{path}:}{Object of class \code{"character"}, it represents the path of the ISA-Tab dataset }
    \item{\code{investigation.filename}:}{Object of class \code{"character"}, it represents the name of the Investigation file }
    \item{\code{investigation.file}:}{Object of class \code{"data.frame"}, it is a data frame with the contents of the Investigation file }
    \item{\code{study.identifiers}:}{Object of class \code{"factor"}, it is the list of study identifiers }
    \item{\code{study.filenames}:}{Object of class \code{"character"}, the names of the study files }
    \item{\code{study.files}:}{Object of class \code{"list"}, a list of data frames wiht the contents of the study files }
    \item{\code{assay.filenames}:}{Object of class \code{"character"}, the names of the assay files }
    \item{\code{assay.filenames.per.study}:}{Object of class \code{"list"}, the names of the assay files according to the study they belong to }
    \item{\code{assay.files}:}{Object of class \code{"list"}, a list of data frames with the contents of the assay files }
    \item{\code{assay.files.per.study}:}{Object of class \code{"list"}, a list of data frames with the contents of the assay files divided per study they belong to }
    \item{\code{assay.technology.types}:}{Object of class \code{"character"}, a list with the technology types corresponding to each assay }
    \item{\code{assay.measurement.types}:}{Object of class \code{"character"}, a list with the names of the data files }
    \item{\code{data.filenames}:}{Object of class \code{"list"},    a list with the names of the samples }
    \item{\code{samples}:}{Object of class \code{"character"} }
    \item{\code{samples.per.study}:}{Object of class \code{"list"}, }
    \item{\code{samples.per.assay.filename}:}{Object of class \code{"list"},the samples classified according to the assay filename they belong to }
    \item{\code{assay.filenames.per.sample}:}{Object of class \code{"character"}, the names of the assay files classified per sample name }
    \item{\code{sample.to.rawdatafile}:}{Object of class \code{"list"}, a list of data frames with samples associated with raw data files}
    \item{\code{sample.to.assayname}:}{Object of class \code{"list"}, maintains the association between samples and assay names }
    \item{\code{rawdatafile.to.sample}:}{Object of class \code{"list"}, mantains the association raw data file and samples }
    \item{\code{assayname.to.sample}:}{Object of class \code{"list"}, mantains the associations between assay names and samples }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "ISAtab")}: ... }
   }
}
\references{
Risa is a component of the ISA software infrastructure \url{http://isa-tools.org/}. The source code and latest version can be found in the GitHub repository \url{https://github.com/ISA-tools/Risa}. Please, submit all 'bugs' and feature requests through \url{https://github.com/ISA-tools/Risa/issues}.
}
\author{
Alejandra Gonzalez-Beltran, \email{isatools@googlegroups.com}
}
\note{
No notes yet.
}
\seealso{
\code{ISAtab-methods}
}
\examples{
showClass("ISAtab")
}
\keyword{classes}