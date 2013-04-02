\name{getAssayRawDataFilenames-methods}
\docType{methods}
\alias{getAssayRawDataFilenames}
\alias{getAssayRawDataFilenames-methods}
\alias{getAssayRawDataFilenames,AssayTab,logical-method}
\alias{getAssayRawDataFilenames,MicroarrayAssayTab,logical-method}
\alias{getAssayRawDataFilenames,MSAssayTab,logical-method}
\alias{getAssayRawDataFilenames,NMRAssayTab,logical-method}
\alias{getAssayRawDataFilenames,SeqAssayTab,logical-method}
\title{ ~~ Methods for Function \code{getAssayRawDataFilenames} in Package \pkg{Risa} ~~}
\description{
 ~~ Methods for function \code{getAssayRawDataFilenames} in package \pkg{Risa} ~~
They retrieve the raw data files for a particular assay type.
}
\section{Methods}{
\describe{

\item{\code{signature(.Object = "AssayTab", full.path = "logical")}}{
Retrieves the raw data filenames for a generic assay.
}

\item{\code{signature(.Object = "MicroarrayAssayTab", full.path = "logical")}}{
Retrieves the raw data filenames for an assay whose technology type is DNA microarray.
}

\item{\code{signature(.Object = "MSAssayTab", full.path = "logical")}}{
Retrieves the raw data filenames for an assay whose technology type is mass spectrometry. These data filenames correspond to those specified in the column 'Raw Spectral Data File'.
}

\item{\code{signature(.Object = "NMRAssayTab", full.path = "logical")}}{
Retrieves the raw data filenames for an assay whose technology type is NMR spectroscopy.
}

\item{\code{signature(.Object = "SeqAssayTab", full.path = "logical")}}{
Retrieves the raw data filenames for an assay whose technology type is nucleotide sequencing.
}
}}
\keyword{methods}
\keyword{raw data filenames}
