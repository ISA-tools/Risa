\name{isatab2bioc}
\alias{isatab2bioc}
\docType{data}
\title{ Build R objects from ISA-tab format files }
\description{
  \code{isatab2bioc} reads a ISA-tab file.
}
\usage{
isatab2bioc(path = getwd(),verbose=FALSE)
}
\arguments{
  \item{path}{ the name of the directory containing ISAtab files. The
    default is the working directory.}
  \item{verbose}{ a boolean indicating to show messages for the different steps, if TRUE, or not to show them, if FALSE}
  }
\value{
  The output is a list composed of R objects corresponding to elements of the ISA-Tab dataset, such as:
    path - the path of the ISA-Tab dataset, 
    investigation.filename - the name of the Investigation file
    investigation.file - a data frame with the contents of the Investigation file
    study.identifiers - the list of study identifiers
    study.filenames - the names of the study files
    study.files - a list of data frames wiht the contents of the study files
    assay.filenames - the names of the assay files
    assay.filenames.per.study - the names of the assay files according to the study they belong to
    assay.files - a list of data frames with the contents of the assay files
    assay.files.per.study - a list of data frames with the contents of the assay files divided per study they belong to
    assay.technology.types - a list with the technology types corresponding to each assay
    data.filenames - a list with the names of the data files
    samples - a list with the names of the samples
    samples.per.assay.filename - the samples classified according to the assay filename they belong to
    assay.filenames.per.sample - the names of the assay files classified per sample name
    sample.to.rawdatafile - a list of data frames with samples associated with raw data files
    sample.to.assayname
    rawdatafile.to.sample
    assayname.to.sample
}
\author{
  Alejandra Gonzalez-Beltran, Audrey Kauffmann / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com>
}
\seealso{ \code{isatab2bioczip} }

\examples{
     ## Example for mass spectrometry dataset
     temp = tempdir()
     datafiles = c(file.path(system.file("cdf/KO", package="faahKO"), grep("CDF",dir(system.file("cdf/KO", package="faahKO")),ignore.case=TRUE,value=TRUE)), file.path(system.file("cdf/WT", package="faahKO"), grep("CDF",dir(system.file("cdf/WT", package="faahKO")),value=TRUE)))
     file.copy(datafiles, temp, recursive=TRUE)
     isafiles = file.path(system.file(package="faahKO"), grep("txt",dir(system.file(package="faahKO")),value=TRUE))
     file.copy(isafiles, temp, recursive=TRUE)
     isaset = isatab2bioc(path = temp)
}
