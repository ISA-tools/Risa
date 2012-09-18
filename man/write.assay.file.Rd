\name{write.assay.file}
\alias{write.assay.file}
\docType{data}
\title{ Write an Assay File from the ISA-Tab dataset }
\description{
  \code{write.assay.file} writes a specific assay file from the ISA-tab dataset into disk.
}
\usage{
write.assay.file(isa, assay.filename, path = getwd())
}
\arguments{
  \item{isa}{ the ISA-Tab object (as retrieved by the \link{readISAtab} function )}
  \item{assay.filename}{ the name of the assay file to be written to disk}
  \item{path}{ the output path in which the study file is going to be written, by default is the working directory}
  }
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com>  
}
\seealso{ \code{write.ISAtab},\code{write.investigation.file}, \code{write.study.file}, \code{readISAtab} }
\examples{
    ### Example of writing the study file of faahKO ISA-Tab dataset into the temp directory
    faahkoISA <- readISAtab(find.package("faahKO"))
    temp = tempdir()
    write.assay.file(faahkoISA, faahkoISA["assay.filenames"][[1]], temp)
}