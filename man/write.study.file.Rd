\name{write.study.file}
\alias{write.study.file}
\docType{data}
\title{ Write a Study File of an ISA-Tab dataset }
\description{
  \code{write.study.file} writes a specific Study file from the ISA-tab dataset.
}
\usage{
write.study.file(isa, study.filename, path = getwd())
}
\arguments{
  \item{isa}{ the ISA-Tab object (as retrieved by the \link{readISAtab} function)}
  \item{study.filename}{ the name of the study file to be written to disk}
  \item{path}{ the output path in which the study file is going to be written, by default is the working directory}
  }
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com>  
}
\seealso{ \code{write.study.file}, \code{write.ISAtab}, \code{readISAtab} }
\examples{
    ### Example of writing the study file of faahKO ISA-Tab dataset into the temp director
    faahkoISA <- readISAtab(find.package("faahKO"))
    temp = tempdir()
    write.study.file(faahkoISA, faahkoISA["study.filenames"][[1]], temp)
}