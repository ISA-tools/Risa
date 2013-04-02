\name{write.ISAtab}
\alias{write.ISAtab}
\docType{data}
\title{ Write an ISA-Tab dataset }
\description{
  \code{write.ISAtab} writes a ISA-Tab dataset.
}
\usage{
write.ISAtab(isa, path = getwd())
}
\arguments{
  \item{isa}{ the ISA-Tab object as built with the function \code{readISAtab}.}
  \item{path}{ the directory where the files are going to be written, by default is the working directory}
  }
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com> 
}
\seealso{ \code{write.investigation.file}, \code{write.study.file}, \code{write.assay.file}, \code{readISAtab} }
\examples{
    ### Example of writing the study file of faahKO ISA-Tab dataset into the temp director
    data.dir <- system.file("extdata", package="Risa")
    isazip <- system.file("extdata","faahKO-metadata.zip", package="Risa")
    faahkoISA <- readISAtab(zipfile = isazip, path = file.path(data.dir,"faahKOISA"), verbose =TRUE)      
    temp = tempdir()
    write.ISAtab(faahkoISA, temp)
}