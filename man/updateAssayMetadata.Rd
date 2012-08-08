\name{updateAssayMetadata}
\alias{updateAssayMetadata}
\docType{data}
\title{Update metadata into the assay file}
\description{
  \code{updateAssayMetadata} Updates metadata into a particular assay file.
}
\usage{
updateAssayMetadata(isa, assay.filename, col.name, values)
}
\arguments{
  \item{isa}{ An isatab object, as retrieved by the \code{isatab2bioc} or the \code{isatab2bioczip} functions.}
  \item{assay.filename}{ the filename of the assay file to be augmented/modified}
  \item{col.name}{ the name of the column of the assay file to be modified}
  \item{values}{ the values to be added to the column of the assay file: it could be a single value, and in this case the value is repeated across the column, or it could be a list of values (whose length must match the number of rows of the assay file) }
  }
\value{
  The updated ISA-Tab object.
}
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com> 
}
\seealso{ \code{isatab2bioc} , \code{isatab2bioczip}}