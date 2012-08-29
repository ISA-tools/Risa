## getters
#' extract parts of ISAtab class
#'
#' @name [
#' @aliases [,ISAtab-method
#' @docType methods
#' @rdname extract-methods
setMethod(f="[",signature="ISAtab", definition=function(x, i,j, drop) {
  if (i=="path") { return(x@path) } else {}
  if (i=="investigation.filename") { return(x@investigation.filename) } else {}
  if (i=="investigation.file") { return(x@investigation.file) } else {}
  if (i=="study.identifiers") { return(x@study.identifiers) } else {}
  if (i=="study.filenames") { return(x@study.filenames) } else {}
  if (i=="study.files") { return(x@study.files) } else {}   
  if (i=="assay.filenames") { return(x@assay.filenames) } else {}
  if (i=="assay.filenames.per.study") { return(x@assay.filenames.per.study) } else {}
  if (i=="assay.files") { return(x@assay.files) } else {}
  if (i=="assay.files.per.study") { return(x@assay.files.per.study) } else {}
  if (i=="assay.technology.types") { return(x@assay.technology.types) } else {}
  if (i=="assay.measurement.types") { return(x@assay.measurement.types) } else {}
  if (i=="data.filenames") { return(x@data.filenames) } else {}
  if (i=="samples") { return(x@samples) } else {}
  if (i=="samples.per.study") { return(x@samples.per.study) } else {}
  if (i=="samples.per.assay.filename") { return(x@samples.per.assay.filename) } else {}
  if (i=="assay.filenames.per.sample") { return(x@assay.filenames.per.sample) } else {}
  if (i=="sample.to.rawdatafile") { return(x@sample.to.rawdatafile) } else {}
  if (i=="sample.to.assayname") { return(x@sample.to.assayname) } else {}
  if (i=="rawdatafile.to.sample") { return(x@rawdatafile.to.sample) } else {}
  if (i=="assayname.to.sample") { return(x@assayname.to.sample) } else {}
}
) 
