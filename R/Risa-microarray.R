getMicroarrayAssayFilenames <- function(isa){
  assay.filenames <- isa["assay.filenames"]
  assay.files <- isa["assay.files"]
  microarray.assay.filenames <- assay.filenames[ sapply(assay.files, function(x) isatab.syntax$hybridization.assay.name %in% names(x)) ]
  return(microarray.assay.filenames)
}

isMicroarrayAssay <- function(isa, assay.filename){
  microarray.assay.filenames <- getMicroarrayAssayFilenames(isa)
  return(assay.filename %in% microarray.assay.filenames)
}

getMIAMEMetadata <- function(isa, assay.filename){
  
  if (isMicroarrayAssay(isa, assay.filename)){
    
    i <- which(names(isa["assay.files"])==assay.filename)
        
    j <- which( lapply(isa["assay.filenames.per.study"], function(x)  (assay.filename %in% x)) == TRUE)
  
    my.desc <- new("MIAME", 
              name = isa@study.identifiers[[j]], 
              lab = isa@study.contacts.affiliations, 
              contact = isa@study.contacts, 
              title = isa@study.titles[[j]],
              abstract = isa@study.descriptions[[j]],
              samples = as.list(isa@samples))

  return(my.desc)
  }else{
    message("The assay is not a Microarray assay, so the MIAME metadata cannot be built")    
  }
}

getMicroarrayDerivedDataFilenames <- function(isa, full.path = TRUE){  
  microarray.assay.filenames <- getMicroarrayAssayFilenames(isa)
  microarray.files <- lapply(isa["data.filenames"][microarray.assay.filenames], function(x) x[isatab.syntax$derived.array.data.file])
  if (full.path)
    microarray.files <- sapply(microarray.files, function(x) sapply(x, function(y) paste(isa["path"], y, sep=.Platform$file.sep)))  
  return(microarray.files)
}

getMicroarrayDerivedDataFilenamesAssay <- function(isa, assay.filename, full.path = TRUE){  
  if (!isMicroarrayAssay(isa, assay.filename))
    stop("The ", assay.filename, " is not a microarray assay")

  microarray.files <-  isa["data.filenames"][[assay.filename]][isatab.syntax$derived.array.data.file]
  if (full.path)
    microarray.files <- sapply(microarray.files, function(x) sapply(x, function(y) paste(isa["path"], y, sep=.Platform$file.sep)))  
  return(microarray.files)
}


getExpressionSet <- function(isa, assay.filename){
  
  if (!isMicroarrayAssay(isa, assay.filename))
    stop("The ", assay.filename, " is not a microarray assay")

  suppressPackageStartupMessages(require("affy"))
  
  pd <- getAnnotatedDataFrameAssay(isa, assay.filename)
  
  miame <- getMIAMEMetadata(isa, assay.filename)
  
  i <- which(names(isa["assay.files"])==assay.filename)
  
  cel.files <- getAssayRawDataFilenames(isa@assay.tabs[[i]], full.path=FALSE)
  
  fnames <- as.vector(cel.files[[1]])
  
  current.path <- getwd()
  
  setwd(isa@path)
  
  eset <- justRMA(filenames=fnames, phenoData=pd, description=miame)
  
  setwd(current.path)  
  
  return(eset)
}