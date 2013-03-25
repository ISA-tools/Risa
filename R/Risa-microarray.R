
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
