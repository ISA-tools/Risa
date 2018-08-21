Sys.setlocale('LC_ALL','C') 

isatab.syntax <- list(
  investigation.prefix="i_",
  study.prefix="s_",
  assay.prefix="a_",
  investigation.identifier="Investigation Identifier",
  study.identifier="Study Identifier",
  study.title="Study Title",
  study.description="Study Description",
  study.person.last.name="Study Person Last Name",
  study.person.first.name="Study Person First Name",
  study.person.mid.initial="Study Person Mid Initial",
  study.person.affiliation="Study Person Affiliation",
  study.file.name="Study File Name",
  study.assay.file.name="Study Assay File Name",
  study.assay.technology.type="Study Assay Technology Type",
  study.assay.measurement.type="Study Assay Measurement Type",
  sample.name="Sample Name",
  assay.name="Assay Name",
  data.file = "Data File",
  raw.data.file="Raw Data File",
  free.induction.decay.data.file="Free Induction Decay Data File",
  array.data.file="Array Data File",
  derived.array.data.file="Derived Array Data File",
  raw.spectral.data.file="Raw Spectral Data File",
  hybridization.assay.name="Hybridization Assay Name",
  factor.name="Factor Name",
  factor.value="Factor Value",
  assay.name="Assay Name"
  )

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

technology.types <- list(
  microarray="DNA microarray",
  ms="mass spectrometry",
  fc="flow cytometry",
  seq="nucleotide sequencing",
  nrm="NMR spectroscopy"
  )

readISAtab = function(path = getwd(), zipfile = NULL, verbose = FALSE)
{
  if (!is.null(zipfile)){
    readISAtabZip(zipfile, path, verbose)
  }
  else{
    readISAtabFiles(path, verbose)
  }
}

## This function only works if the zip file does not contain a directory (but the ISA-TAB files themselves)
readISAtabZip = function(zip, path = getwd(), verbose=FALSE)
{ 
  if (verbose)
    message("Unzipping file in directory ",path)
  d = unzip(zipfile = zip, exdir = extract <- path)  
  if (verbose)
    message("Unzipped files: ",d)
  isaobj = readISAtabFiles(path) 
  return(isaobj)
}##end function readISAtabZip

readISAtabFiles = function(path = getwd(), verbose=FALSE)
{
  if (verbose)
    message("Converting ISA-Tab dataset at ",path," into R objects...")
  isaobject <- new(Class="ISATab",path=path)
  if (verbose)
    message("... done.")
  return(isaobject) 
}##end function readISAtabFiles

updateAssayMetadata = function(isa, assay.filename, col.name, values){  
  assay.file <- isa["assay.files"][[ assay.filename ]]
  if (length(values)==1){
    values <- c(rep(values,nrow(assay.file)))
  }else if (length(values)!=nrow(assay.file)){
    stop("Wrong number of values to be added to the assay file")
  }
  ###update column of the assay.file
  assay.file [ colnames(assay.file) == col.name ] <- values
  #### update the isa object with modified assay.file
  isa <- setAssayFile(isa,assay.filename, assay.file)
  return(isa)
}

write.ISAtab = function(isa, path = getwd()){
  write.investigation.file(isa, path)
  for(i in seq_len(length(isa["study.filenames"]))){
    write.study.file(isa, isa["study.filenames"][[i]], path)
  }
  for(i in seq_len(length(isa["assay.filenames"]))){
    write.assay.file(isa, isa["assay.filenames"][[i]], path)
  }  
}

write.investigation.file = function(isa, path = getwd()){
  write.table(isa["investigation.file"], 
              file=file.path(path,isa["investigation.filename"]), 
              row.names=FALSE, col.names=FALSE, 
              quote=TRUE, sep="\t", na="\"\"")
}

write.study.file = function(isa, study.filename, path = getwd()){
  i <- which(isa["study.filenames"]==study.filename)
  study.file <- isa["study.files"][[ i ]]
  write.table(study.file, 
              file=file.path(path,isa["study.filenames"][[i]]), 
              row.names=FALSE, col.names=TRUE, 
              quote=TRUE, sep="\t", na="\"\"")
}

write.assay.file = function(isa, assay.filename, path = getwd()){
  i <- which(names(isa["assay.files"])==assay.filename)
  assay.file <- isa["assay.files"][[assay.filename ]]
  write.table(assay.file, 
              file=file.path(path,isa["assay.filenames"][[i]]), 
              row.names=FALSE, col.names=TRUE, 
              quote=TRUE, sep="\t", na="\"\"")
}

getStudyFilename <- function(isa, assay.filename){
  j <- which( lapply(isa["assay.filenames.per.study"], function(x)  (assay.filename %in% x)) == TRUE)
 return(isa@study.filenames[[j]])
}

getStudyFilenameIndex <- function(isa, assay.filename){
  j <- which( lapply(isa["assay.filenames.per.study"], function(x)  (assay.filename %in% x)) == TRUE)
  return(j)
}

#AnnnotatedDataFrame - previous phenoData object
getAnnotatedDataFrameAssay <- function(isa, assay.filename)
{
  i <- which(names(isa["assay.files"])==assay.filename)
      
  dataf <- as.data.frame(isa@factors[[i]])
  
  colnames(dataf) <-  names(isa@factors[[i]])
  row.names(dataf) <- isa@samples 
  
  dataf.desc <-  as.data.frame(names(isa@factors[[i]]))
  
  pdata <- new("AnnotatedDataFrame", data = dataf, varMetadata = dataf.desc)
  
  return(pdata)
}

## Check whether all the files exist
checkFilesExist = function(files){
  all(sapply(files, files))
}