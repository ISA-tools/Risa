Sys.setlocale('LC_ALL','C') 

isatab.syntax <- list(
  investigation.prefix="i_",
  study.prefix="s_",
  assay.prefix="a_",
  study.identifier="Study Identifier",
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
  raw.spectral.data.file="Raw Spectral Data File",
  factor.name="Factor Name",
  factor.value="Factor Value",
  assay.name="Assay Name"
  )

technology.types <- list(
  microarray="DNA microarray",
  ms="mass spectrometry",
  fc="flow cytometry"
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
  isaobj = readISAtabFiles(path) 
  return(isaobj)
}##end function readISAtabZip

readISAtabFiles = function(path = getwd(), verbose=FALSE)
{
  if (verbose)
    message("Converting ISA-Tab dataset at ",path," into R objects...")
  isaobject <- new(Class="ISAtab",path=path)
  if (verbose)
    message("... done.")
  return(isaobject) 
}##end function readISAtabFiles



### ADD COMMENT - written with R
### ADD validation for samples
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

### TODO fix quotes when writing files
### ADD COMMENT - written with R 
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

processAssayType = function(isa)
{
  for(i in seq_len(length(isa["assay.filenames"])))
  {
      #############################################################################
      if (isa["assay.technology.types"][i] == technology.types$microarray)
      {
      #  ## Raw and processed data filenames
      #  rawfilenames = if ("Array.Data.File" %in% colnames(isa["data.filenames"][[i]])) isa["data.filenames"][[i]][,"Array.Data.File"] else NULL
      #  procfilenames = if("Derived.Array.Data.File" %in% colnames(isa["data_files"][[i]])) isa["data.filenames"][[i]][,"Derived.Array.Data.File"] else NULL
			   
        ## URL for ADF (Array Design Format) file
      #  urladf = paste("http://www.ebi.ac.uk/microarray-as/ae/files/", unique(isa["asay_files"][[i]][,"Array.Design.REF"]), "/", unique(assay.files[[i]][,"Array.Design.REF"]), ".adf.txt", sep="")
      #  adffilename = file.path(path,unique(isa["assay.files"][[i]][,"Array.Design.REF"]))
      #  adf_download = download.file(urladf, adffilename, mode="wb")

        ## List containing rawfiles, sdrf, idf, adf & directory containing the files
        ## as required by {ArrayExpress} magetab2bioc function
       # files = list(path = path,
       #            rawfiles = rawfilenames,
       #          procfile = procfilenames,
       #           sdrf = isa["assay.filenames"][[i]],
       #            idf = isa["investigation.filename"],
       #            adf = basename(adffilename))
			   
       #if (is.null(dim(dfilenames[[i]])[2]))
       #     ## No processed files
       #    isa[[i]] = try(ae2bioc(files)) 
       #   else {
       #        raw = try(ae2bioc(files))
       #        ## TODO more testing
              ## The following issues an R CMD check warning,
              ## no visible binding for global variable ‘procol’
              ## likely to be a true issue:
              #cn = getcolproc(files)
              #procol = cn[1]
      #        proc = try(procset(files, procol = procol))
      #        isa[[i]] = list(raw=raw, processed=proc)}
      }## end microarray
      #############################################################################
      
      #############################################################################
      #else if (isa["assay.technology.types"][i] == technology.types$fc)
      #{
      #    pd = try(read.AnnotatedDataFrame(file.path(path, afilenames[i]),row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
      #    sampleNames(pd) = pd$Raw.Data.File
      #                                  #if(all(dfiles[[i]] %in% dir(path)) && all(dfiles[[i]] %in% sampleNames(pd)))
      #    isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]], phenoData=pd, path=path)))
      #                                  #if(!all(dfiles[[i]] %in% dir(path)))
      #                                  #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% dir(path)], phenoData=pd, path=path)))				
      #                                  #if(!all(dfiles[[i]] %in% sampleNames(pd)))
      #                                  #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% sampleNames(pd)], phenoData=pd, path=path)))				
      #}## end flow cytometry
      #############################################################################
      
      
      #############################################################################
      else if (isa["assay.technology.types"][i] == technology.types$ms)
      {
          if ("Raw.Spectral.Data.File" %in% colnames(isa["data.filenames"][[i]]))
          {
              #mass spectrometry files
              msfiles = isa["data.filenames"][[i]]$Raw.Spectral.Data.File
              
              pd = try(read.AnnotatedDataFrame(file.path(isa["path"], isa["assay.filenames"][i]),
                row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                varMetadata.char = "$", quote="\""))
              
              sampleNames(pd) = pd$Raw.Spectral.Data.File

              if (length(grep("Factor.Value", colnames(isa["assay.files"][[i]]))) != 0) {
                ## If there are explicit factors, use them
                sclass = isa["assay.files"][[i]][ which(isa["assay.files"][[i]]$Sample.Name %in% pd$Sample.Name), grep("Factor.Value", colnames(isa["assay.files"][[i]]))[1]]
                
                wd <- getwd()
                setwd(isa["path"])
                xset = xcmsSet(files=msfiles, sclass=sclass)
                setwd(wd)
              } else {
                  wd <- getwd()
                  setwd(isa["path"])
                  ## Otherwise just use what was there
                  xset = try(xcmsSet(msfiles, phenoData=pData(pd)))
                  setwd(wd)
              }
              
             

          }# end Raw.Spectral.Data.File			
        }## end mass spectrometry
      #############################################################################
      else{
        stop("Study Assay Technology Type '", isa["assay.technology.types"][i], "' not yet supported in the Risa package")
      }
  }## end for on dfiles

		
  #names(isa) = do.call(paste, list("isa", seq_len(length(isa)), sep=""))
  #isaobj = list(metadata, isa)

  #names(isaobj) = c("metadata","data")
  return(isa)

}##end function processAssayType

## Check whether all the files exist
checkFilesExist = function(files){
  all(sapply(files, files))
}