
isatab.syntax <- list(
  investigation.prefix="i_",
  study.prefix="s_",
  assay.prefix="a_",
  study.identifier="Study Identifier",
  study.file.name="Study File Name",
  study.assay.file.name="Study Assay File Name",
  study.assay.technology.type="Study Assay Technology Type",
  sample.name="Sample Name",
  raw.data.file="Raw Data File",
  free.induction.decay.data.file="Free Induction Decay Data File",
  array.data.file="Array.Data.File",
  raw.spectral.data.file="Raw Spectral Data File",
  factor.name="Factor Name"
  )

technology_types <- list(
  microarray="DNA microarray",
  ms="mass spectrometry",
  fc="flow cytometry"
  )

## This function only works if the zip file does not contain a directory (but the ISA-TAB files themselves)
isatab2bioczip = function(zip, path = getwd(), verbose=FALSE)
{
  
  if (verbose)
    writeLines("Unzipping file...")
  d = unzip(zipfile = zip, exdir = extract <- path)
  if (verbose)
    writeLines("Converting ISA-Tab dataset into R objects")
  isaobj = isatab2bioc(path)
  return(isaobj)
}##end function isatab2bioczip

isatab2bioc = function(path = getwd(), verbose=FALSE)
{
  #### Parse ISATab files
  d = dir(path)

  ## Investigation filename
  ifilename = grep(isatab.syntax$investigation.prefix, d, value=TRUE)
  if (length(ifilename)==0)
    stop("Did not find any investigation file at folder ", path)
  else if (!file.exists(file.path(path, ifilename)))
    stop("Did not find investigation file: ", ifilename)
  
  ## Reading in investigation file into a data frame
  ifile = read.table(file.path(path, ifilename), sep="\t", fill=TRUE, na.strings = "NA")
  row.names(ifile) <- ifile[[1]]
  #ifile <- ifile[,2:length(ifile)]

  ## Study Identifiers  - as a list of strings
  sidentifiers = ifile[grep(isatab.syntax$study.identifier, ifile[,1], useBytes=TRUE),][2][[1]]
                 #ifile[isatab.syntax$study.identifier,]
  
  ## Study filenames (one or more)
  sfilenames = unlist(sapply(ifile[grep(isatab.syntax$study.file.name, ifile[,1], useBytes=TRUE),], function(i) grep(isatab.syntax$study.prefix, i, value=TRUE, useBytes=TRUE)))
  if (length(sidentifiers)!=length(sfilenames))
    stop("There are study files with no identifier assigned")
  ## Assign sidentifiers as names of the list sfilenames
  names(sfilenames) <- sidentifiers
     
  ## TODO pretty printing sfilenames
  ## Validation of existance of study files
  if (!all(sapply(sfilenames, function(i) file.exists(file.path(path, i)))))
    stop("Did not find some of the study files: ", sfilenames)
  
  ## Reading study files into a list of data frames
  sfiles = lapply(sfilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))
  
  ## List of assay filenames 
  #afilenames is a list with all the assay filenames (without association to studies)
  afilenames = unlist(sapply(ifile[grep(isatab.syntax$study.assay.file.name, ifile[,1], useBytes=TRUE),], function(i) grep(isatab.syntax$assay.prefix, i, value=TRUE, useBytes=TRUE)))
  
  #getting afilenames associated with studies
  afilenames.df = ifile[grep(isatab.syntax$study.assay.file.name, ifile[,1], useBytes=TRUE),]
  afilenames.matrix = apply(afilenames.df,c(1,2),function(row) grep(isatab.syntax$assay.prefix,row, value=TRUE))  
  afilenames.lists = split(afilenames.matrix, row(afilenames.matrix, as.factor=TRUE))
  afilenames_per_study = lapply(seq_len(length(afilenames.lists)), function(i) Filter(function(j) !identical(character(0), j), afilenames.lists[[i]]))
  names(afilenames_per_study) <- sidentifiers
  
  ## Reading in assay files 
  # afiles is a list of data frames (containing all the assay files)
  afiles <- lapply(afilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE,  check.names=FALSE))
  names(afiles) <- afilenames
  # afiles_per_study is a list (one element per study) of lists (one element per assay) 
  afiles_per_study = lapply(seq_len(length(afilenames_per_study)), 
                            function(j) (lapply(seq_len(length(afilenames_per_study[[j]])),
                                    function(i) read.table(file.path(path,afilenames_per_study[[j]][[i]]), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))))
  names(afiles_per_study) <- sidentifiers

  ## Assay technology types
  #data frame with types
  assay_tech_types = ifile[isatab.syntax$study.assay.technology.type, ] 
  #remove empty types - results in a list of types
  assay_tech_types = na.omit(assay_tech_types[assay_tech_types != ""])
  #remove headers
  assay_tech_types = assay_tech_types[ assay_tech_types != isatab.syntax$study.assay.technology.type]

  ## Validate number of assay technology types == number of afiles
  if (length(assay_tech_types)!=length(afiles)){
    stop("The number of assay files mismatches the number of assay types")
  }
  
  ## List of data filenames with assay filenames as keys
  dfilenames_per_assay = lapply(afiles, function(i) i[,grep("Data.File", colnames(i))])
  
  ## Identifying what sample is studied in which assay
  ## assays is a list of data frames (one for each assay file)
  assays = lapply(seq_len(length(sfiles)), 
                          function(j) (lapply(seq_len(length(afiles)), 
                                              function(i) sfiles[[j]]$Sample.Name %in% afiles[[i]]$Sample.Name)))
  
  samples = unlist(lapply(sfiles, function(i) i[,grep(isatab.syntax$sample.name, colnames(i))]))
  
  samples.per.assay.filename = lapply(seq_len(length(afiles)), 
                                            function(i) afiles[[i]][[isatab.syntax$sample.name]])
  names(samples.per.assay.filename) <- afilenames
  
  samples_per_study <- lapply(seq_len(length(sfiles)),
                                function(i) sfiles[[i]][[isatab.syntax$sample.name]])
  names(samples_per_study) <- sidentifiers
  
  assay.filenames.per.sample <- unlist(lapply(seq_len(length(samples)), 
                             function(j) lapply(seq_len(length(afilenames)), 
                                    function(i)   if (samples[[j]] %in% afiles[[i]][[isatab.syntax$sample.name]]) {
                                                          afilenames[[i]]
                                                  }
                                                )))
  
  data_col_names = lapply(seq_len(length(afiles)),
                      function(i) if (isatab.syntax$raw.data.file %in% colnames(afiles[[i]])){
                                     isatab.syntax$raw.data.file
                                  }else if (isatab.syntax$free.induction.decay.data.file %in% colnames(afiles[[i]])){
                                    isatab.syntax$free.induction.decay.data.file
                                  }else if (isatab.syntax$array.data.file %in% colnames(afiles[[i]])){
                                    isatab.syntax$array.data.file
                                  }else if (isatab.syntax$raw.spectral.data.file %in% colnames(afiles[[i]])){
                                    isatab.syntax$raw.spectral.data.file
                                  })
                                    
  
  sample.to.rawdatafile <- lapply( seq_len(length(afiles)), 
                                  function(i) afiles[[i]][,c('Sample Name',data_col_names[[i]])] )
  sample.to.rawdatafile <- lapply(seq_len(length(afiles)), function(i)
     merge(sample.to.rawdatafile[[i]][ !duplicated(sample.to.rawdatafile[[i]]$'Sample Name'), ], sample.to.rawdatafile[[i]][ duplicated(sample.to.rawdatafile[[i]]$'Sample Name'), ], all=TRUE))  
  
  sample.to.assayname <-lapply( afiles,
                                function(i) i[,c('Sample Name',grep('Assay Name', colnames(i), value=TRUE))])
  sample.to.assayname <- lapply(seq_len(length(afiles)), function(i)
          merge(sample.to.assayname[[i]][ !duplicated(sample.to.assayname[[i]]$'Sample Name'), ], sample.to.assayname[[i]][ duplicated(sample.to.assayname[[i]]$'Sample Name'), ], all=TRUE))
  
  rawdatafile.to.sample <- lapply( seq_len(length(afiles)), 
                                   function(i) afiles[[i]][,c(data_col_names[[i]],'Sample Name')] )
  rawdatafile.to.sample <- lapply(seq_len(length(afiles)), function(i)
         merge(rawdatafile.to.sample[[i]][ !duplicated(rawdatafile.to.sample[[i]][[data_col_names[[i]]]]), ], rawdatafile.to.sample[[i]][ duplicated(rawdatafile.to.sample[[i]][[data_col_names[[i]]]]), ], all=TRUE))
  
  assayname.to.sample <- lapply( afiles,
                                 function(i) i[,c(grep('Assay Name', colnames(i), value=TRUE),'Sample Name')])
  assayname.to.sample <- lapply(seq_len(length(afiles)), function(i)
          merge(assayname.to.sample[[i]][ !duplicated(assayname.to.sample[[i]][,c(grep('Assay Name', colnames(assayname.to.sample[[i]]), value=TRUE))]), ], 
                assayname.to.sample[[i]][  duplicated(assayname.to.sample[[i]][,c(grep('Assay Name', colnames(assayname.to.sample[[i]]), value=TRUE))]), ], 
                all=TRUE))
  
 

  ## Adding the study file content to the isa object
  ## metadata kept into a data.frame - maintains study files and assay files info
  #metadata = cbind(sfiles, assays)
	
  isaobject <- list(
    path=path,
    investigation.filename=ifilename,
    investigation.file=ifile,
    study.identifiers=sidentifiers,
    study.filenames=sfilenames,
    study.files=sfiles,
    assay.filenames=afilenames,
    assay.filenames.per.study=afilenames_per_study,
    assay.files=afiles,
    assay.files.per.study=afiles_per_study,
    assay.technology.types=assay_tech_types,
    data.filenames=dfilenames_per_assay,
    samples=samples,
    samples.per.assay.filename=samples.per.assay.filename,
    assay.filenames.per.sample=assay.filenames.per.sample,
    sample.to.rawdatafile=sample.to.rawdatafile,
    sample.to.assayname=sample.to.assayname,
    rawdatafile.to.sample=rawdatafile.to.sample,
    assayname.to.sample=assayname.to.sample
    )
  return(isaobject)
  
}##end function isatab2bioc



### specific function to deal with assays whose technology type is mass spectrometry using the xcms package
### it returns an xcmsSet
processAssayXcmsSet = function(isa, assay.filename, ...){
  for(i in seq_len(length(isa$assay.filenames))){
    
    if (isa$assay.filenames[[i]]==assay.filename){
      
      if ("Raw Spectral Data File" %in% colnames(isa$data.filenames[[i]]))
      {
        #mass spectrometry files
        msfiles = isa$data.filenames[[i]][[ isatab.syntax$raw.spectral.data.file ]]
        
        pd = try(read.AnnotatedDataFrame(file.path(isa$path, isa$assay.filenames[i]),
                                         row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                                         varMetadata.char = "$", quote="\""))
        
        sampleNames(pd) = pd$Raw.Spectral.Data.File
        
        if (length(grep("Factor.Value", colnames(isa$assay.files[[i]]))) != 0) {
          ## If there are explicit factors, use them
          sclass = isa$assay.files[[i]][ which(isa$assay.files[[i]][[isatab.syntax$sample.name]] %in% pd$Sample.Name), grep("Factor.Value", colnames(isa$assay.files[[i]]))[1]]
          
          wd <- getwd()
          setwd(isa$path)
          xset = xcmsSet(files=msfiles, sclass=sclass, ...)
          setwd(wd)
        } else {
          wd <- getwd()
          setwd(isa$path)
          ## Otherwise just use what was there
          xset = try(xcmsSet(msfiles, phenoData=pData(pd), ...))
          setwd(wd)
        }
        return(xset)
    }#if
    
  }#if 
  }#for
}#get.assay.xcmsSet


### ADD COMMENT - written with R
### ADD validation for samples
addAssayMetadata = function(isa, assay.filename, col.name, values){
  
  assay.file <- isa$assay.files [[ assay.filename ]]
  if (length(values)==1){
    values <- c(rep(values,nrow(assay.file)))
  }else if (length(values)!=nrow(assay.file)){
    stop("Wrong number of values to be added to the assay file")
  }
  assay.file [ colnames(assay.file) == col.name ] <- values
  isa$assay.files [[ assay.filename ]] <- assay.file
  return(isa)
}

### TODO fix quotes when writing files
### ADD COMMENT - written with R 
write.isatab = function(isa, path = getwd()){
  write.investigation.file(isa, path)
  for(i in seq_len(length(isa$study.filenames))){
    write.study.file(isa, isa$study.files[[i]], path)
  }
  for(i in seq_len(length(isa$assay.filenames))){
    write.assay.file(isa, isa$assay.files[[i]], path)
  }
  
}

write.investigation.file = function(isa, path = getwd()){
  write.table(isa$investigation.file, 
              file=file.path(path,isa$investigation.filename), 
              row.names=FALSE, col.names=FALSE, 
              quote=TRUE, sep="\t", na="\"\"")
}

write.study.file = function(isa, study.filename, path = getwd()){
  i <- which(isa$study.filenames==study.filename)
  study.file <- isa$study.files[[ i ]]
  write.table(study.file, 
              file=file.path(path,isa$study.filenames[[i]]), 
              row.names=FALSE, col.names=TRUE, 
              quote=TRUE, sep="\t", na="\"\"")
}

write.assay.file = function(isa, assay_filename, path = getwd()){
  i <- which(names(isa$assay.files)==assay_filename)
  assay_file <- isa$assay.files[[assay_filename ]]
  write.table(assay_file, 
              file=file.path(path,isa$assay.filenames[[i]]), 
              row.names=FALSE, col.names=TRUE, 
              quote=TRUE, sep="\t", na="\"\"")
}

processAssayType = function(isa)
{
  for(i in seq_len(length(isa$assay.filenames)))
  {
      #############################################################################
      if (isa$assay.technology.types[i] == technology_types$microarray)
      {
      #  ## Raw and processed data filenames
      #  rawfilenames = if ("Array.Data.File" %in% colnames(isa$data.filenames[[i]])) isa$data.filenames[[i]][,"Array.Data.File"] else NULL
      #  procfilenames = if("Derived.Array.Data.File" %in% colnames(isa$data_files[[i]])) isa$data.filenames[[i]][,"Derived.Array.Data.File"] else NULL
			   
        ## URL for ADF (Array Design Format) file
      #  urladf = paste("http://www.ebi.ac.uk/microarray-as/ae/files/", unique(isa$asay_files[[i]][,"Array.Design.REF"]), "/", unique(assay.files[[i]][,"Array.Design.REF"]), ".adf.txt", sep="")
      #  adffilename = file.path(path,unique(isa$assay.files[[i]][,"Array.Design.REF"]))
      #  adf_download = download.file(urladf, adffilename, mode="wb")

        ## List containing rawfiles, sdrf, idf, adf & directory containing the files
        ## as required by {ArrayExpress} magetab2bioc function
       # files = list(path = path,
       #            rawfiles = rawfilenames,
       #          procfile = procfilenames,
       #           sdrf = isa$assay.filenames[[i]],
       #            idf = isa$investigation.filename,
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
      #else if (isa$assay.technology.types[i] == technology_types$fc)
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
      else if (isa$assay.technology.types[i] == technology_types$ms)
      {
          if ("Raw.Spectral.Data.File" %in% colnames(isa$data.filenames[[i]]))
          {
              #mass spectrometry files
              msfiles = isa$data.filenames[[i]]$Raw.Spectral.Data.File
              
              pd = try(read.AnnotatedDataFrame(file.path(isa$path, isa$assay.filenames[i]),
                row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                varMetadata.char = "$", quote="\""))
              
              sampleNames(pd) = pd$Raw.Spectral.Data.File

              if (length(grep("Factor.Value", colnames(isa$assay.files[[i]]))) != 0) {
                ## If there are explicit factors, use them
                sclass = isa$assay.files[[i]][ which(isa$assay.files[[i]]$Sample.Name %in% pd$Sample.Name), grep("Factor.Value", colnames(isa$assay.files[[i]]))[1]]
                
                wd <- getwd()
                setwd(isa$path)
                xset = xcmsSet(files=msfiles, sclass=sclass)
                setwd(wd)
              } else {
                  wd <- getwd()
                  setwd(isa$path)
                  ## Otherwise just use what was there
                  xset = try(xcmsSet(msfiles, phenoData=pData(pd)))
                  setwd(wd)
              }
              
              isa$preprocessing[[i]] <- xset

          }# end Raw.Spectral.Data.File			
        }## end mass spectrometry
      #############################################################################
      else{
        stop("Study Assay Technology Type '", isa$assay_tech_types[i], "' not yet supported in the Risa package")
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