

isa_syntax <- list(
  investigation_prefix="i_",
  study_prefix="s_",
  assay_prefix="a_",
  study_identifier="Study Identifier",
  study_assay_filename="Study Assay File Name"
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

  ## TODO any use case with multiple investigation files?
  ## Investigation filename
  ifilename = grep("i_", d, value=TRUE)
  if (length(ifilename)==0)
    stop("Did not find any investigation file at folder ", path)
  else if (!file.exists(file.path(path, ifilename)))
    stop("Did not find investigation file: ", ifilename)
  
  ## Reading in investigation file into a data frame
  ifile = read.table(file.path(path, ifilename), sep="\t", fill=TRUE, na.strings = "NA")
  #row.names(ifile) <- ifile[[1]]
  #ifile <- ifile[,2:length(ifile)]

  ## Study Identifiers  - as a list of strings
  sidentifiers = ifile[grep("Study Identifier", ifile[,1], useBytes=TRUE),][2][[1]]
                 #ifile["Study Identifier",]
  
  ## Study filenames (one or more)
  sfilenames = unlist(sapply(ifile[grep("Study File Name", ifile[,1], useBytes=TRUE),], function(i) grep("s_", i, value=TRUE, useBytes=TRUE)))
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
  afilenames = unlist(sapply(ifile[grep("Study Assay File Name", ifile[,1], useBytes=TRUE),], function(i) grep("a_", i, value=TRUE, useBytes=TRUE)))
  
  #getting afilenames associated with studies
  afilenames.df = ifile[grep("Study Assay File Name", ifile[,1], useBytes=TRUE),]
  afilenames.matrix = apply(afilenames.df,c(1,2),function(row) grep("a_",row, value=TRUE))  
  afilenames.lists = split(afilenames.matrix, row(afilenames.matrix, as.factor=TRUE))
  afilenames_per_study = lapply(seq_len(length(afilenames.lists)), function(i) Filter(function(j) !identical(character(0), j), afilenames.lists[[i]]))
  names(afilenames_per_study) <- sidentifiers
  
  ## Reading in assay files 
  # afiles is a list of data frames (containing all the assay files)
  afiles <- lapply(afilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE,  check.names=FALSE))
  
  #afiles_headers <- lapply(afilenames, function(i) read.table(file.path(path, i), check.names=FALSE, sep='\t', header=TRUE, 
  #                                 stringsAsFactors=FALSE))
  #afiles_headers <- lapply(afiles_headers, function(i) colnames(i))
  #afiles_headers
  
  names(afiles) <- afilenames
  # afiles_per_study is a list (one element per study) of lists (one element per assay) 
  afiles_per_study = lapply(seq_len(length(afilenames_per_study)), 
                            function(j) (lapply(seq_len(length(afilenames_per_study[[j]])),
                                    function(i) read.table(file.path(path,afilenames_per_study[[j]][[i]]), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))))
  names(afiles_per_study) <- sidentifiers

  ## Assay technology types
  #data frame with types
  assay_tech_types = ifile[grep("Study Assay Technology Type$", ifile[,1], useBytes=TRUE),]    
  #remove empty types - results in a list of types
  assay_tech_types = na.omit(assay_tech_types[assay_tech_types != ""])
  #remove headers
  assay_tech_types = assay_tech_types[ assay_tech_types != "Study Assay Technology Type"]

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
  
  samples = unlist(lapply(sfiles, function(i) i[,grep("Sample.Name", colnames(i))]))
  
  samples_per_assay_filename = lapply(seq_len(length(afiles)), 
                                            function(i) afiles[[i]][["Sample Name"]])
  names(samples_per_assay_filename) <- afilenames
  
  samples_per_study <- lapply(seq_len(length(sfiles)),
                                function(i) sfiles[[i]][["Sample Name"]])
  names(samples_per_study) <- sidentifiers
  
  assay_filenames_per_sample <- unlist(lapply(seq_len(length(samples)), 
                             function(j) lapply(seq_len(length(afilenames)), 
                                    function(i)   if (samples[[j]] %in% afiles[[i]][["Sample Name"]]) {
                                                          afilenames[[i]]
                                                  }
                                                )))
  
  data_col_names = lapply(seq_len(length(afiles)),
                      function(i) if ('Raw Data File' %in% colnames(afiles[[i]])){
                                  'Raw Data File'
                                  }else if ('Free Induction Decay Data File' %in% colnames(afiles[[i]])){
                                    'Free Induction Decay Data File'
                                  }else if ('Array.Data.File' %in% colnames(afiles[[i]])){
                                    'Array Data File'
                                  }else if ('Raw Spectral Data File' %in% colnames(afiles[[i]])){
                                    'Raw Spectral Data File'
                                  })
                                    
  
  sample_to_rawdatafile <- lapply( seq_len(length(afiles)), 
                                  function(i) afiles[[i]][,c('Sample Name',data_col_names[[i]])] )
  sample_to_rawdatafile <- lapply(seq_len(length(afiles)), function(i)
     merge(sample_to_rawdatafile[[i]][ !duplicated(sample_to_rawdatafile[[i]]$'Sample Name'), ], sample_to_rawdatafile[[i]][ duplicated(sample_to_rawdatafile[[i]]$'Sample Name'), ], all=TRUE))  
  
  sample_to_assayname <-lapply( afiles,
                                function(i) i[,c('Sample Name',grep('Assay Name', colnames(i), value=TRUE))])
  sample_to_assayname <- lapply(seq_len(length(afiles)), function(i)
          merge(sample_to_assayname[[i]][ !duplicated(sample_to_assayname[[i]]$'Sample Name'), ], sample_to_assayname[[i]][ duplicated(sample_to_assayname[[i]]$'Sample Name'), ], all=TRUE))
  
  rawdatafile_to_sample <- lapply( seq_len(length(afiles)), 
                                   function(i) afiles[[i]][,c(data_col_names[[i]],'Sample Name')] )
  rawdatafile_to_sample <- lapply(seq_len(length(afiles)), function(i)
         merge(rawdatafile_to_sample[[i]][ !duplicated(rawdatafile_to_sample[[i]][[data_col_names[[i]]]]), ], rawdatafile_to_sample[[i]][ duplicated(rawdatafile_to_sample[[i]][[data_col_names[[i]]]]), ], all=TRUE))
  
  assayname_to_sample <- lapply( afiles,
                                 function(i) i[,c(grep('Assay Name', colnames(i), value=TRUE),'Sample Name')])
  assayname_to_sample <- lapply(seq_len(length(afiles)), function(i)
          merge(assayname_to_sample[[i]][ !duplicated(assayname_to_sample[[i]][,c(grep('Assay Name', colnames(assayname_to_sample[[i]]), value=TRUE))]), ], 
                assayname_to_sample[[i]][  duplicated(assayname_to_sample[[i]][,c(grep('Assay Name', colnames(assayname_to_sample[[i]]), value=TRUE))]), ], 
                all=TRUE))
  
 

  ## Adding the study file content to the isa object
  ## metadata kept into a data.frame - maintains study files and assay files info
  #metadata = cbind(sfiles, assays)
	
  isaobject <- list(
    path=path,
    investigation_filename=ifilename,
    investigation_file=ifile,
    study_identifiers=sidentifiers,
    study_filenames=sfilenames,
    study_files=sfiles,
    assay_filenames=afilenames,
    assay_filenames_per_study=afilenames_per_study,
    assay_files=afiles,
    assay_files_per_study=afiles_per_study,
    assay_technology_types=assay_tech_types,
    data_filenames=dfilenames_per_assay,
    samples=samples,
    samples_per_assay_filename=samples_per_assay_filename,
    assay_filenames_per_sample=assay_filenames_per_sample,
    sample_to_rawdatafile=sample_to_rawdatafile,
    sample_to_assayname=sample_to_assayname,
    rawdatafile_to_sample=rawdatafile_to_sample,
    assayname_to_sample=assayname_to_sample
    )
  return(isaobject)
  
}##end function isatab2bioc



### specific function to deal with assays whose technology type is mass spectrometry using the xcms package
### it returns an xcmsSet
processAssayXcmsSet = function(isa, assay_filename, ...){
  for(i in seq_len(length(isa$assay_filenames))){
    
    if (isa$assay_filenames[[i]]==assay_filename){
      
      if ("Raw Spectral Data File" %in% colnames(isa$data_filenames[[i]]))
      {
        #mass spectrometry files
        msfiles = isa$data_filenames[[i]][["Raw Spectral Data File" ]]
        
        pd = try(read.AnnotatedDataFrame(file.path(isa$path, isa$assay_filenames[i]),
                                         row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                                         varMetadata.char = "$", quote="\""))
        
        sampleNames(pd) = pd$Raw.Spectral.Data.File
        
        if (length(grep("Factor.Value", colnames(isa$assay_files[[i]]))) != 0) {
          ## If there are explicit factors, use them
          sclass = isa$assay_files[[i]][ which(isa$assay_files[[i]][["Sample Name"]] %in% pd$Sample.Name), grep("Factor.Value", colnames(isa$assay_files[[i]]))[1]]
          
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
addAssayMetadata = function(isa, assay_filename, col_name, values){
  
  assay_file <- isa$assay_files [[ assay_filename ]]
  if (length(values)==1){
    values <- c(rep(values,nrow(assay_file)))
  }else if (length(values)!=nrow(assay_file)){
    stop("Wrong number of values to be added to the assay file")
  }
  assay_file [ colnames(assay_file) == col_name ] <- values
  isa$assay_files [[ assay_filename ]] <- assay_file
  return(isa)
}

### TODO fix quotes when writing files
### ADD COMMENT - written with R 
write.isatab = function(isa){
  write.investigation.file(isa)
  for(i in seq_len(length(isa$study_filenames))){
    write.table(isa$study_files[[i]], 
                file=isa$study_filenames[[i]], 
                row.names=FALSE, col.names=FALSE, 
                quote=FALSE, sep="\t", na="\"\"")
  }
  for(i in seq_len(length(isa$assay_filenames))){
    write.assay.file(isa, isa$assay_files[[i]])
  }
  
}

write.investigation.file = function(isa){
  write.table(isa$investigation_file, 
              file=isa$investigation_filename, 
              row.names=FALSE, col.names=FALSE, 
              quote=TRUE, sep="\t", na="\"\"")
}

write.assay.file = function(isa, assay_filename){
  i <- which(names(isa$assay_files)==assay_filename)
  assay_file <- isa$assay_files[[assay_filename ]]
  write.table(assay_file, 
              file=isa$assay_filenames[[i]], 
              row.names=FALSE, col.names=TRUE, 
              quote=TRUE, sep="\t", na="\"\"")
}

processAssayType = function(isa)
{
  for(i in seq_len(length(isa$assay_filenames)))
  {
      #############################################################################
      if (isa$assay_technology_types[i] == technology_types$microarray)
      {
      #  ## Raw and processed data filenames
      #  rawfilenames = if ("Array.Data.File" %in% colnames(isa$data_filenames[[i]])) isa$data_filenames[[i]][,"Array.Data.File"] else NULL
      #  procfilenames = if("Derived.Array.Data.File" %in% colnames(isa$data_files[[i]])) isa$data_filenames[[i]][,"Derived.Array.Data.File"] else NULL
			   
        ## URL for ADF (Array Design Format) file
      #  urladf = paste("http://www.ebi.ac.uk/microarray-as/ae/files/", unique(isa$asay_files[[i]][,"Array.Design.REF"]), "/", unique(assay_files[[i]][,"Array.Design.REF"]), ".adf.txt", sep="")
      #  adffilename = file.path(path,unique(isa$assay_files[[i]][,"Array.Design.REF"]))
      #  adf_download = download.file(urladf, adffilename, mode="wb")

        ## List containing rawfiles, sdrf, idf, adf & directory containing the files
        ## as required by {ArrayExpress} magetab2bioc function
       # files = list(path = path,
       #            rawfiles = rawfilenames,
       #          procfile = procfilenames,
       #           sdrf = isa$assay_filenames[[i]],
       #            idf = isa$investigation_filename,
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
      #else if (isa$assay_technology_types[i] == technology_types$fc)
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
      else if (isa$assay_technology_types[i] == technology_types$ms)
      {
          if ("Raw.Spectral.Data.File" %in% colnames(isa$data_filenames[[i]]))
          {
              #mass spectrometry files
              msfiles = isa$data_filenames[[i]]$Raw.Spectral.Data.File
              
              pd = try(read.AnnotatedDataFrame(file.path(isa$path, isa$assay_filenames[i]),
                row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                varMetadata.char = "$", quote="\""))
              
              sampleNames(pd) = pd$Raw.Spectral.Data.File

              if (length(grep("Factor.Value", colnames(isa$assay_files[[i]]))) != 0) {
                ## If there are explicit factors, use them
                sclass = isa$assay_files[[i]][ which(isa$assay_files[[i]]$Sample.Name %in% pd$Sample.Name), grep("Factor.Value", colnames(isa$assay_files[[i]]))[1]]
                
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
        stop("Study Assay Technology Type '", isa$assay_tech_types[i], "' not yet supported in the ISA2R package")
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