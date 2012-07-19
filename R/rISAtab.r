isa_syntax <- list(
  investigation_prefix="i_",
  study_prefix="s_",
  assay_prefix="a_",
  study_identifier="Study Identifier",
  study_assay_filename="Study Assay File Name"
  )

## This function only works if the zip file does not contain a directory (but the ISA-TAB files themselves)
isatab2bioczip = function(zip, path = getwd())
{
  d = unzip(zipfile = zip, exdir = extract <- path)
  isaobj = isatab2bioc(path)
  return(isaobj)
}##end function isatab2bioczip


isatab2bioc = function(path = getwd())
{
  #### Parse ISATab files
  d = dir(path)

  ## TODO any use case with multiple investigation files?
  ## Investigation filename
  ifilename = grep("i_", d, value=TRUE)
  if (length(ifilename)==0)
    stop("Did not find any investigation file.")
  else if (!file.exists(file.path(path, ifilename)))
    stop("Did not find investigation file: ", ifilename)
  
  ## Reading in investigation file into a data frame
  ifile = read.table(file.path(path, ifilename), sep="\t", fill=TRUE)

  ## Study Identifiers  - as a list of strings
  sidentifiers = ifile[grep("Study Identifier", ifile[,1], useBytes=TRUE),][2][[1]]
  
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
  sfiles = lapply(sfilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE))
  
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
  afiles = lapply(afilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE))
  names(afiles) <- afilenames
  # afiles_per_study is a list (one element per study) of lists (one element per assay) 
  afiles_per_study = lapply(seq_len(length(afilenames_per_study)), 
                            function(j) (lapply(seq_len(length(afilenames_per_study[[j]])),
                                    function(i) read.table(file.path(path,afilenames_per_study[[j]][[i]]), sep="\t", header=TRUE, stringsAsFactors=FALSE))))
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
  
  ####TODO build two data structures relating
  ### a sample name with all the data files
  ### a data file with all the sample names
  
  ## List of data filenames with assay filenames as keys
  dfilenames_per_assay = lapply(afiles, function(i) i[,grep("Data.File", colnames(i))])
  
  ## Identifying what sample is studied in which assay
  ## assays is a list of data frames (one for each assay file)
  assays = lapply(seq_len(length(sfiles)), 
                          function(j) (lapply(seq_len(length(afiles)), 
                                              function(i) sfiles[[j]]$Sample.Name %in% afiles[[i]]$Sample.Name)))
  
  samples = unlist(lapply(sfiles, function(i) i[,grep("Sample.Name", colnames(i))]))
  
  samples_per_assay_filename = lapply(seq_len(length(afiles)), 
                                            function(i) afiles[[i]]$Sample.Name)
  names(samples_per_assay) <- afilenames
  
  samples_per_study <- lapply(seq_len(length(sfiles)),
                                function(i) sfiles[[i]]$Sample.Name)
  names(samples_per_study) <- sidentifiers
  
  assayfilename_per_sample <- unlist(lapply(seq_len(length(samples)), 
                             function(j) lapply(seq_len(length(afilenames)), 
                                    function(i)   if (samples[[j]] %in% afiles[[i]]$Sample.Name) {
                                                          afilenames[[i]]
                                                  }
                                                )))
  
  sample_to_rawdatafile <- afiles[[1]][,c('Sample.Name','Raw.Data.File')]
  sample_to_rawdatafile <- merge(sample_to_rawdatafile[ !duplicated(sample_to_rawdatafile$'Sample.Name'), ], sample_to_rawdatafile[ duplicated(sample_to_rawdatafile$'Sample.Name'), ], all=TRUE)  
  
  sample_to_assayname <-afiles[[1]][,c('Sample.Name','Assay.Name')]
  sample_to_assayname <- merge(sample_to_assayname[ !duplicated(sample_to_assayname$'Sample.Name'), ], sample_to_assayname[ duplicated(sample_to_assayname$'Sample.Name'), ], all=TRUE)
  
  rawdatafile_to_sample <- afiles[[1]][,c('Raw.Data.File','Sample.Name')]
  rawdatafile_to_sample <- merge(rawdatafile_to_sample[ !duplicated(rawdatafile_to_sample$'Raw.Data.File'), ], rawdatafile_to_sample[ duplicated(rawdatafile_to_sample$'Raw.Data.File'), ], all=TRUE)
  
  assayname_to_sample <- afiles[[1]][,c('Assay.Name','Sample.Name')]
  assayname_to_sample <- merge(assayname_to_sample[ !duplicated(assayname_to_sample$'Assay.Name'), ], assayname_to_sample[ duplicated(assayname_to_sample$'Assay.Name'), ], all=TRUE)
  
 

  ## Adding the study file content to the isa object
  ## metadata kept into a data.frame - maintains study files and assay files info
  #metadata = cbind(sfiles, assays)
	
  isaobject <- list(
    investigation_filename=ifilename,
    investigation_file=ifile,
    study_identifiers=sidentifiers,
    study_filenames=sfilenames,
    study_files=sfiles,
    assay_filenames=afilenames,
    assay_filenames_per_study=afilenames_per_study,
    assay_technology_types=assay_tech_types,
    data_filenames_per_assay=dfilenames_per_assay,
    samples_per_assay=samples_per_assay,
    assays_per_sample=assays_per_sample,
    sample_to_rawdatafile=sample_to_rawdatafile,
    sample_to_assayname=sample_to_assayname,
    rawdatafile_to_sample=rawdatafile_to_sample,
    assayname_to_sample=assayname_to_sample
    )
  return(isaobject)
  
}##end function isatab2bioc

processAssayType = function(isa)
{
  for(i in seq_len(length(dfilenames)))
  {
      #############################################################################
      if (assay_types[i] == "DNA microarray")
      {
        ## Raw and processed data filenames
        rawfilenames = if ("Array.Data.File" %in% colnames(dfilenames[[i]])) dfilenames[[i]][,"Array.Data.File"] else NULL
        procfilenames = if("Derived.Array.Data.File" %in% colnames(dfiles[[i]])) dfilenames[[i]][,"Derived.Array.Data.File"] else NULL
			   
        ## URL for ADF (Array Design Format) file
        urladf = paste("http://www.ebi.ac.uk/microarray-as/ae/files/", unique(afiles[[i]][,"Array.Design.REF"]), "/", unique(afiles[[i]][,"Array.Design.REF"]), ".adf.txt", sep="")
        adffilename = file.path(path,unique(afiles[[i]][,"Array.Design.REF"]))
        adf_download = download.file(urladf, adffilename, mode="wb")

        ## List containing rawfiles, sdrf, idf, adf & directory containing the files
        ## as required by {ArrayExpress} magetab2bioc function
        files = list(path = path,
                  rawfiles = rawfilenames,
                  procfile = procfilenames,
                  sdrf = afilenames[[i]],
                  idf = ifilename,
                  adf = basename(adffilename))
			   
        if (is.null(dim(dfilenames[[i]])[2]))
            ## No processed files
            isa[[i]] = try(ae2bioc(files)) 
          else {
              raw = try(ae2bioc(files))
              ## TODO more testing
              ## The following issues an R CMD check warning,
              ## no visible binding for global variable ‘procol’
              ## likely to be a true issue:
              #cn = getcolproc(files)
              #procol = cn[1]
              proc = try(procset(files, procol = procol))
              isa[[i]] = list(raw=raw, processed=proc)}
      }## end microarray
      #############################################################################
      
      #############################################################################
      else if (assay_types[i] == "flow cytometry")
      {
          pd = try(read.AnnotatedDataFrame(file.path(path, afilenames[i]),row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
          sampleNames(pd) = pd$Raw.Data.File
                                        #if(all(dfiles[[i]] %in% dir(path)) && all(dfiles[[i]] %in% sampleNames(pd)))
          isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]], phenoData=pd, path=path)))
                                        #if(!all(dfiles[[i]] %in% dir(path)))
                                        #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% dir(path)], phenoData=pd, path=path)))				
                                        #if(!all(dfiles[[i]] %in% sampleNames(pd)))
                                        #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% sampleNames(pd)], phenoData=pd, path=path)))				
      }## end flow cytometry
      #############################################################################
      
      
      #############################################################################
      else if (assay_types[i] == "mass spectrometry")
        {
          if ("Raw.Spectral.Data.File" %in% colnames(dfilenames[[i]]))
          {
              #mass spectrometry files
              msfiles = dfilenames[[i]]$Raw.Spectral.Data.File
              
              pd = try(read.AnnotatedDataFrame(file.path(path, afilenames[i]),
                row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                varMetadata.char = "$", quote="\""))
              
              sampleNames(pd) = pd$Raw.Spectral.Data.File

              if (length(grep("Factor.Value", colnames(metadata))) != 0) {
              
                ## If there are explicit factors, use them
                sclass = metadata[ which(metadata$Sample.Name %in% pd$Sample.Name), grep("Factor.Value", colnames(metadata))[1]]
                isa[[i]] = xcmsSet(files=msfiles, sclass=sclass)
              } else {
                  ## Otherwise just use what was there
                  isa[[i]] = try(xcmsSet(msfiles, phenoData=pData(pd)))
              }

          }# end Raw.Spectral.Data.File			
        }## end mass spectrometry
      #############################################################################
      else{
        stop("Study Assay Technology Type '", types[i], "' not yet supported in the Risatab package")
      }
  }## end for on dfiles

		
  names(isa) = do.call(paste, list("isa", seq_len(length(isa)), sep=""))
  isaobj = list(metadata, isa)

  names(isaobj) = c("metadata","data")
  return(isaobj)

}##end function processAssayType

## Check whether all the files exist
checkFilesExist = function(files){
  all(sapply(files, files))
}