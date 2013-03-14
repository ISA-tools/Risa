bioconductor.version <- 2.12


getPackagesInBiocView <- function(view, 
                                  reposUrl = c("BiocSoftware", "BiocAnnotationData", "BiocExperimentData"),
                                  biocVersion = bioconductor.version) {
        
    suppressPackageStartupMessages(require("biocViews"))
    data(biocViewsVocab)
    
    reposUrl <- match.arg(reposUrl)
    
    biocMirror <- getOption("BioC_mirror", "http://bioconductor.org")
    
    #builds the paths for the URL
    biocPaths <- switch(reposUrl,         
                        BiocSoftware = "bioc",
                        BiocAnnotationData = "data/annotation", 
                        BiocExperimentData = "data/experiment",
                        )   
    
    reposUrl <- paste(biocMirror,
                 "packages",
                 biocVersion,
                 biocPaths, 
                 sep = "/") 
    
    bv <- getBiocViews(reposUrl, biocViewsVocab, "NoViewProvided")
    
    if (!view %in% names(bv)) {
      warning("BiocView ", view, " not found.")
      return(NULL)
    }  
    return(bv[[view]])
  }

suggestBiocPackage <- function(isa){
  
  #suppressPackageStartupMessages(require("biocViews"))
  #data(biocViewsVocab)
  
  #reposUrl = "http://www.bioconductor.org/packages/2.11/bioc"
  #biocViews <- getBiocViews(reposUrl, biocViewsVocab, "Software")
  #print(biocViews[1:2])

  #getPackageNEWS(prevRepos="http://www.bioconductor.org/packages/2.10/bioc",
  #             currRepos="http://www.bioconductor.org/packages/2.11/bioc",
  #             srcDir)

  #getPacksAndViews(reposUrl, biocViews, "NoViewProvided", false)

  #assayTechnologies <- getBiocSubViews(reposUrl, biocViewsVocab, "AssayTechnologies")
  #assayTechnologies$MassSpectrometry
  
  mapping.file <- system.file("extdata","ISA-BiocViews-Mapping.csv", package="Risa")
  
  mappping <- read.csv(mapping.file, header=FALSE, fill=TRUE, na.strings = "NA", row.names=NULL)
  
  assay.technology.types <- isa["assay.technology.types"]
  
  assay.measurement.types <- isa["assay.measurement.types"]
  
  
}