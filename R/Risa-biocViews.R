bioconductor.version <- 2.12

getPackagesInBiocView <- function(view, 
                                  reposUrl = c("BiocSoftware", "BiocAnnotationData", "BiocExperimentData"),
                                  biocVersion = bioconductor.version) {
        
#    suppressPackageStartupMessages(require("biocViews"))
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

suggestBiocPackage <- function(isa, bioc.version = bioconductor.version ){
    
  mapping.file <- system.file("extdata","ISA-BiocViews-Mapping.csv", package="Risa")
  
  mapping <- read.csv(mapping.file, header=TRUE, fill=TRUE, na.strings = "NA", row.names=NULL)
  
  assay.measurement.types <- isa["assay.measurement.types"]
  
  matching.measurement.types <- mapping[["ISA.measurement.type"]] == assay.measurement.types
  
  measurement.types.views <- mapping[["BiocViews"]][ matching.measurement.types ]
  
  assay.technology.types <- isa["assay.technology.types"]
  
  reduced.mapping <- mapping[which(matching.measurement.types),]
  
  technology.types.views <- reduced.mapping[["BiocViews.1"]][   reduced.mapping[["ISA.technology.type"]] == assay.technology.types ]
  
  views <- setdiff(union(technology.types.views, measurement.types.views),"")
  
  package.list <- lapply(views, getPackagesInBiocView)
  
  return(package.list)
}