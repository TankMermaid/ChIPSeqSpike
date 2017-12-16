.validateParameters <- function(type, reverse, outputFolder = NULL, 
        bigWigfilePath = NA, firstPart = TRUE){
    
    if(!is.null(outputFolder)){
        
        if(!file.exists(outputFolder))
            stop("The specified output folder does not exist.")
        
        if(isTRUE(all.equal(str_sub(outputFolder, -1), "/")))
            stop("The path to the output folder should not end by '/'")
    }
    
    if(firstPart){
        
        if(!(isTRUE(all.equal(type, "endo")) || 
                    isTRUE(all.equal(type, "exo"))))
            stop("Accepted types are endo and exo.")
        
        if(isTRUE(all.equal(type, "exo")) && reverse)
            stop("Exogenous scaling factor cannot be reverted.")
    }else{
        
        if(reverse && !length(grep("RPM", bigWigfilePath)))
            stop("RPM normalization must be performed before reverting it")
        
        if(reverse && !length(grep("BGSub", bigWigfilePath)))
            stop("Input subtraction should be performed before ",
                            "reverting RPM normalization")
        
        if(isTRUE(all.equal(type, "exo")) && !length(grep("reverted", 
                        bigWigfilePath)))
            stop("Exogenous scaling factor should be applied when RPM ",
                            "normalization has been reverted.")
    }
}


.retrieveSF <- function(object, type){
    
    if(isTRUE(all.equal(type, "endo")))
        scaling_factor <- getScalingFactor(object)
    else
        scaling_factor <- getExogenousScalingFactor(object)
    
    return(scaling_factor)
}


.modifyBigWigName <- function(object, name, outputFolder = NULL){
    
    if(is.null(outputFolder)){
        output_folder_bigwig <- paste0(dirname(getBigWigFile(object)), "/")
        output_bigWig <- strsplit(basename(getBigWigFile(object)), "\\.")
        output_bigWig <- paste0(output_folder_bigwig, output_bigWig[[1]][1], 
                "-", name, ".", output_bigWig[[1]][2])
    }else{
        output_bigWig <- strsplit(basename(getBigWigFile(object)), "\\.")
        output_bigWig <- paste0(outputFolder, "/", output_bigWig[[1]][1], "-", 
                name, ".", output_bigWig[[1]][2])
    }
    
    return(output_bigWig)
}


.computeScaling <- function(object, outputFolder = NULL, verbose = TRUE, 
        reverse = FALSE, type = "endo"){
    
    scaling_factor <- .retrieveSF(object, type)
    
    if(verbose) message("\t Reading experiment bigWig file.")
    
    bigWig_file <- import(getBigWigFile(object), format="BigWig")
    
    if(reverse)
    {
        if(verbose) message("\t Reverse RPM.")
        score(bigWig_file) <- score(bigWig_file) / scaling_factor
        if(verbose) message("\t Output RPM reverted bigWig file")
        output_bigWig <- .modifyBigWigName(object, "reverted", outputFolder)
    }else{
        
        if(verbose) message("\t Apply scaling factor")
        score(bigWig_file) <- score(bigWig_file) * scaling_factor
        
        if(verbose) message("\t Output bigWig file")
        if(isTRUE(all.equal(type, "endo")))
            output_bigWig <- .modifyBigWigName(object, "RPM", outputFolder)
        else
            output_bigWig <- .modifyBigWigName(object, "spiked", outputFolder)
    }
    
    export(bigWig_file, con = output_bigWig, format="BigWig")
    
    return(output_bigWig)
}


.computeScalingBoost <- function(object, verbose, reverse, type){
    
    scaling_factor <- .retrieveSF(object, type)
    currentData <- getLoadedData(object)
    
    if(reverse){
        if(verbose) message("\t Reverse RPM.")
        score(currentData) <- score(currentData)/scaling_factor
        
    }else{
        if(verbose) message("\t Apply scaling factor")
        score(currentData) <- score(currentData) * scaling_factor
    }
    
    return(currentData)
}


setMethod(
        
        f = "scaling",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, reverse = FALSE, type = "endo",
                verbose = TRUE, outputFolder = NULL){
            
            .validateParameters(type, reverse, outputFolder)
            
            if(!reverse && isTRUE(all.equal(type, "endo"))){
                if(verbose) message("Processing input")
                
                bigWigFile(theObject) <- 
                        .computeScaling(theObject, outputFolder, verbose, 
                                reverse, type)
            }
            
            
            experimentList(theObject) <- 
                    lapply(getExperimentList(theObject), function(experiment){
                                
                   .validateParameters(type, reverse, outputFolder, 
                           getBigWigFile(experiment), firstPart = FALSE)
              if(verbose)
                  message("Processing ", getExpName(experiment))
                                
              bigWigFile(experiment) <- .computeScaling(experiment, 
                      outputFolder, verbose, reverse, type)
              
              return(experiment)
              })
            
            return(theObject)
        })


setMethod(
        
        f = "scaling",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, reverse = FALSE, type = "endo",
                verbose = TRUE, outputFolder = NULL){
            
            .validateParameters(type, reverse)
            
            if(!reverse && isTRUE(all.equal(type, "endo"))){
                if(verbose) message("Processing input")
                
                loadedData(theObject) <- .computeScalingBoost(theObject, 
                        verbose, reverse, type)
                
                bigWigFile(theObject) <- .modifyBigWigName(theObject, "RPM", 
                                                            outputFolder)
            }
            
            
            experimentList(theObject) <- 
                  lapply(getExperimentList(theObject), function(experiment){
                                
                     .validateParameters(type, reverse, outputFolder,
                             getBigWigFile(experiment), firstPart = FALSE)
                                
                             if(verbose)
                                message("Processing ", getExpName(experiment))
                                
                             loadedData(experiment) <- .computeScalingBoost(
                                     experiment, verbose, reverse, type)
                             
                             bigWigFile(experiment) <- 
                                     .modifyBigWigName(experiment, 
                                     if(reverse) "reverted" else{
                                                 
                                            if(isTRUE(all.equal(type, "endo")))
                                                "RPM" 
                                            else
                                                "spiked"
                                             }, outputFolder)
                                return(experiment)
                            })
            
            return(theObject)
        })


.loadScalingOnList <- function(theObject, reverse, type, verbose, 
        outputFolder){
    
    datasetList(theObject) <- lapply(getDatasetList(theObject), 
            function(object){
                return(scaling(object,reverse,type,verbose, outputFolder))
            })
    
    return(theObject)
}


setMethod(
        
        f = "scaling",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, reverse = FALSE, type = "endo",
                verbose = TRUE, outputFolder = NULL){
            
            .loadScalingOnList(theObject, reverse, type, verbose, 
                    outputFolder)
        }
)


setMethod(
        
        f = "scaling",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, reverse = FALSE, type = "endo",
                verbose = TRUE, outputFolder = NULL){
            
            .loadScalingOnList(theObject, reverse, type, verbose, 
                    outputFolder)
        }
)
