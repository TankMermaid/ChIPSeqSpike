.computeBamCount <- function(bam_file, paired){
    
    result <- countBam(bam_file, 
            param=ScanBamParam(scanBamFlag(isPaired= paired, 
                            isUnmappedQuery =FALSE)))$records;
    
    return(result);
};


.computeScalingFactor <- function(bam_count){
    
    return(1/(bam_count/1000000));
};


setMethod(
        
        f = "estimateScalingFactors",
        
        signature = "Experiment",
        
        definition = function(theObject, paired = FALSE, 
                verbose = TRUE){
            
            if(verbose) message("\t Computing endogenous scaling factor.");
            
            endo_count <- .computeBamCount(getBam(theObject), paired); 
            
            if(!endo_count) 
                stop("Endogenous bam file contains no aligned reads.");
            
            endo_SF <- .computeScalingFactor(endo_count);
            
            if(verbose) message("\t Computing exogenous scaling factor.");
            
            exo_count <- .computeBamCount(getExogenousBam(theObject),
                    paired);
            
            if(!exo_count) 
                stop("Exogenous bam file contains no aligned reads.");
            
            exo_SF <- .computeScalingFactor(exo_count);
            
            scalingFactor(theObject) <- endo_SF;
            exogenousScalingFactor(theObject) <- exo_SF;
            count(theObject) <- endo_count;
            exoCount(theObject) <- exo_count;
            
            return(theObject);
        });


setMethod(
        
        f = "estimateScalingFactors",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject, paired = FALSE, 
                verbose = TRUE){
            
            return(callNextMethod(theObject));
        });


setMethod(
        
        f = "estimateScalingFactors",
                
                signature = "ChIPSeqSpikeDataset",
                
                definition = function(theObject, paired = FALSE, 
                        verbose = TRUE){
                    
                    if(verbose) message("Computing scaling factor for input");
                    count(theObject) <- 
                            .computeBamCount(getBam(theObject), 
                                    paired);
                    
                    if(!getCount(theObject)) 
                        stop("input bam file contains no aligned reads.");
                    
                    scalingFactor(theObject) <- 
                            .computeScalingFactor(getCount(
                                            theObject));
                    
                    theObject@experimentList <-
                            lapply(theObject@experimentList, 
                                    function(experiment){
                                        
                                        if(verbose)
                                            message("Processing ", 
                                                    getExpName(experiment));
                                        return(
                                 estimateScalingFactors(experiment, paired, 
                                         verbose));
                                    });
                    
                    return(theObject);
                });


setMethod(
        
        f = "estimateScalingFactors",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, paired = FALSE, 
                verbose = TRUE){
            
            if(verbose) message("Computing scaling factor for input");
            count(theObject) <- 
                    .computeBamCount(getBam(theObject), 
                            paired);
            
            if(!getCount(theObject)) 
                stop("input bam file contains no aligned reads.");
            
            scalingFactor(theObject) <- 
                    .computeScalingFactor(getCount(
                                    theObject));
            
            theObject@experimentListLoaded <-
                    lapply(theObject@experimentListLoaded, 
                            function(experiment){
                                
                                if(verbose)
                                    message("Processing ", 
                                            getExpName(experiment));
                                return(
                                      estimateScalingFactors(experiment, 
                                              paired, verbose));
                            });
            
            return(theObject);
        });



.loadOnList <- function(object, paired, verbose){
    
    object@datasetList <- lapply(object@datasetList, 
            function(object2){
                
                return(estimateScalingFactors(object2, paired, verbose)); 
            });
    
    return(object);
};

setMethod(f = "estimateScalingFactors",
        signature = "ChIPSeqSpikeDatasetList",
        definition = function(theObject, paired = FALSE, verbose = TRUE){
            .loadOnList(theObject, paired, verbose);
        });

setMethod(f = "estimateScalingFactors",
        signature = "ChIPSeqSpikeDatasetListBoost",
        definition = function(theObject, paired = FALSE, verbose = TRUE){
            .loadOnList(theObject, paired, verbose);
        });
