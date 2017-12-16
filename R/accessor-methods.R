## Methods to access slots

setMethod(
        
        f = "getBam",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@endogenousBam)
        })


setMethod(
        
        f = "getBam",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getBam",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject){
            
            return(theObject@inputBam)
        })


setMethod(
        
        f = "getBam",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getBam",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getExogenousBam",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@exogenousBam)
        })


setMethod(
        
        f = "getBigWigFile",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@bigWigFile)
        })


setMethod(
        
        f = "getBigWigFile",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getBigWigFile",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject){
            
            return(theObject@inputBigWig)
        })


setMethod(
        
        f = "getBigWigFile",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getBigWigFile",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getBigWigFile",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject){
            
            bigwig_vec <- unlist(lapply(theObject@datasetList, 
                    function(dataset){
                        
                        bw_vec <- getBigWigFile(dataset)
                        bw_vec <- c(bw_vec, unlist(
                                        lapply(dataset@experimentList, 
                                                function(exp){
                                                    return(getBigWigFile(exp))
                                                    }
                            )))
                        
                        return(bw_vec)
                        
                    }))
    
            return(bigwig_vec)
        })


setMethod(
        
        f = "getExperimentListBigWigs",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            bigWigFiles_vec <- unlist(lapply(theObject@experimentList, 
                            function(x){return(getBigWigFile(x))}))
            
            return(bigWigFiles_vec)
        })


setMethod(
        
        f = "getExperimentListBigWigs",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            bigWigFiles_vec <- unlist(lapply(theObject@experimentListLoaded, 
                            function(x){return(getBigWigFile(x))}))
            
            return(bigWigFiles_vec)
        })


setMethod(
        
        f = "getExpName",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@expName)
        })


setMethod(
        
        f = "getExpName",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            result <- unlist(lapply(theObject@experimentList, function(exp){
                                return(getExpName(exp))
                            }))
            return(as.character(result))
        })


setMethod(
        
        f = "getExpName",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            result <- unlist(lapply(theObject@experimentListLoaded, 
                            function(exp){
                                return(getExpName(exp))
                            }))
            return(as.character(result))
        })


setMethod(
        
        f = "getExpName",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getScalingFactor",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@endogenousScalingFactor)
        })


setMethod(
        
        f = "getScalingFactor",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getScalingFactor",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject){
            
            return(theObject@inputScalingFactor)
        })

setMethod(
        
        f = "getScalingFactor",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })

setMethod(
        
        f = "getScalingFactor",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })



setMethod(
        
        f = "getExogenousScalingFactor",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@exogenousScalingFactor)
        })


setMethod(
        
        f = "getExogenousScalingFactor",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getCount",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@endoCount)
        })


setMethod(
        
        f = "getCount",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getCount",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject){
            
            return(theObject@inputCount)
        })

setMethod(
        
        f = "getCount",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })

setMethod(
        
        f = "getCount",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getExoCount",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            return(theObject@exoCount)
        })


setMethod(
        
        f = "getExoCount",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getAverageBindingValues",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject){
            
            return(theObject@plotSetArrayList)
        })


setMethod(
        
        f = "getAverageBindingValues",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getAverageBindingValues",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getMatBindingValues",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject){
            
            return(theObject@matBindingValList)
        })


setMethod(
        
        f = "getMatBindingValues",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getMatBindingValues",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })


setMethod(
        
        f = "getLoadedData",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(theObject@inputBigWigLoaded)
        })


setMethod(
        
        f = "getLoadedData",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(theObject@loadedBigWigFile)
        })


setMethod(
        
        f = "getDatasetList",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject){
            
            return(theObject@datasetList)
        })


setMethod(
        
        f = "getDatasetList",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject){
            
            return(theObject@datasetList)
        })


setMethod(
        
        f = "getExperimentList",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            return(theObject@experimentList)
        })


setMethod(
        
        f = "getExperimentList",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            return(theObject@experimentListLoaded)
        })


## Summary of counts and scaling factors of Experiment object

setMethod(
        
        f = "spikeSummary",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            result <- matrix(c(getScalingFactor(theObject),
                               getExogenousScalingFactor(theObject),
                               getCount(theObject),
                               getExoCount(theObject)), nrow=1)
            
            rownames(result) <- getExpName(theObject)
            colnames(result) <- c("endoScalFact", "exoScalFact", 
                                  "endoCount",    "exoCount")
                          
            return(result)
        })


setMethod(
        
        f = "spikeSummary",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            return(callNextMethod(theObject))
        })

## Summary of counts and scaling factors of ChIPSeqSpikeDataset

setMethod(
        
        f = "spikeSummary",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            stat_list <- mapply(function(x){
                        
                        return(spikeSummary(x))
                        
                    }, theObject@experimentList, SIMPLIFY = FALSE)
            
            input <- c(getScalingFactor(theObject), NA, 
                    getCount(theObject), NA)
            
            return(rbind(do.call(rbind, stat_list), input))
        })


## Summary of counts and scaling factors of ChIPSeqSpikeDatasetList and 
## ChIPSeqSpikeDatasetListBoost


.summaryDatasetList <- function(theObject){
    
    result <- lapply(theObject@datasetList, function(dataset){
                
                return(spikeSummary(dataset))
            })
    
    return(do.call(rbind, result))
}

setMethod(f = "spikeSummary",
        signature = "ChIPSeqSpikeDatasetList", .summaryDatasetList)

setMethod(f = "spikeSummary",
        signature = "ChIPSeqSpikeDatasetListBoost", .summaryDatasetList)


## Summary of counts and scaling factors of ChIPSeqSpikeDatasetBoost

setMethod(
        
        f = "spikeSummary",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            stat_list <- mapply(function(x){
                        
                        return(spikeSummary(x))
                        
                    }, theObject@experimentListLoaded, SIMPLIFY = FALSE)
             
            input <- c(getScalingFactor(theObject), NA, 
                    getCount(theObject), NA)
            
            return(rbind(do.call(rbind, stat_list), input))
        })


## Checking the endogenous dna ratio


.checkPercentages <- function(ratio_list, names_vec){
    
    invisible(mapply(function(ratio, expname){
                
                if(ratio < 2)
                    warning(expname, 
                            " contains less than 2% of endogenous dna:",
                            " The scaling might not work.", 
                            immediate. = TRUE)
                
                if(ratio > 25)
                    warning(expname, " contains more than ",
                            "25% of endogenous DNA.", immediate. = TRUE)
            }, ratio_list, names_vec))
}


setMethod(
        
        f = "getRatio",
        
        signature = "Experiment",
        
        definition = function(theObject){
            
            total_count <- countBam(getBam(theObject), param=ScanBamParam(
                            scanBamFlag(isPaired= FALSE, 
                                    isUnmappedQuery =NA)))$records
            exo_count <- getExoCount(theObject)
            percent_exo <- round((exo_count*100)/total_count, 
                    digits=1)
            
            return(as.numeric(percent_exo))
        })


setMethod(
        
        f = "getRatio",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject){
            
            callNextMethod(theObject)
        })


setMethod(
        
        f = "getRatio",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject){
            
            ratio_list <- lapply(theObject@experimentList, function(exp){
                        
                        return(getRatio(exp))
                    })
            
            names_vec <- names(ratio_list)
            .checkPercentages(ratio_list, names_vec)
            result <- do.call(rbind, ratio_list)
            colnames(result) <- c("Percentage Exo")
            return(result)
})


setMethod(
        
        f = "getRatio",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject){
            
            ratio_list <- lapply(theObject@experimentListLoaded, function(exp){
                        
                        return(getRatio(exp))
                    })
            
            names_vec <- names(ratio_list)
            .checkPercentages(ratio_list, names_vec)
            result <- do.call(rbind, ratio_list)
            colnames(result) <- c("Percentage Exo")
            return(result)
        })


.callRatio <- function(list_object){
    
    ratioResultList <- lapply(list_object, function(dataset){
                
                return(getRatio(dataset))
            })
    
    return(do.call(rbind,ratioResultList))
}


setMethod(
        
        f = "getRatio",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject){
            
            .callRatio(theObject@datasetList)
        })


setMethod(
        
        f = "getRatio",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject){
            
            .callRatio(theObject@datasetList)
        })


## Export all files for the object obtained in boost mode

setMethod(
        
        f = "exportBigWigs",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, verbose = TRUE){
            
            if(verbose)
                message("Exporting RPM scaled input.")
            
            export(getLoadedData(theObject), con = getBigWigFile(theObject), 
                    format="BigWig")
            
            invisible(lapply(theObject@experimentListLoaded, function(exp){
                        
                        if(verbose)
                            message("Writing spiked ", getExpName(exp))
                        
             export(getLoadedData(exp), con = getBigWigFile(exp), 
                                format="BigWig")
                    }))
        })


setMethod(
        
        f = "exportBigWigs",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, verbose = TRUE){
            
            invisible(lapply(theObject@datasetList, function(dataset){
                        
                        exportBigWigs(dataset, verbose)
                    }))
        })
