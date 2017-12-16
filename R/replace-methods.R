## Setters for Experiment objects


setReplaceMethod(
        
        f = "scalingFactor",
        
        signature = "Experiment",
        
        definition = function(theObject, value){
            
            if(!is.numeric(value)) stop("Scaling factor should be numeric.")
            
            theObject@endogenousScalingFactor <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "exogenousScalingFactor",
        
        signature = "Experiment",
        
        definition = function(theObject, value){
            
            if(!is.numeric(value)) stop("Scaling factor should be numeric.")
            
            theObject@exogenousScalingFactor <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "count",
        
        signature = "Experiment",
        
        definition = function(theObject, value){
            
            if(!is.numeric(value)) stop("Counts should be numeric.")
            
            theObject@endoCount <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "exoCount",
        
        signature = "Experiment",
        
        definition = function(theObject, value){
            
            if(!is.numeric(value)) stop("Counts should be numeric.")
            
            theObject@exoCount <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "bigWigFile",
        
        signature = "Experiment",
        
        definition = function(theObject, value){
            
            theObject@bigWigFile <- value
            validObject(theObject)
            return(theObject)
        })


## Setters for ExperimentLoaded objects


setReplaceMethod(
        
        f = "bigWigFile",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject, value){
            
            theObject@bigWigFile <- value
            return(theObject)
        })


## Setters for ChIPSeqSpikeCore objects


setReplaceMethod(
        
        f = "scalingFactor",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject, value){
            
            if(!is.numeric(value)) stop("Scaling factor should be numeric.")
            
            theObject@inputScalingFactor <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "count",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject, value){
            
            if(!is.numeric(value)) stop("Counts should be numeric.")
            
            theObject@inputCount <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "bigWigFile",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject, value){
            
            theObject@inputBigWig <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "averageBindingValues",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject, value){
            
            if(!inherits(value, "list"))
            {
                stop("SetArrayList should be of type list.")
            }
            
            if(!length(value))
            {
                stop("The list of binding values is empty.")
            }
            
            theObject@plotSetArrayList <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "matBindingValues",
        
        signature = "ChIPSeqSpikeCore",
        
        definition = function(theObject, value){
            
            if(!inherits(value, "list"))
            {
                stop("matBinding should be of type list.")
            }
            
            if(!length(value))
            {
                stop("The list of binding values is empty.")
            }
            
            theObject@matBindingValList <- value
            validObject(theObject)
            return(theObject)
        })


## Setters for ChIPSeqSpikeDataset object


setReplaceMethod(
        
        f = "scalingFactor",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject,value))
        })


setReplaceMethod(
        
        f = "count",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject,value))
        })


setReplaceMethod(
        
        f = "bigWigFile",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject, value))
        })


setReplaceMethod(
        
        f = "averageBindingValues",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject, value))
        })


setReplaceMethod(
        
        f = "matBindingValues",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject, value))
        })


setReplaceMethod(
        
        f = "experimentList",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, value){
            
            theObject@experimentList <- value
            validObject(theObject)
            return(theObject)
        })


## Setters for ChIPSeqSpikeDatasetBoost object


setReplaceMethod(
        
        f = "scalingFactor",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject,value))
        })


setReplaceMethod(
        
        f = "count",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject,value))
        })


setReplaceMethod(
        
        f = "loadedData",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, value){
            
            theObject@inputBigWigLoaded <- value
            validObject(theObject)
            return(theObject)
        })


setReplaceMethod(
        
        f = "loadedData",
        
        signature = "ExperimentLoaded",
        
        definition = function(theObject, value){
            
            theObject@loadedBigWigFile <- value
            return(theObject)
        })


setReplaceMethod(
        
        f = "bigWigFile",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, value){
            
            theObject@inputBigWig <- value
            return(theObject)
        })


setReplaceMethod(
        
        f = "averageBindingValues",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject, value))
        })


setReplaceMethod(
        
        f = "matBindingValues",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, value){
            
            return(callNextMethod(theObject, value))
        })


setReplaceMethod(
        
        f = "experimentList",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, value){
            
            theObject@experimentListLoaded <- value
            return(theObject)
        })


## Setters for ChIPSeqSpikeDatasetList object


setReplaceMethod(
        
        f = "datasetList",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, value){
            
            theObject@datasetList <- value
            validObject(theObject)
            return(theObject)
        })


## Setters for ChIPSeqSpikeDatasetListBoost object


setReplaceMethod(
        
        f = "datasetList",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, value){
            
            theObject@datasetList <- value
            validObject(theObject)
            return(theObject)
        })
