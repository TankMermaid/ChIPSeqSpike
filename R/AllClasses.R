## The ChIPSeqSpikeCore class contains information about the input.
## The slots are:
##  1. Path to the input bam file.
##  2. Path to the input BigWig file.
##  3. Input scaling factor (only endogenous, used for input subtraction).
##  4. Input count (used to compute the scaling factor).
##  5. A slot that will contains binding values for performing plots.

ChIPSeqSpikeCore <- setClass(
        
        Class = "ChIPSeqSpikeCore",
        
        slots = c(inputBam = "character",
                inputBigWig = "character",
                inputScalingFactor = "numeric",
                inputCount = "numeric",
                plotSetArrayList = "list",
                matBindingValList = "list"),
        
        validity = function(object){
            
            .validatePath(object@inputBam)
            .validatePath(object@inputBigWig)
            .validateBigWig(object@inputBigWig)
            .validateSFAndCount(object@inputScalingFactor)
            .validateSFAndCount(object@inputCount)
            .validatePlotSetArrayList(object@plotSetArrayList)
            .validateMatBindingValList(object@matBindingValList)
        })


## The ChIPSeqSpikeDatasetList class intend to handle the case where different
## input files are used in the info file. It consists in a list of 
## ChIPSeqSpikeDataset objects. One object per input file is created by the 
## function spikeDataset(). see constructors.R

ChIPSeqSpikeDatasetList <- setClass(
        
        Class = "ChIPSeqSpikeDatasetList",
        
        slots = c(datasetList = "list"),
        
        validity = function(object){
            
            if(length(object@datasetList) < 1)
                stop("must contain at least one ChIPSeqSpikeDataset object.")
            else if(!all(vapply(object@datasetList, function(x) 
                                inherits(x, "ChIPSeqSpikeDataset"), logical(1))
            ))
               stop("all objects in list must be of ChIPSeqSpikeDataset", 
                       " type.")
            
        })


## The ChIPSeqSpikeDataset class inherits from ChIPSeqSpikeCore and adds a list
## of Experiment objects.
## The slots are (in addition of ChIPSeqSpikeCore slots):
##  1. List of Experiment objects.

ChIPSeqSpikeDataset <- setClass(
        
        Class = "ChIPSeqSpikeDataset",
        
        slots = c(experimentList = "list"),
        
        contains = "ChIPSeqSpikeCore",
        
        validity = function(object){
            
            if(length(object@experimentList) < 1)
                stop("must contain at least one Experiment object.")
            else if(!all(vapply(object@experimentList, function(x) 
                                inherits(x, "Experiment"), logical(1))))
                stop("all objects in list must be of Experiment type.")
            
        })


## The ChIPSeqSpikeDatasetListBoost class intend to handle the case where 
## different input files are used in the info file and where the user whishes 
## to perform the analysis in boost mode. It consists in a list of 
## ChIPSeqSpikeDatasetBoost objects. One object per input file is created by 
## the function spikeDataset(). see constructors.R

ChIPSeqSpikeDatasetListBoost <- setClass(
        
        Class = "ChIPSeqSpikeDatasetListBoost",
        
        slots = c(datasetList = "list"),
        
        validity = function(object){
            
            if(length(object@datasetList) < 1)
              stop("must contain at least one ChIPSeqSpikeDatasetBoost ",
                      "object.")
            else if(!all(vapply(object@datasetList, function(x) 
                                inherits(x, "ChIPSeqSpikeDatasetBoost"), 
                            logical(1))))
          stop("all objects in list must be of ChIPSeqSpikeDatasetBoost ",
                  "type.")
            
        })


## The ChIPSeqSpikeDatasetBoost class inherits from the ChIPSeqSpikeCore 
## class. The aims is to not have to read/write bigWigs files at each 
## operation. This offers the possibility of faster processing if the user has 
## sufficiant memory power. It contains the loaded input file and a list of 
## ExperimentLoaded which itself inherits from the Experiment class.

ChIPSeqSpikeDatasetBoost <- setClass(
        
        Class = "ChIPSeqSpikeDatasetBoost",
        
        slots = c(experimentListLoaded = "list",
                inputBigWigLoaded = "GRanges"),
        
        contains = "ChIPSeqSpikeCore",
        
        validity = function(object){
            
            if(length(object@experimentListLoaded) < 1)
                stop("Must contain at least one ExperimentLoaded object.")
            else if(!all(vapply(object@experimentListLoaded, function(x) 
                                inherits(x, "ExperimentLoaded"), logical(1))))
                stop("All objects in list must be of ExperimentLoaded type.")
        })


## This class defines the Experiment object which is the basis of the analysis.
## It contains all elements required to start scaling by spike-in:
##  1. Path to bam file endogenous
##  2. Path to bam file exogenous
##  3. Path to original fixed step wig file 
##  4. Endogenous scaling factor
##  5. Exogenous scaling factor
##  6. Endogenous count (used to compute the scaling factor).
##  7. Exogenous count (used to compute the scaling factor).

Experiment <- setClass(
        
        Class = "Experiment",
        
        slots = c(
                endogenousBam = "character",
                exogenousBam = "character",
                bigWigFile = "character",
                expName = "character",
                endogenousScalingFactor = "numeric",
                exogenousScalingFactor = "numeric",
                endoCount = "numeric",
                exoCount = "numeric"),
        
        validity = function(object){
            
            .validatePath(object@endogenousBam)
            .validatePath(object@exogenousBam)
            .validatePath(object@bigWigFile)
            .validateBigWig(object@bigWigFile)
            
            .validateSFAndCount(object@endogenousScalingFactor)
            .validateSFAndCount(object@exogenousScalingFactor)
            .validateSFAndCount(object@endoCount)
            .validateSFAndCount(object@exoCount)
        }
)



.validatePath <- function(value){
    
    if(!file.exists(value))
        stop(value, " is not a valid path")
    
    if(length(gregexpr(pattern="\\.", basename(value))[[1]]) > 1)
        stop(value, " should contain only one point for the extension.")
    
}

.validateSFAndCount <- function(value){
    
    if(!identical(value, numeric(0)) && value < 0)
        stop(deparse(substitute(value)), " should be a positive number or 0.")
    
}

.validateBigWig <- function(value){
    
    extension <- strsplit(basename(value), "\\.")[[1]][2]
    
    if(!identical(extension, "bw"))
        stop("Wig files should be in bigWig format.")
}

.validatePlotSetArrayList <- function(value){
    
    if(length(value)){
        
        if(!all(vapply(value, function(x)inherits(x, "PlotSetArray"), 
                        logical(1))))
            stop("All objects in setArray list must be of type PlotSetArray.")
    }
}

.validateMatBindingValList <- function(value){
    
    if(length(value))
        if(!all(vapply(value, function(x) inherits(x, "matrix"), logical(1))))
            stop("All objects in matBindingValList must be of type matrix.")
}

## The ExperimentLoaded class inherits from Experiment. It just adds the loaded
## values in the form of a GRanges object.

ExperimentLoaded <- setClass(
        
        Class = "ExperimentLoaded",
        
        slots = c(loadedBigWigFile = "GRanges"),
        
        contains = "Experiment")
