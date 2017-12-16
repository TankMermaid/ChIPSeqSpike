#####################
## getter methods
#####################


setGeneric(
        
        name = "getBam",
        
        def = function(theObject){
            
            standardGeneric("getBam")
        },
        signature = "theObject")


setGeneric(
        
        name = "getExogenousBam",
        
        def = function(theObject){
            
            standardGeneric("getExogenousBam")
        },
        signature = "theObject")

setGeneric(
        
        name = "getBigWigFile",
        
        def = function(theObject){
            
            standardGeneric("getBigWigFile")
        },
        signature = "theObject")


setGeneric(
        
        name = "getExperimentListBigWigs",
        
        def = function(theObject){
            
            standardGeneric("getExperimentListBigWigs")
        },
        signature = "theObject")


setGeneric(
        
        name = "getExpName",
        
        def = function(theObject){
            
            standardGeneric("getExpName")
        },
        signature = "theObject")


setGeneric(
        
        name = "getScalingFactor",
        
        def = function(theObject){
            
            standardGeneric("getScalingFactor")
        },
        signature = "theObject")


setGeneric(
        
        name = "getExogenousScalingFactor",
        
        def = function(theObject){
            
            standardGeneric("getExogenousScalingFactor")
        },
        signature = "theObject")


setGeneric(
        
        name = "getCount",
        
        def = function(theObject){
            
            standardGeneric("getCount")
        },
        signature = "theObject")


setGeneric(
        
        name = "getExoCount",
        
        def = function(theObject){
            
            standardGeneric("getExoCount")
        },
        signature = "theObject")


setGeneric(
        
        name = "spikeSummary",
        
        def = function(theObject){
            
            standardGeneric("spikeSummary")
        },
        signature = "theObject")


setGeneric(
        
        name = "getAverageBindingValues",
        
        def = function(theObject){
            
            standardGeneric("getAverageBindingValues")
        },
        signature = "theObject")


setGeneric(
        
        name = "getMatBindingValues",
        
        def = function(theObject){
            
            standardGeneric("getMatBindingValues")
        },
        signature = "theObject")


setGeneric(
        
        name = "getLoadedData",
        
        def = function(theObject){
            
            standardGeneric("getLoadedData")
        },
        signature = "theObject")


setGeneric(
        
        name = "getRatio",
        
        def = function(theObject){
            
            standardGeneric("getRatio")
        },
        signature = "theObject")


setGeneric(
        
        name = "getDatasetList",
        
        def = function(theObject){
            
            standardGeneric("getDatasetList")
        },
        signature = "theObject")


setGeneric(
        
        name = "getExperimentList",
        
        def = function(theObject){
            
            standardGeneric("getExperimentList")
        },
        signature = "theObject")


#####################
## setter methods
#####################


setGeneric(
        
        name = "scalingFactor<-",
        
        def = function(theObject, value){
            
            standardGeneric("scalingFactor<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "exogenousScalingFactor<-",
        
        def = function(theObject, value){
            
            standardGeneric("exogenousScalingFactor<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "count<-",
        
        def = function(theObject, value){
            
            standardGeneric("count<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "exoCount<-",
        
        def = function(theObject, value){
            
            standardGeneric("exoCount<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "bigWigFile<-",
        
        def = function(theObject, value){
            
            standardGeneric("bigWigFile<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "averageBindingValues<-",
        
        def = function(theObject, value){
            
            standardGeneric("averageBindingValues<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "matBindingValues<-",
        
        def = function(theObject, value){
            
            standardGeneric("matBindingValues<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "loadedData<-",
        
        def = function(theObject, value){
            
            standardGeneric("loadedData<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "datasetList<-",
        
        def = function(theObject, value){
            
            standardGeneric("datasetList<-")
        },
        signature = "theObject")


setGeneric(
        
        name = "experimentList<-",
        
        def = function(theObject, value){
            
            standardGeneric("experimentList<-")
        },
        signature = "theObject")


#####################
## Other methods
#####################


setGeneric(
        
        name = "plotProfile",
        
        def = function(theObject, legends = FALSE, colVec = NULL, 
                notScaled = FALSE){
            
            standardGeneric("plotProfile")
        },
        signature = "theObject")


setGeneric(
        
        name = "plotTransform",
        
        def = function(theObject, legends = FALSE, colVec = NULL, 
                separateWindows = FALSE){
            
            standardGeneric("plotTransform")
        },
        signature = "theObject")


setGeneric(
        
        name = "plotCor",
        
        def = function(theObject, rawFile = FALSE, rpmFile = FALSE, 
                bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, main = "", 
                add_contour = FALSE, method_cor = "spearman", nlevels = 10, 
                color_contour = "black", show_cor = TRUE, allOnPanel = TRUE, 
                method_scale = "none", method_corrplot = "circle", 
                heatscatterplot = TRUE, type_corrplot = "upper", 
                diag_corrplot = FALSE, separateWindows = FALSE, 
                verbose = FALSE, ...){
            
            standardGeneric("plotCor")
        },
        signature = "theObject")


setGeneric(
        
        name = "boxplotSpike",
        
        def = function(theObject, col = NULL, rawFile = FALSE, rpmFile = FALSE,
                bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, ylab = NULL,
                outline = TRUE, violinPlot = FALSE, notch = TRUE, 
                mean_with_sd = FALSE, mean = FALSE, median = FALSE, 
                boxplot = FALSE, jitter = FALSE, plot = TRUE, verbose = FALSE){
            
            standardGeneric("boxplotSpike")
        },
        signature = "theObject")


setGeneric(
        
        name = "plotHeatmaps",
        
        def = function(theObject, location = "start", transformType = "spiked",
                legend = TRUE, plot_scale = "no", sort_rows = "decreasing", 
                nb_of_groups = 1, clustering_method = "none", 
                include_exp_vec = NULL, auto_scale = FALSE, 
                raster_value = FALSE, col_value = "blue", ...){
            
            standardGeneric("plotHeatmaps")
        },
        signature = "theObject")


setGeneric(
        
        name = "estimateScalingFactors",
        
        def = function(theObject, paired = FALSE, verbose = TRUE){
            
            standardGeneric("estimateScalingFactors")
        },
        signature = "theObject")


setGeneric(
        
        name = "scaling",
        
        def = function(theObject, reverse = FALSE, type = "endo", 
                verbose = TRUE, outputFolder = NULL){
            
            standardGeneric("scaling")
        },
        signature = "theObject")



setGeneric(
        
        name = "inputSubtraction",
        
        def = function(theObject, verbose = TRUE){
            
            standardGeneric("inputSubtraction")
        },
        signature = "theObject")


setGeneric(
        
        name = "exportBigWigs",
        
        def = function(theObject, verbose = TRUE){
            
            standardGeneric("exportBigWigs")
        },
        signature = "theObject")


setGeneric(
        
        name = "extractBinding",
        
        def = function(theObject, gff_vec, genome, binsize = 50, before = 2000,
                after=2000, mean_or_median = "mean", 
                interpolation_number = 100, interpolation_average = 10000, 
                ignore_strand = FALSE, verbose = FALSE){
            
            standardGeneric("extractBinding")
        },
        signature = "theObject")


setGeneric(
        
        name = "buildMeanMatrix",
        
        def = function(theObject, ...){
            
            standardGeneric("buildMeanMatrix")
        },
        signature = "theObject")


#####################
## Methods defined from existing generics
#####################

setMethod(
        
        f = "[[",
        
        signature = c("ChIPSeqSpikeDatasetList", "ANY", "ANY"),
        
        definition = function(x, i, j, ...){
            
            return(x@datasetList[[i]])
        }
)


setMethod(
        
        f = "[[",
        
        signature = c("ChIPSeqSpikeDataset", "ANY", "ANY"),
        
        definition = function(x, i, j, ...){
            
            return(x@experimentList[[i]])
        }
)


setMethod(
        
        f = "[[",
        
        signature = c("ChIPSeqSpikeDatasetListBoost", "ANY", "ANY"),
        
        definition = function(x, i, j, ...){
            
            return(x@datasetList[[i]])
        }
)


setMethod(
        
        f = "[[<-",
        
        signature = c("Experiment", "ANY", "ANY"),
        
        definition=function(x, i, j, value) {
            
            x <- value
            return(x)
        }
)


setMethod(
        
        f = "[[<-",
        
        signature = c("ChIPSeqSpikeDataset", "ANY", "ANY"),
        
        definition=function(x, i, j, value) {
            
            x@experimentList[[i]] <- value
            return(x)
        }
)


setMethod(
        
        f = "[[<-",
        
        signature = c("ChIPSeqSpikeDatasetList", "ANY", "ANY"),
        
        definition=function(x, i, j, value) {
            
            x@datasetList[[i]] <- value
            return(x)
        }
)
