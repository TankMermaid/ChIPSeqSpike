#####################
## getter methods
#####################


setGeneric(
        
        name = "getBam",
        
        def = function(theObject){
            
            standardGeneric("getBam");
            
        });


setGeneric(
        
        name = "getExogenousBam",
        
        def = function(theObject){
            
            standardGeneric("getExogenousBam");
            
        });

setGeneric(
        
        name = "getBigWigFile",
        
        def = function(theObject){
            
            standardGeneric("getBigWigFile");
        });


setGeneric(
        
        name = "getExperimentList",
        
        def = function(theObject){
            
            standardGeneric("getExperimentList");
        });


setGeneric(
        
        name = "getExpName",
        
        def = function(theObject){
            
            standardGeneric("getExpName");
        });


setGeneric(
        
        name = "getScalingFactor",
        
        def = function(theObject){
            
            standardGeneric("getScalingFactor");
        });


setGeneric(
        
        name = "getExogenousScalingFactor",
        
        def = function(theObject){
            
            standardGeneric("getExogenousScalingFactor");
        });


setGeneric(
        
        name = "getCount",
        
        def = function(theObject){
            
            standardGeneric("getCount");
        });


setGeneric(
        
        name = "getExoCount",
        
        def = function(theObject){
            
            standardGeneric("getExoCount");
        });


setGeneric(
        
        name = "spikeSummary",
        
        def = function(theObject){
            
            standardGeneric("spikeSummary");
        });


setGeneric(
        
        name = "getAverageBindingValues",
        
        def = function(theObject){
            
            standardGeneric("getAverageBindingValues");
        });


setGeneric(
        
        name = "getMatBindingValues",
        
        def = function(theObject){
            
            standardGeneric("getMatBindingValues");
        });


setGeneric(
        
        name = "getLoadedData",
        
        def = function(theObject){
            
            standardGeneric("getLoadedData");
        });


setGeneric(
        
        name = "getRatio",
        
        def = function(theObject){
            
            standardGeneric("getRatio");
        });


#####################
## setter methods
#####################


setGeneric(
        
        name = "scalingFactor<-",
        
        def = function(theObject, value){
            
            standardGeneric("scalingFactor<-");
        });


setGeneric(
        
        name = "exogenousScalingFactor<-",
        
        def = function(theObject, value){
            
            standardGeneric("exogenousScalingFactor<-");
        });


setGeneric(
        
        name = "count<-",
        
        def = function(theObject, value){
            
            standardGeneric("count<-");
        });


setGeneric(
        
        name = "exoCount<-",
        
        def = function(theObject, value){
            
            standardGeneric("exoCount<-");
        });


setGeneric(
        
        name = "bigWigFile<-",
        
        def = function(theObject, value){
            
            standardGeneric("bigWigFile<-");
        });


setGeneric(
        
        name = "averageBindingValues<-",
        
        def = function(theObject, value){
            
            standardGeneric("averageBindingValues<-");
        });


setGeneric(
        
        name = "matBindingValues<-",
        
        def = function(theObject, value){
            
            standardGeneric("matBindingValues<-");
        });


setGeneric(
        
        name = "loadedData<-",
        
        def = function(theObject, value){
            
            standardGeneric("loadedData<-");
        });


#####################
## Other methods
#####################


setGeneric(
        
        name = "plotProfile",
        
        def = function(theObject, legends = FALSE, colVec = NULL, 
                notScaled = FALSE){
            
            standardGeneric("plotProfile");
        });


setGeneric(
        
        name = "plotTransform",
        
        def = function(theObject, legends = FALSE, colVec = NULL, 
                separateWindows = FALSE){
            
            standardGeneric("plotTransform");
        });


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
            
            standardGeneric("plotCor");
        });


setGeneric(
        
        name = "boxplotSpike",
        
        def = function(theObject, col = NULL, rawFile = FALSE, rpmFile = FALSE,
                bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, ylab = NULL,
                outline = TRUE, violinPlot = FALSE, notch = TRUE, 
                mean_with_sd = FALSE, mean = FALSE, median = FALSE, 
                boxplot = FALSE, jitter = FALSE, plot = TRUE, verbose = FALSE){
            
            standardGeneric("boxplotSpike");
        });


setGeneric(
        
        name = "plotHeatmaps",
        
        def = function(theObject, location = "start", transformType = "spiked",
                legend = TRUE, plot_scale = "no", sort_rows = "decreasing", 
                nb_of_groups = 1, clustering_method = "none", 
                include_exp_vec = NULL, auto_scale = FALSE, 
                raster_value = TRUE, col_value = "blue", ...){
            
            standardGeneric("plotHeatmaps");
        });


setGeneric(
        
        name = "estimateScalingFactors",
        
        def = function(theObject, paired = FALSE, verbose = TRUE){
            
            standardGeneric("estimateScalingFactors");
        });


setGeneric(
        
        name = "scaling",
        
        def = function(theObject, reverse = FALSE, type = "endo", 
                verbose = TRUE, outputFolder = NULL){
            
            standardGeneric("scaling");
        });



setGeneric(
        
        name = "inputSubtraction",
        
        def = function(theObject, verbose = TRUE){
            
            standardGeneric("inputSubtraction");
        });


setGeneric(
        
        name = "exportBigWigs",
        
        def = function(theObject, verbose = TRUE){
            
            standardGeneric("exportBigWigs");
        });


setGeneric(
        
        name = "extractBinding",
        
        def = function(theObject, gff_vec, genome, binsize = 50, before = 2000,
                after=2000, mean_or_median = "mean", 
                interpolation_number = 100, interpolation_average = 10000, 
                ignore_strand = FALSE, verbose = FALSE){
            
            standardGeneric("extractBinding");
        });


setGeneric(
        
        name = "buildMeanMatrix",
        
        def = function(theObject, ...){
            
            standardGeneric("buildMeanMatrix");
        });


#####################
## Methods defined from existing generics
#####################

setMethod(
        
        f = "[[",
        
        signature = c("ChIPSeqSpikeDatasetList", "ANY", "ANY"),
        
        definition = function(x, i, j, ...){
            
            return(x@datasetList[[i]]);
        }
);


setMethod(
        
        f = "[[",
        
        signature = c("ChIPSeqSpikeDataset", "ANY", "ANY"),
        
        definition = function(x, i, j, ...){
            
            return(x@experimentList[[i]]);
        }
);


setMethod(
        
        f = "[[",
        
        signature = c("ChIPSeqSpikeDatasetListBoost", "ANY", "ANY"),
        
        definition = function(x, i, j, ...){
            
            return(x@datasetList[[i]]);
        }
);


setMethod(
        
        f = "[[<-",
        
        signature = c("Experiment", "ANY", "ANY"),
        
        definition=function(x, i, j, value) {
            
            x <- value;
            return(x);
        }
);


setMethod(
        
        f = "[[<-",
        
        signature = c("ChIPSeqSpikeDataset", "ANY", "ANY"),
        
        definition=function(x, i, j, value) {
            
            x@experimentList[[i]] <- value;
            return(x);
        }
);


setMethod(
        
        f = "[[<-",
        
        signature = c("ChIPSeqSpikeDatasetList", "ANY", "ANY"),
        
        definition=function(x, i, j, value) {
            
            x@datasetList[[i]] <- value;
            return(x);
        }
);
