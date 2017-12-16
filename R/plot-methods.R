################
## Common plot methods
################



.removingNABigWig <- function(input_mat, verbose){
    
    ## Removing lines containing NA due to missing values for 
    ## instance at the end of chromosomes if using genes 
    ## annotations
    
    keep <-  complete.cases(input_mat)
    
    if(sum(!keep)){
        
        input_mat <- input_mat[keep,]
        
        if(verbose)
            message(sum(!keep), " element(s) removed ",
                    "due to missing values in the bigWig files.")
    }
    
    return(input_mat)
}


setMethod(
        
        f = "buildMeanMatrix",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, rawFile, rpmFile, bgsubFile, revFile,
                spiked, verbose, is_list = FALSE){
            
            matBindingValues_list <- getMatBindingValues(theObject)
            
            mat_mean_list <- lapply(getExperimentList(theObject), 
                    function(exp){
                        
                        exp_name <- getExpName(exp)
                        
                        ## Reading the different files to plot
                        if(verbose)
                            message("Reading ", exp_name)
                        
                        
                        ## Preparing the different file paths
                        name_vec <- vector()
                        
                        if(rawFile)
                            name_vec <- paste0(exp_name, "-raw")
                        
                        if(rpmFile)
                            name_vec <- c(name_vec, paste0(exp_name, "-RPM"))
                        
                        if(bgsubFile)
                           name_vec <- c(name_vec, paste0(exp_name, "-BGSub"))
                        
                        if(revFile)
                            name_vec <- c(name_vec, 
                                    paste0(exp_name, "-rev"))
                        
                        if(spiked)
                            name_vec <- c(name_vec, 
                                    paste0(exp_name, "-spike"))
                        
                        ## Retrieving mean binding values 
                        if(verbose)
                            message("\t Preparing mean value matrix.")
                        
                        mean_list <- lapply(name_vec, function(current_exp){
                                    
                                    if(verbose)
                                        message("\t\t", current_exp)
                                    
                                    current_mat <- 
                                          matBindingValues_list[[current_exp]]
                                    
                                    
                                    current_mean <- 
                                          apply(current_mat, MARGIN = 1, mean)
                                    
                                    return(current_mean)
                                })
                        
                        mean_mat <- do.call(cbind, mean_list)
                        colnames(mean_mat) <- name_vec
                        
                        return(mean_mat)
                    })
            
            input_mat <- do.call(cbind,mat_mean_list)
            
            if(!is_list)
                return(.removingNABigWig(input_mat, verbose))
            else
                return(input_mat)
            
            }
)


setMethod(
        
        f = "buildMeanMatrix",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, verbose, is_list = FALSE){
            
            matBindingValues_list <- getMatBindingValues(theObject)
            
            mat_mean_list <- lapply(getExperimentList(theObject), 
                    function(exp){
                        
                        exp_name <- getExpName(exp)
                        
                        ## Reading the different files to plot
                        if(verbose)
                            message("Reading ", exp_name, 
                                    " and computing means")
                        
                        current_mean <- 
                                apply(matBindingValues_list[[exp_name]], 
                                        MARGIN = 1, mean)
                        
                        return(current_mean)
                    })
            
            input_mat <- do.call(cbind,mat_mean_list)
            
            if(!is_list)
                return(.removingNABigWig(input_mat, verbose))
            else
                return(input_mat)
            
        }
)


################
## Plotting average profiles
################


.plotAverageProfiles <- function(legends, arrayList, indexes_track, title_vec,
        colVec){
    
    type_count_vec <- seq_len(length(arrayList))
    
    if(!legends)
        par(mfrow=c(2,2))
    else
        par(mfrow=c(3,2))
    
    invisible(mapply( function(bindingValues, type_count, indexes_track, 
                            main_title){
                        
                        plotAverage(if(!is.na(indexes_track[1])) 
                                            bindingValues[indexes_track] else 
                                            bindingValues,
                                keepratio = FALSE,
                                type = "full",
                                cex.lab = 9,
                                main = main_title[type_count],
                                legend=FALSE,
                                colvec = if(!is.null(colVec)) colVec
                                        else NULL,
                                cex.main = 10,
                                cex.axis = 10,
                                las = 1)
                    }, arrayList, type_count_vec, 
                    MoreArgs = list(indexes_track, title_vec)))
}


.plotLegends <- function(legends, colVec, experiment_number, expNamesVec){
    
    if(legends){
        
        plot.new()
        
        if(is.null(colVec))
            cols <- c("darkblue", "darkgreen", "darkred", "darkmagenta",
                    "darkgray","darkorange", "darkcyan", "black",
                    rainbow(experiment_number-8))
        else cols <- colVec
        
        legend("topleft", expNamesVec, fill = cols)
    }
}


.convertTypes <- function(list_location){
    
    lt <- c(pf="start", mf="midpoint", ef="end", af="composite")
    type_names <- lt[names(list_location)]
    
    return(type_names)
}


.plotSpiked <- function(arrayList, legends, expNamesVec, colVec, 
        experiment_number, spiked_indexes = NA){
    
    type_names <- .convertTypes(arrayList)
    
    .plotAverageProfiles(legends, arrayList, spiked_indexes, type_names, 
            colVec)
    
    .plotLegends(legends, colVec, experiment_number, expNamesVec)
    
}


.retrieveIndexes <- function(experiment_number, notScaled = FALSE){
    
    spiked_indexes <- if(notScaled){
                
                if(experiment_number > 5)
                    c(seq(from = 5, to = (5*experiment_number), by = 5),
                            seq(from = 3, to = (5*experiment_number), by = 5))
                else c(5,3)
            }else{
                if(experiment_number > 5)
                    seq(from = 5, to = (5*experiment_number), by = 5)
                else 5
            }
    return(spiked_indexes)
}


setMethod(
        
        f = "plotProfile",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, legends = FALSE, colVec = NULL){
            
            arrayList <- getAverageBindingValues(theObject)
            experiment_number <- length(getExperimentList(theObject))
            expNamesVec <- unlist(lapply(getExperimentList(theObject), 
                            getExpName), use.names = FALSE)
            
            .plotSpiked(arrayList, legends, expNamesVec, colVec, 
                    experiment_number)
            
        }
)


setMethod(
        
        f = "plotProfile",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, legends = FALSE, colVec = NULL, 
                notScaled = FALSE){
            
            arrayList <- getAverageBindingValues(theObject)
            experiment_number <- arrayList[[1]]$ntracks()/5
            
            spiked_indexes <- .retrieveIndexes(experiment_number, notScaled)
            expNamesVec <- unlist(lapply(getExperimentList(theObject), 
                            getExpName), use.names = FALSE)
            
            if(notScaled)
                expNamesVec <- c(paste0(expNamesVec,"-spiked"), expNamesVec)
            
           .plotSpiked(arrayList, legends, expNamesVec, colVec, 
                   experiment_number, spiked_indexes)
            
        }
)


setMethod(
        
        f = "plotProfile",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, legends = FALSE, colVec = NULL, 
                notScaled = FALSE){
            
            list_arrayList <- lapply(getDatasetList(theObject), 
                    function(object){
                        return(getAverageBindingValues(object))
                    })
            
            expNamesVec <- as.character(unlist(lapply(
                                    getDatasetList(theObject), 
                                    function(object){
                                        exp_name <- getExpName(object)
                                        if(notScaled)
                                            return(c(paste0(exp_name, 
                                                                    "-spiked"),
                                                            exp_name))
                                        else
                                            return(exp_name)
                                    })))
            
            ## Creating a list:
            ## Each element correspond to a location (pf, mf, ef, af).
            ## Each element contains a list of PlotSetPairs object for spiked
            ## files and, if notScaled = TRUE, BGSub files.
            
            
            plot_list <- lapply(names(list_arrayList[[1]]), function(location){
                        
                        unlist(lapply(list_arrayList, 
                                        function(dataset){
                                            dataset <- dataset[[location]]
                                            exp_nb <- dataset$ntracks()
                                            ind <- .retrieveIndexes(
                                                    exp_nb, notScaled)
                                            return(unlist(vapply(ind,
                                                function(x) list(dataset[[x]]),
                                                c(new("PlotSetPair")))
                                                    ))
                                        }))
                    })
            
            
            names(plot_list) <- names(list_arrayList[[1]])
            
            type_names <- .convertTypes(plot_list)
            
            type_count_vec <- seq_len(length(plot_list))
            
            if(!legends)
                par(mfrow=c(2,2))
            else
                par(mfrow=c(3,2))
            
            invisible(mapply(function(location_list, type_count){
                                
                                plotAverage(location_list,
                                        keepratio = FALSE,
                                        type = "full",
                                        cex.lab = 9,
                                        main = type_names[type_count],
                                        legend=FALSE,
                                        colvec = if(!is.null(colVec)) colVec
                                                else NULL,
                                        cex.main = 10,
                                        cex.axis = 10,
                                        las = 1)
                            }, plot_list, type_count_vec))
            
            .plotLegends(legends, colVec, length(expNamesVec), 
                    expNamesVec)
        }
)


setMethod(
        
        f = "plotProfile",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, legends = FALSE, colVec = NULL){
            
            list_arrayList <- lapply(getDatasetList(theObject), 
                    function(object){
                        return(getAverageBindingValues(object))
                    })
            
            expNamesVec <- as.character(unlist(lapply(
                                    getDatasetList(theObject), getExpName)))
            
            ## Creating a list:
            ## Each element correspond to a location (pf, mf, ef, af).
            ## Each element contains a list of PlotSetPairs object for spiked
            ## files and, if notScaled = TRUE, BGSub files.
            
            plot_list <- lapply(names(list_arrayList[[1]]), function(location){
                        
                        unlist(lapply(list_arrayList, 
                                        function(dataset){
                                            dataset <- dataset[[location]]
                                            exp_nb <- dataset$ntracks()
                                            return(unlist(
                                                        vapply(seq_len(exp_nb),
                                                                    function(x)
                                                           list(dataset[[x]]),
                                                       c(new("PlotSetPair"))
                                                            )))
                                        }))})
            
            
            names(plot_list) <- names(list_arrayList[[1]])
            type_names <- .convertTypes(plot_list)
            type_count_vec <- seq_len(length(plot_list))
            
            if(!legends)
                par(mfrow=c(2,2))
            else
                par(mfrow=c(3,2))
            
            invisible(mapply(function(location_list, type_count){
                                
                                plotAverage(location_list,
                                        keepratio = FALSE,
                                        type = "full",
                                        cex.lab = 9,
                                        main = type_names[type_count],
                                        legend=FALSE,
                                        colvec = if(!is.null(colVec)) colVec
                                                else NULL,
                                        cex.main = 10,
                                        cex.axis = 10,
                                        las = 1)
                            }, plot_list, type_count_vec))
            
            .plotLegends(legends, colVec, length(expNamesVec), 
                    expNamesVec)
        }
)


################
## Plotting transformation steps
################


setMethod(
        
        f = "plotTransform",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, legends = FALSE, colVec = NULL, 
                separateWindows = FALSE){
            
            arrayList <- getAverageBindingValues(theObject)
            experiment_number <- arrayList[[1]]$ntracks()/5
            
            list_indexes <- split(seq_len(arrayList[[1]]$ntracks()), 
                    rep(seq_len(experiment_number), each = 5))
            
            expNamesVec <- unlist(lapply(getExperimentList(theObject), 
                            getExpName), use.names = FALSE)
            
            lt <- c(pf="start", mf="midpoint", ef="end", af="composite")
            type_location <- lt[names(arrayList)]
            
            type_transform <- c("noSpike", "RPM", "RPM_BGSub", 
                    "RPM_BGSub_reverted", "spiked")
            
            exp_count_vec <- seq_len(length(list_indexes))
            
            invisible(mapply(function(indexes_track, name_track, exp_count, 
                                    locations, transform, binding_list, 
                                    nb_exp){
                        
                        .plotAverageProfiles(legends, binding_list, 
                                indexes_track, paste(name_track, locations, 
                                        sep="-"), colVec)
                        
                        .plotLegends(legends, colVec, experiment_number, 
                                transform)
                        
                        if(separateWindows){
                            
                            if(exp_count < experiment_number){
                                
                                dev.new()
                            }
                        }
                   }, list_indexes, expNamesVec, exp_count_vec, MoreArgs = 
                           list(type_location, type_transform, arrayList, 
                                   experiment_number)))
        }
)


setMethod(
        
        f = "plotTransform",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, legends = FALSE, colVec = NULL, 
                separateWindows = FALSE){
            
            nb_element <- length(getDatasetList(theObject))
            count_vec <- seq_len(length(getDatasetList(theObject)))
            
            invisible(mapply(function(dataset, count){
                                
                                plotTransform(dataset, legends, colVec, 
                                        separateWindows)
                                
                                if(separateWindows)
                                    if(count < nb_element)
                                        dev.new()
                                }, getDatasetList(theObject), count_vec))
        }
)


################
## Plotting heatmaps
################


.checkPlotHeatmapsParam <- function(location, transformType = NULL){
    
    if(!(isTRUE(all.equal(location, "start")) ||
                isTRUE(all.equal(location, "end")) ||
                isTRUE(all.equal(location, "midpoint")) ||
                isTRUE(all.equal(location, "composite"))))
        stop("location should be start, end, midpoint or composite.")
    
    if(!is.null(transformType)){
        
        if(!(isTRUE(all.equal(transformType, "spiked")) ||
                    isTRUE(all.equal(transformType, "reverse")) ||
                    isTRUE(all.equal(transformType, "BGSub")) ||
                    isTRUE(all.equal(transformType, "RPM")) ||
                    isTRUE(all.equal(transformType, "raw"))))
            stop("transformType should be raw, RPM, BGSub, reverse or ",
                    "spiked.")
    }
}



.retrieveIndexesFromAllTypes <- function(transformType, experiment_number){
    
    result <- if(experiment_number > 5)
                switch(transformType,
                   spiked = seq(from = 5, to = (5*experiment_number), by = 5),
                   reverse = seq(from = 4, to = (5*experiment_number), by = 5),
                   BGSub = seq(from = 3, to = (5*experiment_number), by = 5),
                   RPM = seq(from = 2, to = (5*experiment_number), by = 5),
                   raw = seq(from = 1, to = (5*experiment_number), by = 5))
            else
                switch(transformType,
                        spiked = 5,
                        reverse = 4,
                        BGSub = 3,
                        RPM = 2,
                        raw = 1)
    return(result)
}


setMethod(
        
        f = "plotHeatmaps",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, location = "start", 
                transformType = "spiked", legend = TRUE, plot_scale = "no",
                sort_rows = "decreasing", nb_of_groups = 1, 
                clustering_method = "none", include_exp_vec = NULL, 
                auto_scale = FALSE, raster_value = FALSE, col_value = "blue",
                ...){
            
            .checkPlotHeatmapsParam(location, transformType)
            
            location <- switch(location, start="pf", midpoint="mf", end="ef",
                    composite="af")
            arrayList <- getAverageBindingValues(theObject)
            experiment_number <- arrayList[[1]]$ntracks()/5
            index_vec <- .retrieveIndexesFromAllTypes(transformType)
            expNames_vec <- unlist(lapply(getExperimentList(theObject), 
                            getExpName), use.names = FALSE)
            expNames_vec <- paste0(expNames_vec,
                            switch(transformType, 
                                    spiked = "-spiked",
                                    reverse = "-rev",
                                    BGSub = "-BGSub",
                                    RPM = "-RPM",
                                    raw = "-raw"))
            
            ## Only the first experiment is used for clustering
            if(is.null(include_exp_vec) && experiment_number > 1)
                include_exp_vec <- c(TRUE, rep(FALSE, experiment_number - 1))
                    
            plotHeatmap(arrayList[[location]][index_vec], 
                    labels = expNames_vec,
                    legend = legend,
                    plotScale = plot_scale, 
                    sortrows = sort_rows,
                    clusters = nb_of_groups, 
                    clstmethod = clustering_method, 
                    include = include_exp_vec, 
                    autoscale = auto_scale, 
                    colvec = if(length(col_value) == 1) 
                                rep(col_value, experiment_number) else 
                                col_value, raster = raster_value, ...)
        }
)



setMethod(
        
        f = "plotHeatmaps",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, location = "start", transformType = 
                        "spiked", legend = TRUE, plot_scale = "no", sort_rows =
                        "decreasing", nb_of_groups = 1, clustering_method = 
                        "none", include_exp_vec = NULL, auto_scale = FALSE, 
                raster_value = FALSE, col_value = "blue", ...){
            
            if(!isTRUE(all.equal(transformType, "spiked")))
                stop("Raw, RPM scaled, background subtracted or RPM reverted ",
                        "experiments are not available in boost mode.")
            
            .checkPlotHeatmapsParam(location)
            
            location <- switch(location, start="pf", midpoint="mf", end="ef",
                    composite="af")
            arrayList <- getAverageBindingValues(theObject)
            experiment_number <- arrayList[[1]]$ntracks()
            
            expNames_vec <- unlist(lapply(getExperimentList(theObject), 
                            getExpName), use.names = FALSE)
            
            ## Only the first experiment is used for clustering
            if(is.null(include_exp_vec) && experiment_number > 1)
                include_exp_vec <- c(TRUE, rep(FALSE, experiment_number - 1))
            
            plotHeatmap(arrayList[[location]], 
                    labels = expNames_vec,
                    legend = legend,
                    plotScale = plot_scale, 
                    sortrows = sort_rows,
                    clusters = nb_of_groups, 
                    clstmethod = clustering_method, 
                    include = include_exp_vec, 
                    autoscale = auto_scale,  
                    colvec = if(length(col_value) == 1) 
                                rep(col_value, experiment_number) else 
                                col_value, raster = raster_value, ...)
        }
)


setMethod(
        
        f = "plotHeatmaps",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, location = "start", 
                transformType = "spiked", legend = TRUE, plot_scale = "no",
                sort_rows = "decreasing", nb_of_groups = 1, 
                clustering_method = "none", include_exp_vec = NULL, 
                auto_scale = FALSE, raster_value = FALSE,
                col_value = "blue", ...){
            
            .checkPlotHeatmapsParam(location, transformType)
            
            list_arrayList <- lapply(getDatasetList(theObject), 
                    function(object){
                        return(getAverageBindingValues(object))
                    })
            
            location <- switch(location, start="pf", midpoint="mf", end="ef",
                    composite="af")
            
            
            plot_list <- unlist(lapply(list_arrayList, function(dataset){
                                
                                dataset <- dataset[[location]]
                                exp_nb <- dataset$ntracks()
                                ind <- .retrieveIndexesFromAllTypes(
                                        transformType, exp_nb)
                                return(unlist(vapply(ind,
                                                        function(x)
                                                          list(dataset[[x]])
                                                        ,c(new("PlotSetPair"))
                                )))
                            }))
            
            expNames_vec <- as.character(unlist(lapply(
                                    getDatasetList(theObject), getExpName)))
            
            expNames_vec <- paste0(expNames_vec,
                    switch(transformType, 
                            spiked = "-spiked",
                            reverse = "-rev",
                            BGSub = "-BGSub",
                            RPM = "-RPM",
                            raw = "-raw"))
            
            experiment_number <- length(expNames_vec)
            
            ## Only the first experiment is used for clustering
            if(is.null(include_exp_vec) && experiment_number > 1)
                include_exp_vec <- c(TRUE, rep(FALSE, experiment_number - 1))
            
            plotHeatmap(plot_list, 
                    labels = expNames_vec,
                    legend = legend,
                    plotScale = plot_scale, 
                    sortrows = sort_rows,
                    clusters = nb_of_groups, 
                    clstmethod = clustering_method, 
                    include = include_exp_vec, 
                    autoscale = auto_scale, 
                    colvec = if(length(col_value) == 1) 
                                rep(col_value, experiment_number) else 
                                col_value, raster = raster_value,...)
        }
)


setMethod(
        
        f = "plotHeatmaps",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, location = "start", transformType = 
                        "spiked", legend = TRUE, plot_scale = "no", sort_rows =
                        "decreasing", nb_of_groups = 1, clustering_method = 
                        "none", include_exp_vec = NULL, auto_scale = FALSE, 
                raster_value = FALSE, col_value = "blue", ...){
            
            if(!isTRUE(all.equal(transformType, "spiked")))
                stop("Raw, RPM scaled, background subtracted or RPM reverted ",
                        "experiments are not available in boost mode.")
            
            .checkPlotHeatmapsParam(location)
            
            list_arrayList <- lapply(getDatasetList(theObject), 
                    function(object){
                        return(getAverageBindingValues(object))
                    })
            
            location <- switch(location, start="pf", midpoint="mf", end="ef",
                    composite="af")
            
            plot_list <- unlist(lapply(list_arrayList, function(dataset){
                                
                                dataset <- dataset[[location]]
                                exp_nb <- dataset$ntracks()
                                return(unlist(
                                                vapply(seq_len(exp_nb), 
                                                        function(x) 
                                                           list(dataset[[x]]),
                                                       c(new("PlotSetPair"))
                                )))
                }))
            
            expNames_vec <- as.character(unlist(lapply(
                                    getDatasetList(theObject), getExpName)))
            
            experiment_number <- length(expNames_vec)
            
            ## Only the first experiment is used for clustering
            if(is.null(include_exp_vec) && experiment_number > 1)
                include_exp_vec <- c(TRUE, rep(FALSE, experiment_number - 1))
            
            plotHeatmap(plot_list, 
                    labels = expNames_vec,
                    legend = legend,
                    plotScale = plot_scale, 
                    sortrows = sort_rows,
                    clusters = nb_of_groups, 
                    clstmethod = clustering_method, 
                    include = include_exp_vec, 
                    autoscale = auto_scale, 
                    colvec = if(length(col_value) == 1) 
                                rep(col_value, experiment_number) else 
                                col_value, raster = raster_value, ...)
        }
)


################
## boxplotSpike
################


.perform_boxplot_violinplot <- function(input_matrix, violinPlot, ylab, colvec,
        outline, expname_vec, notch, indicate_mean_with_sd, indicate_mean, 
        indicate_median, indicate_boxplot, indicate_jitter, plot){
    
    if(!is.null(colvec))
        cols <- colvec
    else
        cols <- c("darkblue", "darkgreen", "darkred", "darkmagenta", 
                "darkgray", "darkorange", "darkcyan", "black", 
                rainbow(ncol(input_matrix)-8))
    
    if(!violinPlot)
        boxplot.matrix(input_matrix, ylab=ylab, col=cols, outline=outline, 
                las=2, notch=notch)
    else{
        
        if(!outline){
            
            lower_upper <- apply(input_matrix, MARGIN = 2, function(x){
                        return(boxplot.stats(x)$stats[c(1,5)])})
            
            col_number_vec <- seq_len(ncol(input_matrix))
            
            filtered_list <- lapply(col_number_vec, function(col_number){
                        
                        result <- which(input_matrix[,col_number] > 
                                        lower_upper[1,col_number] & 
                                        input_matrix[,col_number] < 
                                        lower_upper[2,col_number])
                        if(length(result)) 
                            return(input_matrix[result,col_number]) 
                        else return(input_matrix[,col_number])
                    })
            names(filtered_list) <- colnames(input_matrix)
            
            nb_rows_vec <- unlist(lapply(filtered_list, length))
            names_category_vec <- rep(names(nb_rows_vec), nb_rows_vec)
            values_vec <- as.numeric(unlist(filtered_list))
            
        }else{
            
            nb_rows_vec <- rep(nrow(input_matrix), ncol(input_matrix))
            names_category_vec <- rep(colnames(input_matrix), nb_rows_vec)
            values_vec <- as.vector(input_matrix)
            
        }
        
        df_fromMatrix <- data.frame(expnames = factor(names_category_vec, 
                        levels = unique(names_category_vec)), 
                values = values_vec)
        
        g <- ggplot(df_fromMatrix, aes_string("expnames", "values", 
                                color="expnames")) + geom_violin(trim=FALSE) + 
                scale_colour_manual(values = cols) + 
                scale_fill_manual(values = rep("white", length(cols))) + 
                theme_classic() + labs(y= ylab)
        
        
        g <- g + theme(panel.background=element_rect(fill='white', 
                        colour='black'), 
                axis.text.x=element_text(angle=90, size=10), 
                legend.position="none")
        
        if(indicate_mean_with_sd)
        {
            g <- g + stat_summary(fun.data = mean_sdl, geom="pointrange", 
                    shape=23, size=1, color="black")
        }
        
        if(indicate_mean)
        {
            g <- g + stat_summary(fun.y = mean, geom="point", 
                    shape=23, size=1, color="black")
        }
        
        if(indicate_median)
        {
            g <- g + stat_summary(fun.y = median, geom="point", shape=23, 
                    size=1)
        }
        
        
        if(indicate_boxplot)
        {
            g <- g + geom_boxplot(width=0.1, outlier.colour="NA")
        }
        
        
        if(indicate_jitter)
        {
            g <- g + geom_jitter(alpha=0.5, 
                    position= position_jitter(width = 0.1))
        }
        
        if(plot)
            g
        
    }
}


.callPlot <- function(input_mat, verbose, col, violinPlot, ylab, outline, 
        notch, mean_with_sd, mean, median, boxplot, jitter, plot){
    
    if(!is.null(col) && !isTRUE(all.equal(length(col), 
                    ncol(input_mat))))
        stop("Number of colors should equal the number of ",
                "experiments")
    
    .perform_boxplot_violinplot(input_mat, violinPlot, ylab, 
            col, outline, colnames(input_mat), notch, 
            mean_with_sd, mean, median, boxplot, jitter, plot)
}



.retrieveMeanBindingValues <- function(matBindingValues_list, verbose){
    
    ## Retrieving mean binding values 
    if(verbose)
        message("\t Preparing mean value matrix.")
    
    input_mat <- mapply(function(current_mat, current_exp){
                
                if(verbose)
                    message("\t\t", current_exp)
                
                current_mean <- 
                        apply(current_mat, MARGIN = 1, mean)
                
                return(current_mean)
            }, matBindingValues_list, 
            names(matBindingValues_list))
    
    return(input_mat)
}


setMethod(
        
        f = "boxplotSpike",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, col = NULL, rawFile = FALSE, 
                rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
                spiked = TRUE, ylab = NULL, outline = TRUE, violinPlot = FALSE,
                notch = TRUE, mean_with_sd = FALSE, mean = FALSE, 
                median = FALSE, boxplot = FALSE, jitter = FALSE, plot = TRUE,
                verbose = FALSE){
            
            
            input_mat <- buildMeanMatrix(theObject, rawFile, rpmFile, 
                    bgsubFile, revFile, spiked, verbose)
            
            .callPlot(input_mat, verbose, col, violinPlot, ylab, outline, 
                    notch, mean_with_sd, mean, median, boxplot, jitter, plot)
            
            }
)


setMethod(
        
        f = "boxplotSpike",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, col = NULL, rawFile = FALSE, 
                rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
                spiked = TRUE, ylab = NULL, outline = TRUE, violinPlot = FALSE,
                notch = TRUE, mean_with_sd = FALSE, mean = FALSE, 
                median = FALSE, boxplot = FALSE, jitter = FALSE, plot = TRUE,
                verbose = FALSE){
            
            input_mat_list <- lapply(getDatasetList(theObject), 
                    function(dataset){
                        
                        return(buildMeanMatrix(dataset, rawFile, rpmFile, 
                                        bgsubFile, revFile, spiked, verbose,
                                        is_list = TRUE))
                    })
            
            input_mat <- do.call(cbind, input_mat_list)
            
            input_mat <- .removingNABigWig(input_mat, verbose)
            
            .callPlot(input_mat, verbose, col, violinPlot, ylab, outline, 
                    notch, mean_with_sd, mean, median, boxplot, jitter, plot)
            
        }
)


setMethod(
        
        f = "boxplotSpike",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, col = NULL, rawFile = FALSE, 
                rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
                spiked = TRUE,ylab = NULL, outline = TRUE, 
                violinPlot = FALSE, notch = TRUE, mean_with_sd = FALSE, 
                mean = FALSE, median = FALSE, boxplot = FALSE, 
                jitter = FALSE, plot = TRUE, verbose = FALSE){
            
            if(rawFile || rpmFile || bgsubFile || revFile)
                stop("Raw, RPM scaled, background subtracted or RPM reverted ",
                        "experiments are not available in boost mode.")
            
            list_matBindingValues_list <- lapply(getDatasetList(theObject), 
                    function(dataset){
                        
                        return(getMatBindingValues(dataset))
                    })
            
            input_mat_list <- lapply(list_matBindingValues_list, 
                    function(matBindingValues_list)
                    return(.retrieveMeanBindingValues(matBindingValues_list, 
                                    verbose)))
            
            input_mat <- do.call(cbind, input_mat_list)
                
            .callPlot(input_mat, verbose, col, violinPlot, ylab, outline, 
                    notch, mean_with_sd, mean, median, boxplot, jitter, plot)
        }
)


setMethod(
        
        f = "boxplotSpike",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, col = NULL, rawFile = FALSE, 
                rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
                spiked = TRUE,ylab = NULL, outline = TRUE, violinPlot = FALSE, 
                notch = TRUE, mean_with_sd = FALSE, mean = FALSE, 
                median = FALSE, boxplot = FALSE, jitter = FALSE, plot = TRUE, 
                verbose = FALSE){
            
            if(rawFile || rpmFile || bgsubFile || revFile)
                stop("Raw, RPM scaled, background subtracted or RPM reverted ",
                        "experiments are not available in boost mode.")
            
            matBindingValues_list <- getMatBindingValues(theObject)
            
            ## Retrieving mean binding values 
            if(verbose)
                message("\t Preparing mean value matrix.")
            
            input_mat <- 
                    .retrieveMeanBindingValues(matBindingValues_list, verbose)
            
            .callPlot(input_mat, verbose, col, violinPlot, ylab, outline, 
                    notch, mean_with_sd, mean, median, boxplot, jitter, plot)
            
        }
)



################
## Correlation plot
################


.perform_heatscatterplot <- function(input_mat, method_scale, allOnPanel, 
        separateWindows, main, show_cor, method_cor, add_contour, nlevels, 
        color_contour, ...){
    
    if(!(isTRUE(all.equal(method_scale, "none")) || 
                isTRUE(all.equal(method_scale, "log")) ||
                isTRUE(all.equal(method_scale, "asinh")) ||
                isTRUE(all.equal(method_scale, "cuberoot")) ||
                isTRUE(all.equal(method_scale, "zscore"))))
        stop("method_scale should be none, log, asinh, cuberoot or ",
                "zscore")
    
    if(allOnPanel && separateWindows)
        stop("allOnPanel and separateWindows cannot be both TRUE.")
    
    name_vec <- colnames(input_mat)
    mat_comb <- combn(seq_len(ncol(input_mat)), 2)
    
    if(isTRUE(all.equal(method_scale, "asinh")))
        input_mat <- asinh(input_mat)
    else if(isTRUE(all.equal(method_scale, "cuberoot")))
        input_mat <- input_mat^(1/3)
    else if(isTRUE(all.equal(method_scale, "zscore")))
        input_mat <- scale(input_mat)
    
    
    if(allOnPanel){
        
        if(isTRUE(all.equal(ncol(mat_comb)%%2,0)))
            par(mfrow=c(ncol(mat_comb)/2, ncol(mat_comb)/2))
        else
            par(mfrow=c(round(ncol(mat_comb)/2)+1, 
                            round(ncol(mat_comb)/2)-1))
    }
    
    
    count_exp_vec <- seq_len(ncol(mat_comb))
    
    for(count_exp in count_exp_vec){
        
        heatscatter(input_mat[,mat_comb[1,count_exp]], 
                input_mat[,mat_comb[2,count_exp]],
                xlab = switch(method_scale,
                        log = paste0("log(", name_vec[mat_comb[1,count_exp]],
                                ")"),
                        asinh = paste0("asinh(", 
                                name_vec[mat_comb[1,count_exp]], ")"),
                        cuberoot = paste0("cuberoot(", 
                                name_vec[mat_comb[1,count_exp]], ")"),
                        zscore = paste0("zscore(", 
                                name_vec[mat_comb[1,count_exp]], ")"),
                        none = name_vec[mat_comb[1,count_exp]]),
                ylab = switch(method_scale,
                        log = paste0("log(",
                                name_vec[mat_comb[2,count_exp]], ")"),
                        asinh = paste0("asinh(",
                                name_vec[mat_comb[2,count_exp]], ")"),
                        cuberoot = paste0("cuberoot(",
                                name_vec[mat_comb[2,count_exp]], ")"),
                        zscore = paste0("zscore(",
                                name_vec[mat_comb[2,count_exp]], ")"),
                        none = name_vec[mat_comb[2,count_exp]]),
                main = main, cor = show_cor, method = method_cor, 
                add.contour = add_contour, nlevels = nlevels, 
                color.contour = color_contour,
                log = ifelse(isTRUE(all.equal(method_scale, 
                                        "log")), "xy", ""),...)
        
        if(separateWindows)
            if(count_exp < ncol(mat_comb))
                dev.new()
    }
}


.loadCorrelationPlot <- function(heatscatterplot, input_mat, method_scale, 
        allOnPanel, separateWindows, main, show_cor, method_cor, add_contour,
        nlevels, color_contour, method_corrplot, type_corrplot, diag_corrplot,
        ...){
    
    if(heatscatterplot){
        
        .perform_heatscatterplot(input_mat, method_scale, allOnPanel, 
                separateWindows, main, show_cor, method_cor, 
                add_contour, nlevels, color_contour, ...)
    }else{
        
        mat_cor <- cor(input_mat, method = method_cor)
        corrplot(mat_cor, method= method_corrplot, type = type_corrplot,
                title = main, diag = diag_corrplot, ...)
        
        return(invisible(mat_cor))
    }
}


setMethod(
        
        f = "plotCor",
        
        signature = "ChIPSeqSpikeDataset",
        
        
        definition = function(theObject, rawFile = FALSE, 
                rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
                spiked = TRUE, main = "", add_contour = FALSE, 
                method_cor = "spearman", nlevels = 10, color_contour = "black",
                show_cor = TRUE, allOnPanel = TRUE, method_scale = "none",
                method_corrplot = "circle", heatscatterplot = TRUE, 
                type_corrplot = "upper", diag_corrplot = FALSE, 
                separateWindows = FALSE, verbose = FALSE, ...){
            
            input_mat <- buildMeanMatrix(theObject, rawFile, rpmFile, 
                    bgsubFile, revFile, spiked, verbose)
            
            .loadCorrelationPlot(heatscatterplot, input_mat, method_scale, 
                    allOnPanel, separateWindows, main, show_cor, method_cor, 
                    add_contour, nlevels, color_contour, method_corrplot, 
                    type_corrplot, diag_corrplot, ...)
            }
)


setMethod(
        
        f = "plotCor",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        
        definition = function(theObject, rawFile = FALSE, rpmFile = FALSE, 
                bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, main = "", 
                add_contour = FALSE, method_cor = "spearman", nlevels = 10, 
                color_contour = "black", show_cor = TRUE, allOnPanel = TRUE, 
                method_scale = "none", method_corrplot = "circle", 
                heatscatterplot = TRUE, type_corrplot = "upper", 
                diag_corrplot = FALSE, separateWindows = FALSE, 
                verbose = FALSE, ...){
            
            if(rawFile || rpmFile || bgsubFile || revFile)
                stop("Raw, RPM scaled, background subtracted or RPM reverted ",
                        "experiments are not available in boost mode.")
            
            input_mat <- buildMeanMatrix(theObject, verbose)
            
            .loadCorrelationPlot(heatscatterplot, input_mat, method_scale, 
                    allOnPanel, separateWindows, main, show_cor, method_cor, 
                    add_contour, nlevels, color_contour, method_corrplot, 
                    type_corrplot, diag_corrplot, ...)
        }
)


setMethod(
        
        f = "plotCor",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        
        definition = function(theObject, rawFile = FALSE, 
                rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
                spiked = TRUE, main = "", add_contour = FALSE, 
                method_cor = "spearman", nlevels = 10, color_contour = "black",
                show_cor = TRUE, allOnPanel = TRUE, method_scale = "none",
                method_corrplot = "circle", heatscatterplot = TRUE, 
                type_corrplot = "upper", diag_corrplot = FALSE, 
                separateWindows = FALSE, verbose = FALSE, ...){
            
            input_mat_list <- lapply(getDatasetList(theObject), 
                    function(dataset){
                        
                        return(buildMeanMatrix(dataset, rawFile, rpmFile, 
                                        bgsubFile, revFile, spiked, verbose,
                                        is_list = TRUE))
                    })
            
            input_mat <- do.call(cbind, input_mat_list)
            input_mat <- .removingNABigWig(input_mat, verbose)
            
            .loadCorrelationPlot(heatscatterplot, input_mat, method_scale, 
                    allOnPanel, separateWindows, main, show_cor, method_cor, 
                    add_contour, nlevels, color_contour, method_corrplot, 
                    type_corrplot, diag_corrplot, ...)
        }
)


setMethod(
        
        f = "plotCor",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        
        definition = function(theObject, rawFile = FALSE, rpmFile = FALSE, 
                bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, main = "", 
                add_contour = FALSE, method_cor = "spearman", nlevels = 10, 
                color_contour = "black", show_cor = TRUE, allOnPanel = TRUE, 
                method_scale = "none", method_corrplot = "circle", 
                heatscatterplot = TRUE, type_corrplot = "upper", 
                diag_corrplot = FALSE, separateWindows = FALSE, 
                verbose = FALSE, ...){
            
            if(rawFile || rpmFile || bgsubFile || revFile)
                stop("Raw, RPM scaled, background subtracted or RPM reverted ",
                        "experiments are not available in boost mode.")
            
            input_mat_list <- lapply(getDatasetList(theObject), 
                    function(dataset){
                        
                        return(buildMeanMatrix(dataset, verbose, 
                                        is_list = TRUE))
                    })
            
            input_mat <- do.call(cbind, input_mat_list)
            input_mat <- .removingNABigWig(input_mat, verbose)
            
            .loadCorrelationPlot(heatscatterplot, input_mat, method_scale, 
                    allOnPanel, separateWindows, main, show_cor, method_cor, 
                    add_contour, nlevels, color_contour, method_corrplot, 
                    type_corrplot, diag_corrplot, ...)
        }
)
