.retrieveAverageValuesList <- function(files_vec, gff_vec, genome, binsize, 
        before, after, mean_or_median, interpolation_average, verbose){
    
    location_type_vec <- c("pf", "mf", "ef", "af");
    setArray_list <- list();
    
    invisible(sapply(location_type_vec, function(location_name, files, gff_vec,
                    genome_version, bin_size, profile_length_before, 
                    profile_length_after, mean_or_median){
                
                if(verbose)
                    message("Extraction at ", location_name);
                
                setArray_list[[location_name]] <<- 
                        getPlotSetArray(tracks = files, 
                                features = gff_vec, 
                                refgenome = genome_version, 
                                bin = bin_size, 
                                xmin = profile_length_before, 
                                xmax = profile_length_after, 
                                xanchored = interpolation_average,
                                type = location_name, 
                                add_heatmap = TRUE, 
                                stat = mean_or_median);
              
            }, files_vec, gff_vec, genome, binsize, before, after, 
            mean_or_median));
    
    return(setArray_list);
};


.retrieveMatBindingList <- function(files_vec, gff_vec, interpolation_number, 
        ignore_strand, verbose){
    
    current_gff <- read.table(gff_vec, stringsAsFactors=FALSE);
    current_gff_GRanges <- GRanges(seqnames=current_gff[,1], 
            ranges=IRanges(start= current_gff[,4], 
                    end = current_gff[,5], 
                    names = current_gff[,2]), 
            strand=current_gff[,7]);
    
    mat_list <- list();
    expname_vec <- names(files_vec);
    
    invisible(mapply(function(file_path, expname){
                
                if(verbose)
                    message("\t Processing ", expname);
                
                mat_list[[expname]] <<- summary(BigWigFile(file_path), 
                                                which = current_gff_GRanges, 
                                                size = interpolation_number, 
                                                type = "mean", as = "matrix");
                
                if (!ignore_strand)
                    mat_list[[expname]][
                            as.character(strand(current_gff_GRanges))=='-',] <-
                            mat_list[[expname]][
                                    as.character(strand(current_gff_GRanges))==
                                            '-', ncol(mat_list[[expname]]):1];
                
            }, files_vec, expname_vec));
    
    return(mat_list);
};


.notFoundMessage <- function(file_name){
    
    if(!file.exists(file_name))
        stop("The file ", file_name, " was not found. Did you perform all ",
                "steps of scaling?");
};


.retrieveAllFiles <- function(experiment_list){
    
    files_vec <- unlist(unname(lapply(experiment_list, 
                            function(exp){
                                
                                if(!length(grep("spiked", getBigWigFile(exp))))
                                    stop("Data must be spiked for profiling");
                                
                                exp_name <- getExpName(exp);
                                spiked_bigWig <- getBigWigFile(exp);
                                noSpike_bigWig <-
                                        strsplit(spiked_bigWig, 
                                            "-RPM-BGSub-reverted-spiked")[[1]];
                                
                                RPM_bigWig <- paste0(noSpike_bigWig[1], "-RPM",
                                        noSpike_bigWig[2]);
                                RPM_BGSub_bigWig <- paste0(noSpike_bigWig[1], 
                                        "-RPM-BGSub", noSpike_bigWig[2]);
                                RPM_BGSub_reverted_bigWig <- 
                                        paste0(noSpike_bigWig[1], 
                                        "-RPM-BGSub-reverted",
                                        noSpike_bigWig[2]);
                                noSpike_bigWig <- paste0(noSpike_bigWig[1], 
                                        noSpike_bigWig[2]);
                                
                                .notFoundMessage(spiked_bigWig);
                                .notFoundMessage(noSpike_bigWig);
                                .notFoundMessage(RPM_bigWig);
                                .notFoundMessage(RPM_BGSub_bigWig);
                                .notFoundMessage(RPM_BGSub_reverted_bigWig);
                                
                                name_vec <- paste(exp_name, c("-raw", "-RPM", 
                                                "-BGSub", "-rev", "-spike"), 
                                        sep="");
                                result_files <- c(noSpike_bigWig, RPM_bigWig, 
                                        RPM_BGSub_bigWig, 
                                        RPM_BGSub_reverted_bigWig, 
                                        spiked_bigWig);
                                names(result_files) <- name_vec;
                                
                                return(result_files);
                            })));
    
    return(files_vec);
};

setMethod(
        
        f = "extractBinding",
        
        signature = "ChIPSeqSpikeDataset",
        
        definition = function(theObject, gff_vec, genome, binsize = 50, 
                before = 2000, after=2000, mean_or_median = "mean", 
                interpolation_number = 100, interpolation_average = 10000,
                ignore_strand = FALSE, verbose = FALSE){
            
            ## Getting all files with different transformation applied to
            
            files_vec <- .retrieveAllFiles(theObject@experimentList);
     
     if(verbose)
         message("Retrieving values for profiles and heatmaps.");
     setArray_list <- .retrieveAverageValuesList(files_vec, gff_vec, genome, 
             binsize, before, after, mean_or_median, interpolation_average,
             verbose);
     
     if(verbose)
         message("Retrieving interpolated matrices:");
     matBinding_list <- .retrieveMatBindingList(files_vec, gff_vec, 
             interpolation_number, ignore_strand, verbose);
     
     averageBindingValues(theObject) <- setArray_list;
     matBindingValues(theObject) <- matBinding_list;
     
     return(theObject);
            
        });


setMethod(
        
        f = "extractBinding",
        
        signature = "ChIPSeqSpikeDatasetBoost",
        
        definition = function(theObject, gff_vec, genome, binsize = 50, 
                before = 2000, after=2000, mean_or_median = "mean", 
                interpolation_number = 100, interpolation_average = 10000,
                ignore_strand = FALSE, verbose = FALSE){
            
            ## Getting all files with different transformation applied to
            
            files_vec <- unlist(lapply(theObject@experimentListLoaded, 
                            function(exp){return(getBigWigFile(exp));}));
            
            if(!unique(sapply(files_vec, file.exists, simplify = TRUE)))
                stop("Use exportBigWigs function before extracting values.");
            
            if(verbose)
                message("Retrieving values for profiles and heatmaps.");
            setArray_list <- .retrieveAverageValuesList(files_vec, gff_vec, 
                    genome, binsize, before, after, mean_or_median, 
                    interpolation_average, verbose);
            
            if(verbose)
                message("Retrieving interpolated matrices:");
            matBinding_list <- .retrieveMatBindingList(files_vec, gff_vec, 
                    interpolation_number, ignore_strand, verbose);
            
            averageBindingValues(theObject) <- setArray_list;
            matBindingValues(theObject) <- matBinding_list;
            
            return(theObject);
            
        });


.loadExtractionOnList <- function(theObject, gff_vec, genome, binsize, before,
        after, mean_or_median, interpolation_number, interpolation_average,
        ignore_strand, verbose){
    
    theObject@datasetList <- lapply(theObject@datasetList, 
            function(object){
                
                return(extractBinding(object, gff_vec, genome, binsize,
                                before, after, mean_or_median, 
                                interpolation_number, 
                                interpolation_average,
                                ignore_strand, verbose));
            });
    
    return(theObject);
};


setMethod(
        
        f = "extractBinding",
        
        signature = "ChIPSeqSpikeDatasetList",
        
        definition = function(theObject, gff_vec, genome, binsize = 50, 
                before = 2000, after=2000, mean_or_median = "mean", 
                interpolation_number = 100, interpolation_average = 10000,
                ignore_strand = FALSE, verbose = FALSE){
            
            .loadExtractionOnList(theObject, gff_vec, genome, binsize, before,
                    after, mean_or_median, interpolation_number, 
                    interpolation_average, ignore_strand, verbose);
            }
);


setMethod(
        
        f = "extractBinding",
        
        signature = "ChIPSeqSpikeDatasetListBoost",
        
        definition = function(theObject, gff_vec, genome, binsize = 50, 
                before = 2000, after=2000, mean_or_median = "mean", 
                interpolation_number = 100, interpolation_average = 10000,
                ignore_strand = FALSE, verbose = FALSE){
            
            .loadExtractionOnList(theObject, gff_vec, genome, binsize, before,
                    after, mean_or_median, interpolation_number, 
                    interpolation_average, ignore_strand, verbose);
        }
);
