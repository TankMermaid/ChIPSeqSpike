spikePipe <- function(infoFile, bamPath, bigWigPath, anno, genome_version, 
                paired = FALSE, binsize = 50, profile_length_before = 2000, 
                profile_length_after= 2000, mean_or_median = "mean", 
                interpolation_number = 100, interpolation_average = 10000,
                ignore_strand = FALSE, verbose = FALSE, boost = FALSE, 
                outputFolder = NULL){
            
            
            csds <- spikeDataset(infoFile, bamPath, bigWigPath, boost, 
                    verbose);
            
            if(verbose)
                message("\n\n\t\t ### Step 1. Computing scaling factors ###");
            
            csds <- estimateScalingFactors(csds, paired, verbose);
            
            if(verbose)
                message("\n\n\t\t ### Step 2. RPM scaling ###");
            
            csds <- scaling(csds, verbose = verbose, 
                            outputFolder = outputFolder);
            
            if(verbose)
                message("\n\n\t\t ### Step 3. Input Subtraction ###");
            
            csds <- inputSubtraction(csds, verbose);
            
            if(verbose)
                message("\n\n\t\t ### Step 4. Reverse RPM scaling ###");
            
            csds <- scaling(csds, reverse = TRUE, verbose = verbose);
            
            if(verbose)
                message("\n\n\t\t ### Step 5. Spike-in scaling ###");
            
            csds <- scaling(csds, type = "exo", verbose = verbose);
            
            if(boost){
                if(verbose)
                    message("\n\n\t\t ### Step 6. Writing spiked files ###");
                
                exportBigWigs(csds, verbose);
            };
            
            if(boost){
                if(verbose)
                    message("\n\n\t\t ### Step 7. Extract values ###");
            }else{
                if(verbose)
                    message("\n\n\t\t ### Step 6. Extract values ###");
            };
            
            csds <- extractBinding(csds, anno, genome_version, binsize, 
                    profile_length_before, profile_length_after, 
                    mean_or_median, interpolation_number, 
                    interpolation_average, ignore_strand, verbose);
            
            return(csds);
                
        };
