## Preparing data for testing


indir <- system.file("extdata", package="airway", mustWork=TRUE);

file_vec <- system.file(
        c("extdata/empty_files/endoBam_empty.bam", 
                "extdata/empty_files/exoBam_empty.bam",
                "extdata/empty_files/wigfile_empty.bw", 
                "extdata/empty_files/wigfile_empty-RPM.bw"), 
        package="ChIPSeqSpike");

csds <- ChIPSeqSpikeDataset(
        endogenousBam_vec = c(file_vec[1], file_vec[1]), 
        exogenousBam_vec = c(file_vec[2], file_vec[2]), 
        bigWigFile_endogenous_vec = c(file_vec[3], file_vec[3]), 
        inputBigWigFile = file_vec[3], 
        inputBamFile = file_vec[1], 
        expnames = c("toto", "toto2"),
        inputSF = 0.7, inputNb = 1000);

csds2 <- ChIPSeqSpikeDataset(
        endogenousBam_vec = c(file_vec[1], file_vec[1]), 
        exogenousBam_vec = c(file_vec[2], file_vec[2]), 
        bigWigFile_endogenous_vec = c(file_vec[4], file_vec[3]), 
        inputBigWigFile = file_vec[3], 
        inputBamFile = file_vec[1], 
        expnames = c("toto", "toto2"),
        inputSF = 0.7, inputNb = 1000);

e1 <- Experiment(endogenousBamFilePath = file_vec[1], 
        exogenousBamFilePath = file_vec[2], 
        bigWigFilePath = file_vec[3], 
        name = "toto",
        endoScalingFactor = 1, exoScalingFactor = 2, 
        endoNb = 3, exoNb = 4);

e2 <- Experiment(endogenousBamFilePath = paste(indir, 
                list.files(indir)[7], sep="/"), 
        exogenousBamFilePath = file_vec[2], 
        bigWigFilePath = file_vec[3], 
        name = "toto",
        endoScalingFactor = 1, exoScalingFactor = 2, 
        endoNb = 3, exoNb = 4);


## Testing estimateScalingFactors and estimateScalingFactors_Experiment

context("Testing estimateScalingFactors");

test_that("Lack of aligned reads in bam files is handled", {
            
            expect_error(estimateScalingFactors(csds),
                    "input bam file contains no aligned reads.");
            
            expect_error(estimateScalingFactors(e1),
                    "Endogenous bam file contains no aligned reads.");
            
            expect_error(estimateScalingFactors(e2, paired=T),
                    "Exogenous bam file contains no aligned reads.");
            
        });


## Testing scaling and .computeScaling

context("Testing scaling and .computeScaling");

test_that("All parameters are valid", {
            
            expect_error(scaling(csds, outputFolder = "/toto/"), 
                    "The specified output folder does not exist.");
            
            expect_error(scaling(csds, outputFolder = "./"), 
                    "The path to the output folder should not end by '/'");
            
            expect_error(scaling(csds, type = "toto"), 
                    "Accepted types are endo and exo.");
            
            expect_error(scaling(csds, reverse = TRUE, type = "exo"), 
                    "Exogenous scaling factor cannot be reverted.");
            
            expect_error(scaling(csds, reverse = TRUE), 
                    "RPM normalization must be performed before reverting it");
            
            expect_error(scaling(csds2, reverse = TRUE), 
                    paste0("Input subtraction should be performed before ", 
                            "reverting RPM normalization"));
            
            expect_error(scaling(csds, type = "exo"),
                    paste0("Exogenous scaling factor should be applied when ",
                            "RPM normalization has been reverted."));
        });


## Testing inputSubtraction

context("Testing inputSubtraction");

test_that("Invalid data and differing number of chromosomes is not permitted", {
            
            expect_error(inputSubtraction(csds), paste0("RPM normalization must"
                                    ," be performed before subtracting the ",
                                    "input"));
        });

