## This file verifies that classes declarations throw expected messages.

## Preparing data

file_vec <- system.file(c("extdata/empty_files/endoBam_empty.bam",
                "extdata/empty_files/exoBam_empty.bam",
                "extdata/empty_files/input.wig_empty.bw",
                "extdata/empty_files/inputwig_empty.wig",
                "extdata/empty_files/inputwig_empty.bw"), 
        package="ChIPSeqSpike")


## ------------------- ChIPSeqSpikeDataset object -------------------

context("Testing the creation of ChIPSeqSpikeDataset object")

test_that("Error when elements have different lengths.", {
            
            expect_error(ChIPSeqSpikeDataset(endogenousBam_vec = c(1,2), 
                            exogenousBam_vec = c(1), 
                            bigWigFile_endogenous_vec = c(1,3), 
                            inputBigWigFile = c(1),
                            inputBamFile = c(1),
                            expnames = c("toto", "toto")), 
"A corresponding exogenous bam is needed for each experiment.")

expect_error(ChIPSeqSpikeDataset(endogenousBam_vec = c(1,2), 
                exogenousBam_vec = c(2,4), 
                bigWigFile_endogenous_vec = c(1), 
                inputBigWigFile = c(1),
                inputBamFile = c(1),
                expnames = c("toto", "toto")), 
        "A bigWig and bam file should be given for each experiment.")

            expect_error(ChIPSeqSpikeDataset(endogenousBam_vec = c(1,2), 
                            exogenousBam_vec = c(1,2), 
                            bigWigFile_endogenous_vec = c(1,2), 
                            inputBigWigFile = c(1),
                            inputBamFile = c(1),
                            expnames = c("toto")), 
        "One name should be provided per experiment.")

            suppressWarnings(expect_error(ChIPSeqSpikeDataset(
                                    endogenousBam_vec = c(1,2),
                                    exogenousBam_vec = c(1,2),
                                    bigWigFile_endogenous_vec = c(1,2),
                                    inputBigWigFile = c(1,2),
                                    inputBamFile = c(1),
                                    expnames = c("toto", "toto")),
                            paste0("Only one bigWig and bam file should be ",
                                    "given for the input experiment.")))

            suppressWarnings(expect_error(ChIPSeqSpikeDataset(
                                    endogenousBam_vec = c(1,2),
                                    exogenousBam_vec = c(1,2),
                                    bigWigFile_endogenous_vec = c(1,2),
                                    inputBigWigFile = c(1),
                                    inputBamFile = c(1,2),
                                    expnames = c("toto", "toto")),
                            paste0("Only one bigWig and bam file should be ",
                                    "given for the input experiment.")))
}) 


test_that("Invalid paths for the input bam and wig files are not accepted", {
            
          expect_error(ChIPSeqSpikeDataset(
                          endogenousBam_vec = c(file_vec[1], file_vec[1]), 
                          exogenousBam_vec = c(file_vec[2], file_vec[2]), 
                          bigWigFile_endogenous_vec = c(file_vec[5], file_vec[5]), 
                          inputBigWigFile = file_vec[3], 
                          inputBamFile = file_vec[1], 
                          expnames = c("toto", "toto2"),
                          inputSF = 0.7, inputNb = 1000), 
                  paste0(file_vec[3], 
                          " should contain only one point for the extension."))
            
          expect_error(ChIPSeqSpikeDataset(
                          endogenousBam_vec = c(file_vec[1], file_vec[1]), 
                          exogenousBam_vec = c(file_vec[2], file_vec[2]), 
                          bigWigFile_endogenous_vec = c(file_vec[5], file_vec[5]), 
                          inputBigWigFile = "toto", 
                          inputBamFile = file_vec[1], 
                          expnames = c("toto", "toto2"),
                          inputSF = 0.7, inputNb = 1000), 
                  "toto is not a valid path")
          
          expect_error(ChIPSeqSpikeDataset(
                          endogenousBam_vec = c(file_vec[1], file_vec[1]), 
                          exogenousBam_vec = c(file_vec[2], file_vec[2]), 
                          bigWigFile_endogenous_vec = c(file_vec[5], file_vec[5]), 
                          inputBigWigFile = file_vec[4], 
                          inputBamFile = file_vec[1], 
                          expnames = c("toto", "toto2"),
                          inputSF = 0.7, inputNb = 1000),
                  "Wig files should be in bigWig format.")
          
          expect_error(ChIPSeqSpikeDataset(
                          endogenousBam_vec = c(file_vec[1], file_vec[1]), 
                          exogenousBam_vec = c(file_vec[2], file_vec[2]), 
                          bigWigFile_endogenous_vec = c(file_vec[4], file_vec[5]), 
                          inputBigWigFile = file_vec[5], 
                          inputBamFile = file_vec[1], 
                          expnames = c("toto", "toto2"),
                          inputSF = 0.7, inputNb = 1000),
                  "Wig files should be in bigWig format.")
          
        })

test_that("PlotSetArrayList is of the right type", {
            
            expect_error(ChIPSeqSpikeDataset(
                            endogenousBam_vec = c(file_vec[1], file_vec[1]), 
                            exogenousBam_vec = c(file_vec[2], file_vec[2]), 
                           bigWigFile_endogenous_vec = c(file_vec[5], file_vec[5]), 
                            inputBigWigFile = file_vec[5], 
                            inputBamFile = file_vec[1], 
                            expnames = c("toto", "toto2"),
                            inputSF = 0.7, inputNb = 1000, SetArrayList = 3),
                    "SetArrayList should be of type list.")
            
            expect_error(ChIPSeqSpikeDataset(
                   endogenousBam_vec = c(file_vec[1], file_vec[1]), 
                   exogenousBam_vec = c(file_vec[2], file_vec[2]), 
                   bigWigFile_endogenous_vec = c(file_vec[5], file_vec[5]), 
                   inputBigWigFile = file_vec[5], 
                   inputBamFile = file_vec[1], 
                   expnames = c("toto", "toto2"),
                   inputSF = 0.7, inputNb = 1000, SetArrayList = list(c(3,2))),
                  "All objects in setArray list must be of type PlotSetArray.")
        })
        
        
        
## ------------------- Experiment object ------------------- 

context("Testing the creation of Experiment object")

test_that("Invalid paths for bam and wig files are not accepted",{
            
            expect_error(
                    Experiment(endogenousBamFilePath="toto", 
                            exogenousBamFilePath="toto2", 
                            bigWigFilePath="toto3", name="a"),
                    "toto is not a valid path")
            
            expect_error(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="toto2", 
                            bigWigFilePath="toto3", name="a"),
                    "toto2 is not a valid path")
            
            expect_error(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="toto3", name="a"),
                    "toto3 is not a valid path")
            
            expect_error(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFile= file_vec[3], name="a"),
                    paste0(file_vec[3], 
                          " should contain only one point for the extension."))
          
        })


test_that("Negative scaling factors and counts are not accepted",{
            
            expect_error(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath= file_vec[5],name="a",
                            endoScalingFactor = -2),
            "object@endogenousScalingFactor should be a positive number or 0.")
            
            expect_error(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                    bigWigFilePath= file_vec[5], name="a",
                    exoScalingFactor = -2),
            "object@exogenousScalingFactor should be a positive number or 0.")
    
            expect_error(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                    bigWigFile= file_vec[5], name="a",
                    endoNb = -2),
            "object@endoCount should be a positive number or 0.")
    
            expect_error(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                    bigWigFilePath= file_vec[5], name="a",
                    exoNb = -2),
            "object@exoCount should be a positive number or 0.")
    
        })


test_that("If files are not strings, not valid",{
            
            expect_warning(
                    Experiment(endogenousBamFilePath=0, 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name="a"),
                    "endogenousBamFilePath should be a string.")
            
            expect_warning(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath=0, 
                            bigWigFilePath="./", name="a"),
                    "exogenousBamFilePath should be a string.")
            
            expect_warning(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath=0, name="a"),
                    "bigWigFilePath should be a string.")
            
            expect_equal(
                    Experiment(endogenousBamFilePath=0, 
                            exogenousBamFilePath="./", 
                            bigWigFile="./", name="a"),
                    NA)
            
            expect_equal(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath=0, 
                            bigWigFile="./", name="a"),
                    NA)
            
            expect_equal(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFile=0, name="a"),
                    NA)
        })


test_that("If scaling factors or counts are not numeric, not valid",{
            
            expect_warning(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name="a", 
                            endoScalingFactor = "2"),
                    "Scaling factors and counts should be numeric.")
            
            expect_warning(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFile="./", name="a",
                            exoScalingFactor = "2"),
                    "Scaling factors and counts should be numeric.")
            
            expect_warning(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name="a",
                            endoNb = "2"),
                    "Scaling factors and counts should be numeric.")
            
            expect_warning(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name="a",
                            exoNb = "2"),
                    "Scaling factors and counts should be numeric.")
            
            expect_equal(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name="a",
                            endoScalingFactor = "2"),
                    NA)
            
            expect_equal(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFile="./", name="a",
                            exoScalingFactor = "2"),
                    NA)
            
            expect_equal(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name="a",
                            endoNb = "2"),
                    NA)
            
            expect_equal(
                    Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name="a",
                            exoNb = "2"),
                    NA)
            
        })


test_that("Experiment should be an alphanumeric", {
            
            expect_warning(Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name= matrix()),
                    "name should be an alphanumeric.")
            
            expect_equal(Experiment(endogenousBamFilePath="./", 
                            exogenousBamFilePath="./", 
                            bigWigFilePath="./", name= matrix()),
                    NA)
            
        })
