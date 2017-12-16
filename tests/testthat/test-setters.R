## Test that changing values of slots by setters is robust

file_vec <- system.file("extdata", 
        c("bam_files/H3K79me2_0_dm3-filtered.bam", 
                "bam_files/H3K79me2_0_hg19-filtered.bam", 
                "bigwig_files/H3K79me2_0-filtered.bw", 
                "bigwig_files/input_0-filtered.bw", 
                "bam_files/input_0_hg19-filtered.bam"),
        package="ChIPSeqSpike")

file_empty <- system.file(c("extdata/empty_files/input.wig_empty.bw"),
        package="ChIPSeqSpike")

e <- Experiment(endogenousBamFilePath = file_vec[2], 
        exogenousBamFilePath = file_vec[1], 
        bigWigFilePath = file_vec[3], 
        name = "H3K79me2_0", endoScalingFactor = 1, exoScalingFactor = 2, 
        endoNb = 3, exoNb = 4)

eboost <- ExperimentLoaded(endogenousBamFilePath = file_vec[2], 
        exogenousBamFilePath = file_vec[1], 
        bigWigFilePath = file_vec[3], 
        name = "H3K79me2_0", endoScalingFactor = 1, exoScalingFactor = 2, 
        endoNb = 3, exoNb = 4)

csds <- ChIPSeqSpikeDataset(endogenousBam_vec = file_vec[2],
        exogenousBam_vec = file_vec[1],
        bigWigFile_endogenous_vec = file_vec[3], 
        inputBigWigFile = file_vec[4], 
        inputBamFile = file_vec[5], 
        expnames = "H3K79me2_0", inputSF = 0.7, inputNb = 1000)

csdsboost <- ChIPSeqSpikeDatasetBoost(endogenousBam_vec = file_vec[2],
        exogenousBam_vec = file_vec[1],
        bigWigFile_endogenous_vec = file_vec[3], 
        inputBigWigFile = file_vec[4], 
        inputBamFile = file_vec[5], 
        expnames = "H3K79me2_0", inputSF = 0.7, inputNb = 1000)


## For ChIPSeqSpikeDataset object

context("Scaling factors of ChIPSeqSpikeDataset objects are changed correctly.")

test_that("Negative value cannot be given as scaling factor or count", {

           expect_error(scalingFactor(csds) <- -2, 
               "object@inputScalingFactor should be a positive number or 0.")
            
           expect_error(count(csds) <- -2,
                   "object@inputCount should be a positive number or 0.")
            
        })


test_that("Non numeric value is not accepted as scaling factor or count",{
            
            expect_error(scalingFactor(csds) <- "toto", 
                    "Scaling factor should be numeric.")
            
            expect_error(count(csds) <- "toto", 
                    "Counts should be numeric.")
            
        })



test_that("Slot is correctly modified if a valid score is given", {
            
            scalingFactor(csds) <- 1
            count(csds) <- 3
            
            expect_equal(getScalingFactor(csds), 1)
            expect_equal(getCount(csds), 3)
            
        })


test_that("The input bigwig file is modified correctly", {
            
            bigWigFile(csds) <- file_vec[4]
            
            expect_equal(getBigWigFile(csds), file_vec[4])
            
            expect_error(bigWigFile(csds) <- file_empty, 
                    paste0(file_empty, " should contain only one point for ",
                            "the extension."))
        })


test_that("The plotSetArrayList is modified correctly",{
            
            expect_error(averageBindingValues(csds) <- 3, 
                    "SetArrayList should be of type list.")
            
            expect_error(averageBindingValues(csds) <- list(), 
                    "The list of binding values is empty.")
            
            expect_error(averageBindingValues(csds) <- list(c(1,2)), 
                 "All objects in setArray list must be of type PlotSetArray.")
            
        })


test_that("The list of matrices is modified correctly", {
            
            expect_error(matBindingValues(csds) <- 3, 
                    "matBinding should be of type list.")
            
            expect_error(matBindingValues(csds) <- list(),
                    "The list of binding values is empty.")
            
            expect_error(matBindingValues(csds) <- list(c(1,2),
                   "All objects in matBindingValList must be of type matrix."))
        })


## For ChIPSeqSpikeDatasetBoost object

context("Scaling factors of ChIPSeqSpikeDatasetBoost objects are changed 
        correctly.")

test_that("Negative value cannot be given as scaling factor or count", {
            
            expect_error(scalingFactor(csdsboost) <- -2, 
                    "object@inputScalingFactor should be a positive number or 0.")
            
            expect_error(count(csdsboost) <- -2,
                    "object@inputCount should be a positive number or 0.")
            
        })


test_that("Non numeric value is not accepted as scaling factor or count",{
            
            expect_error(scalingFactor(csdsboost) <- "toto", 
                    "Scaling factor should be numeric.")
            
            expect_error(count(csdsboost) <- "toto", 
                    "Counts should be numeric.")
            
        })



test_that("Slot is correctly modified if a valid score is given", {
            
            scalingFactor(csdsboost) <- 1
            count(csdsboost) <- 3
            
            expect_equal(getScalingFactor(csdsboost), 1)
            expect_equal(getCount(csdsboost), 3)
            
        })


test_that("The input bigwig file is modified correctly", {
            
            bigWigFile(csdsboost) <- file_vec[4]
            
            expect_equal(getBigWigFile(csdsboost), file_vec[4])
            
        })


test_that("The plotSetArrayList is modified correctly",{
            
            expect_error(averageBindingValues(csdsboost) <- 3, 
                    "SetArrayList should be of type list.")
            
            expect_error(averageBindingValues(csdsboost) <- list(), 
                    "The list of binding values is empty.")
            
            expect_error(averageBindingValues(csdsboost) <- list(c(1,2)), 
                    "All objects in setArray list must be of type PlotSetArray.")
            
        })


test_that("The list of matrices is modified correctly", {
            
            expect_error(matBindingValues(csds) <- 3, 
                    "matBinding should be of type list.")
            
            expect_error(matBindingValues(csds) <- list(),
                    "The list of binding values is empty.")
            
            expect_error(matBindingValues(csds) <- list(c(1,2),
                            "All objects in matBindingValList must be of type matrix."))
        })


## For Experiment object

context("Scaling factors of Experiment objects are changed correctly.")

test_that("Negative value cannot be given as scaling factor or count", {
            
            expect_error(scalingFactor(e) <- -2, 
           "object@endogenousScalingFactor should be a positive number or 0")
            
            expect_error(exogenousScalingFactor(e) <- -2, 
           "object@exogenousScalingFactor should be a positive number or 0")
   
            expect_error(count(e) <- -2, 
           "object@endoCount should be a positive number or 0")
   
            expect_error(exoCount(e) <- -2, 
           "object@exoCount should be a positive number or 0")
   
        })

test_that("Non numeric value is not accepted as scaling factor or count",{
        
        expect_error(scalingFactor(e) <- "toto", 
                "Scaling factor should be numeric.")
        
        expect_error(exogenousScalingFactor(e) <- "toto", 
                "Scaling factor should be numeric.")
    
        expect_error(count(e) <- "toto", 
                "Counts should be numeric.")
    
        expect_error(exoCount(e) <- "toto", 
                "Counts should be numeric.")
    })


test_that("Slot is correctly modified if a valid score is given", {
            
            scalingFactor(e) <- 1
            exogenousScalingFactor(e) <- 2
            count(e) <- 3
            exoCount(e) <- 4
            
            expect_equal(getScalingFactor(e), 1)
            expect_equal(getExogenousScalingFactor(e), 2)
            expect_equal(getCount(e), 3)
            expect_equal(getExoCount(e), 4)
            
        })


test_that("The bigwig file is modified correctly", {
            
            bigWigFile(e) <- file_vec[3]
            
            expect_equal(getBigWigFile(e), file_vec[3])
            
            expect_error(bigWigFile(e) <- file_empty, 
                    paste0(file_empty, " should contain only one point for ",
                            "the extension."))
        })



## For ExperimentLoaded object

context("Scaling factors of Experiment objects are changed correctly.")

test_that("Negative value cannot be given as scaling factor or count", {
            
            expect_error(scalingFactor(eboost) <- -2, 
                    "object@endogenousScalingFactor should be a positive number or 0")
            
            expect_error(exogenousScalingFactor(eboost) <- -2, 
                    "object@exogenousScalingFactor should be a positive number or 0")
            
            expect_error(count(eboost) <- -2, 
                    "object@endoCount should be a positive number or 0")
            
            expect_error(exoCount(eboost) <- -2, 
                    "object@exoCount should be a positive number or 0")
            
        })

test_that("Non numeric value is not accepted as scaling factor or count",{
            
            expect_error(scalingFactor(eboost) <- "toto", 
                    "Scaling factor should be numeric.")
            
            expect_error(exogenousScalingFactor(eboost) <- "toto", 
                    "Scaling factor should be numeric.")
            
            expect_error(count(eboost) <- "toto", 
                    "Counts should be numeric.")
            
            expect_error(exoCount(eboost) <- "toto", 
                    "Counts should be numeric.")
        })


test_that("Slot is correctly modified if a valid score is given", {
            
            scalingFactor(eboost) <- 1
            exogenousScalingFactor(eboost) <- 2
            count(eboost) <- 3
            exoCount(eboost) <- 4
            
            expect_equal(getScalingFactor(eboost), 1)
            expect_equal(getExogenousScalingFactor(eboost), 2)
            expect_equal(getCount(eboost), 3)
            expect_equal(getExoCount(eboost), 4)
            
        })


test_that("The bigwig file is modified correctly", {
            
            bigWigFile(eboost) <- file_vec[3]
            
            expect_equal(getBigWigFile(eboost), file_vec[3])
            
        })
