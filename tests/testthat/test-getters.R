## Preparing data
file_vec <- system.file("extdata", 
        c("bam_files/H3K79me2_0_dm3-filtered.bam", 
                "bam_files/H3K79me2_0_hg19-filtered.bam", 
                "bigwig_files/H3K79me2_0-filtered.bw", 
                "bigwig_files/input_0-filtered.bw", 
                "bam_files/input_0_hg19-filtered.bam"),
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
        endoNb = 3, exoNb = 4);

csds <- ChIPSeqSpikeDataset(endogenousBam_vec = file_vec[2],
        exogenousBam_vec = file_vec[1],
        bigWigFile_endogenous_vec = file_vec[3], 
        inputBigWigFile = file_vec[4], 
        inputBamFile = file_vec[5], 
        expnames = "H3K79me2_0", inputSF = 0.7, inputNb = 1000);

csdsboost <- ChIPSeqSpikeDatasetBoost(endogenousBam_vec = file_vec[2],
        exogenousBam_vec = file_vec[1],
        bigWigFile_endogenous_vec = file_vec[3], 
        inputBigWigFile = file_vec[4], 
        inputBamFile = file_vec[5], 
        expnames = "H3K79me2_0", inputSF = 0.7, inputNb = 1000);

spikesumtest <- matrix(c(0,0.7,0,NA,0,1000,0,NA), ncol=4, nrow=2);
colnames(spikesumtest) <- c("endoScalFact", "exoScalFact", "endoCount", 
        "exoCount");
rownames(spikesumtest) <- c("H3K79me2_0", "input");

ratiotest <- matrix(0);
colnames(ratiotest) <- "Percentage Exo";
rownames(ratiotest) <- "H3K79me2_0";

## Testing accessors to all types of object

context("Testing proper accession to objects");

test_that("Accessors of Experiment object are valid", {
            
            expect_equal(getBam(e), file_vec[2]);
            expect_equal(getExogenousBam(e), file_vec[1]);
            expect_equal(getBigWigFile(e), file_vec[3]);
            expect_equal(getExpName(e), "H3K79me2_0");
            expect_equal(getScalingFactor(e), 1);
            expect_equal(getExogenousScalingFactor(e), 2);
            expect_equal(getCount(e), 3);
            expect_equal(getExoCount(e), 4);
        });

test_that("Accessors of ExperimentLoaded object are valid", {
            
            expect_equal(getBam(eboost), file_vec[2]);
            expect_equal(getExogenousBam(eboost), file_vec[1]);
            expect_equal(getBigWigFile(eboost), file_vec[3]);
            expect_equal(getExpName(eboost), "H3K79me2_0");
            expect_equal(getScalingFactor(eboost), 1);
            expect_equal(getExogenousScalingFactor(eboost), 2);
            expect_equal(getCount(eboost), 3);
            expect_equal(getExoCount(eboost), 4);
        });

test_that("Accessors to ChIPSeqSpikeDataset are valid", {
            
            expect_equal(getBam(csds), file_vec[5]);
            expect_equal(getBigWigFile(csds), file_vec[4]);
            expect_equal(as.character(getExperimentList(csds)), file_vec[3]);
            expect_equal(getExpName(csds), "H3K79me2_0");
            expect_equal(getScalingFactor(csds), 0.7);
            expect_equal(getCount(csds), 1000);
            expect_equal(spikeSummary(csds), spikesumtest);
            suppressWarnings(expect_equal(getRatio(csds), ratiotest))
        })

test_that("Accessors to ChIPSeqSpikeDatasetBoost are valid", {
            
            expect_equal(getBam(csdsboost), file_vec[5]);
            expect_equal(getBigWigFile(csdsboost), file_vec[4]);
            expect_equal(as.character(getExperimentList(csdsboost)), 
                    file_vec[3]);
            expect_equal(getExpName(csdsboost), "H3K79me2_0");
            expect_equal(getScalingFactor(csdsboost), 0.7);
            expect_equal(getCount(csdsboost), 1000);
            expect_equal(spikeSummary(csdsboost), spikesumtest);
            suppressWarnings(expect_equal(getRatio(csds), ratiotest));
        })
