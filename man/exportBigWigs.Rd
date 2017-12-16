\docType{methods}
\name{exportBigWigs}
\alias{exportBigWigs}
\alias{exportBigWigs,ChIPSeqSpikeDatasetBoost-method}
\alias{exportBigWigs,ChIPSeqSpikeDatasetListBoost-method}

\title{Export bigwig files from values contained in a boost mode object}

\description{
    In boost mode, binding values of experiments are stored in the form of 
    GRanges tables in the object. exportBigWigs output these values giving a 
    proper bigwig file name corresponding to the transformations already 
    performed.
}

\usage{
exportBigWigs(theObject, verbose = TRUE)
            
            
\S4method{exportBigWigs}{ChIPSeqSpikeDatasetBoost}(theObject, verbose = TRUE)

\S4method{exportBigWigs}{ChIPSeqSpikeDatasetListBoost}(theObject, 
verbose = TRUE)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset in boost mode (see 
                     ?spikeDataset)}
 
    \item{verbose}{If FALSE, do not output processing messages. Default is 
                   TRUE}
}

\value{
    Output bigwig files of binding values.
    
    The suffix of the bigwig file reflects the transformation steps performed 
    on the object. If all steps were performed, the file name will be of the 
    form: 'expName-RPM-BGSub-reverse-spiked.bw'. The suffixes 'RPM', 'BGSub', 
    'reverse' and 'spiked' stands for RPM scaling, input subtraction, RPM 
    scaling reversal and exogenous scaling respectively.
}

\section{Methods (by class)}{
\itemize{

\item \code{ChIPSeqSpikeDatasetBoost}: Method for signature theObject=
                                                     'ChIPSeqSpikeDatasetBoost'

\item \code{ChIPSeqSpikeDatasetListBoost}: Method for signature theObject=
                                                 'ChIPSeqSpikeDatasetListBoost'
}}

\examples{
## Mock example on a restricted number of reads
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike")
bam_path <- system.file(c("extdata/bam_files"), package="ChIPSeqSpike")
bigwig_path <- system.file(c("extdata/bigwig_files"), package="ChIPSeqSpike")

csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path, boost = TRUE)

## Creating test folder
dir.create("./test_chipseqspike")

csds <- estimateScalingFactors(csds)

## Apply RPM scaling
csds <- scaling(csds, outputFolder = "test_chipseqspike")

## output RPM scaled files
exportBigWigs(csds)

## Apply input subtraction
csds <- inputSubtraction(csds)

## output input subtracted files
exportBigWigs(csds)

## Delete all files generated in this example
unlink("test_chipseqspike/", recursive = TRUE)
}


\author{Nicolas Descostes}
