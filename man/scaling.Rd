\docType{methods}
\name{scaling}
\alias{scaling}
\alias{scaling,ChIPSeqSpikeDataset-method}
\alias{scaling,ChIPSeqSpikeDatasetBoost-method}
\alias{scaling,ChIPSeqSpikeDatasetList-method}
\alias{scaling,ChIPSeqSpikeDatasetListBoost-method}

\title{Applies different type of scaling/normalization procedures}

\description{
    Applies or reverse different type of scaling/normalization to bigwig files 
contained in 'ChIPSeqSpikeDataset', 'ChIPSeqSpikeDatasetBoost', 
'ChIPSeqSpikeDatasetList' or 'ChIPSeqSpikeDatasetListBoost' objects
}

\usage{

scaling(theObject, reverse = FALSE, type = "endo", verbose = TRUE, 
        outputFolder = NULL)
            
\S4method{scaling}{ChIPSeqSpikeDataset}(theObject, reverse = FALSE, 
type = "endo", verbose = TRUE, outputFolder = NULL)

\S4method{scaling}{ChIPSeqSpikeDatasetBoost}(theObject, reverse = FALSE, 
type = "endo", verbose = TRUE, outputFolder = NULL)

\S4method{scaling}{ChIPSeqSpikeDatasetList}(theObject, reverse = FALSE, 
type = "endo", verbose = TRUE, outputFolder = NULL)

\S4method{scaling}{ChIPSeqSpikeDatasetListBoost}(theObject, reverse = FALSE, 
type = "endo", verbose = TRUE, outputFolder = NULL)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}

    \item{reverse}{If TRUE, reverse the scaling applied previously (see 
                   details). Default is FALSE}

    \item{type}{Should be 'endo' or 'exo' to apply the endogenous or exogenous 
                scaling factor respectively.}
 
    \item{verbose}{If FALSE, do not output processing messages. Default is 
                   TRUE}

    \item{outputFolder}{Define the folder where scaled bigwig are output. 
                        Default is NULL (see details).}
}

\value{
    Return an object of the same class of the input object and output scaled 
    bigwig files if not in boost mode (see ?ChIPSeqSpikeDatasetBoost).
    
    '-RPM', 'reverse' or 'spiked' suffixes are added to the bigwig file name if
    performing RPM scaling, reversing RPM scaling or applying exogenous scaling
    factors respectively.
}

\details{
  Scaling is performed after estimating scaling factors with 
  '?estimateScalingFactors'.
  
  Different type of scaling can be performed according to the steps of spike-in
  normalization. The first type is to transform the data by performing a RPM 
  scaling. The second type is to reverse the RPM scaling after having performed
  input subtraction with '?inputSubtraction'. The third and last type is to 
  apply the exogenous scaling factor which finalizes the spike-in normalization
  (see example below for parameter settings).
  
  To perform all steps in one call and in the right order, one can consider 
  using the '?spikePipe' function.
  
  If not in boost mode, the bigwig files are written to the folder containing 
  the currently processed bigwig files. In boost mode, use the method 
  ?exportBigWigs to output the transformed files.
  
  If outputFolder is not NULL, the original bigwig files should be copied to 
  this folder before performing the analysis. This parameter was created to 
  test the package with the provided files in extdata/.
}

\section{Methods (by class)}{
\itemize{

\item \code{ChIPSeqSpikeDataset}: Method for signature theObject=
                                                          'ChIPSeqSpikeDataset'

\item \code{ChIPSeqSpikeDatasetBoost}: Method for signature theObject=
                                                     'ChIPSeqSpikeDatasetBoost'

\item \code{ChIPSeqSpikeDatasetList}: Method for signature theObject=
                                                      'ChIPSeqSpikeDatasetList'

\item \code{ChIPSeqSpikeDatasetListBoost}: Method for signature theObject=
                                                 'ChIPSeqSpikeDatasetListBoost'
}}

\examples{

## Mock example on a restricted number of reads
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike");
bam_path <- system.file(c("extdata/bam_files"), package="ChIPSeqSpike");
bigwig_path <- system.file(c("extdata/bigwig_files"), package="ChIPSeqSpike");

csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path);

## Copying test files to the current folder
originalBW_vec <- as.character(getBigWigFile(csds));
dir.create("./test_chipseqspike");
result <- file.copy(originalBW_vec, "test_chipseqspike");

csds <- estimateScalingFactors(csds);

## Apply RPM scaling
csds <- scaling(csds, outputFolder = "test_chipseqspike");

## Apply input subtraction
csds <- inputSubtraction(csds);

## Reverse RPM scaling after input subtraction
csds <- scaling(csds, reverse = TRUE);

## Apply exogenous scaling factors
csds <- scaling(csds, type = "exo");

## Delete all files generated in this example
unlink("test_chipseqspike/", recursive = TRUE)
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link{ChIPSeqSpikeDatasetBoost}}
    \code{\link{spikePipe}}
    \code{\link{inputSubtraction}}
    \code{\link{exportBigWigs}}
}

\author{Nicolas Descostes}
