\docType{methods}
\name{inputSubtraction}
\alias{inputSubtraction}
\alias{inputSubtraction,ChIPSeqSpikeDataset-method}
\alias{inputSubtraction,ChIPSeqSpikeDatasetBoost-method}
\alias{inputSubtraction,ChIPSeqSpikeDatasetList-method}
\alias{inputSubtraction,ChIPSeqSpikeDatasetListBoost-method}

\title{Substracts binding scores of input DNA to experiment binding scores}

\description{
    Substracts binding scores of input DNA to experiment binding scores. This 
    step enables to remove artifactual signal.
}

\usage{
inputSubtraction(theObject, verbose = TRUE)
            
            
\S4method{inputSubtraction}{ChIPSeqSpikeDataset}(theObject, verbose = TRUE)

\S4method{inputSubtraction}{ChIPSeqSpikeDatasetBoost}(theObject, 
verbose = TRUE)

\S4method{inputSubtraction}{ChIPSeqSpikeDatasetList}(theObject, verbose = TRUE)

\S4method{inputSubtraction}{ChIPSeqSpikeDatasetListBoost}(theObject, 
verbose = TRUE)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}
 
    \item{verbose}{If FALSE, do not output processing messages. Default is 
                   TRUE}
}

\value{
    Return an object of the same class of the input object with subtracted 
    experiment scores.
    
    A 'BGSub' suffix is added to the bigwig file name.
}

\details{
  When immunoprecipitating (IP) DNA bound by a given protein, a control is 
  needed to distinguish background noise from true signal. This is typically 
  achieved by performing a mock IP omitting the use of antibody. After 
  sequencing, one can notice peaks of signal above background. These peaks have
  to be removed from analysis since they represent false positives.
  
  The inputSubtraction function reads bigwig files into GRanges objects that 
  are used to perform the subtraction. In boost mode (ChIPSeqSpikeDatasetBoost 
  and ChIPSeqSpikeDatasetListBoost), The reading/writing steps are omitted.
  
  If not in boost mode, the input DNA subtracted bigwig files are written to 
  the folder containing the currently processed bigwig files. In boost mode, 
  use the method exportBigWigs to output the transformed files. 
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

## Delete all files generated in this example
unlink("test_chipseqspike/", recursive = TRUE);
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link{exportBigWigs}}
    \code{\link{spikePipe}}
    \code{\link{scaling}}
}

\author{Nicolas Descostes}
