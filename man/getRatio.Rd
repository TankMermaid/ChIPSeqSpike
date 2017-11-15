\docType{methods}
\name{getRatio}
\alias{getRatio}
\alias{getRatio,Experiment-method}
\alias{getRatio,ExperimentLoaded-method}
\alias{getRatio,ChIPSeqSpikeDataset-method}
\alias{getRatio,ChIPSeqSpikeDatasetBoost-method}
\alias{getRatio,ChIPSeqSpikeDatasetList-method}
\alias{getRatio,ChIPSeqSpikeDatasetListBoost-method}

\title{Output the percentage of exogenous DNA compared to that of endogenous 
DNA}

\description{
    Output the percentage of exogenous DNA compared to that of endogenous DNA 
}

\usage{

    getRatio(theObject)

\S4method{getRatio}{Experiment}(theObject)
\S4method{getRatio}{ExperimentLoaded}(theObject)
\S4method{getRatio}{ChIPSeqSpikeDataset}(theObject)
\S4method{getRatio}{ChIPSeqSpikeDatasetList}(theObject)
\S4method{getRatio}{ChIPSeqSpikeDatasetBoost}(theObject)
\S4method{getRatio}{ChIPSeqSpikeDatasetListBoost}(theObject)
}

\arguments{
    \item{theObject}{ChIPSeqSpike dataset (see ?spikeDataset)}
}

\value{
    A numeric matrix
}

\details{
  The rows represent, for each experiment of the dataset, the percentage of 
exogenous DNA defined as the number of aligned exogenous reads compared to the 
total number of reads.

  The method 'getRatio' will throw a warning if the percentage of exogenous DNA
 represents less than 2 percent or more than 25 percent of the endogenous DNA. 
Less than 2 percent of exogenous DNA does not guarantee a proper scaling. Large
 amount of exogenous DNA should not impact the scaling procedure but is worth 
notifying to the user[1].
}

\section{Methods (by class)}{
\itemize{
\item \code{Experiment}: Method for signature theObject = 'Experiment'

\item \code{ExperimentLoaded}: Method for signature theObject = 
                                                             'ExperimentLoaded'

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
## Mock example on files samples
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike");
bam_path <- system.file(c("extdata/bam_files"), package="ChIPSeqSpike");
bigwig_path <- system.file(c("extdata/bigwig_files"), package="ChIPSeqSpike");
gff_vec <- system.file("extdata/test_coord.gff", package="ChIPSeqSpike");
genome_name <- "hg19";

csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path);
csds <- estimateScalingFactors(csds);
getRatio(csds);

## Results on the complete files
data("ratio");
print(ratio);
}

\seealso{
    \code{\link{spikeSummary}}
}

\references{
  [1] Orlando et al, "Quantitative ChIP-Seq normalization reveals global 
  modulation of the epigenome", Cell Rep, 2014.
  
}

\author{Nicolas Descostes}
