\docType{methods}
\name{estimateScalingFactors}
\alias{estimateScalingFactors}
\alias{estimateScalingFactors,Experiment-method}
\alias{estimateScalingFactors,ExperimentLoaded-method}
\alias{estimateScalingFactors,ChIPSeqSpikeDataset-method}
\alias{estimateScalingFactors,ChIPSeqSpikeDatasetBoost-method}
\alias{estimateScalingFactors,ChIPSeqSpikeDatasetList-method}
\alias{estimateScalingFactors,ChIPSeqSpikeDatasetListBoost-method}

\title{Compute scaling factors to perform spike-in normalization}

\description{
    Compute scaling factors for endogenous and exogenous experiment from 
'Experiment', 'ExperimentLoaded', 'ChIPSeqSpikeDataset', 
'ChIPSeqSpikeDatasetBoost', 'ChIPSeqSpikeDatasetList', and 
'ChIPSeqSpikeDatasetListBoost'
}

\usage{

estimateScalingFactors(theObject, paired = FALSE, verbose = TRUE)

\S4method{estimateScalingFactors}{Experiment}(theObject, paired = FALSE, 
verbose = TRUE)

\S4method{estimateScalingFactors}{ExperimentLoaded}(theObject, paired = FALSE, 
verbose = TRUE)

\S4method{estimateScalingFactors}{ChIPSeqSpikeDataset}(theObject, 
paired = FALSE, verbose = TRUE)

\S4method{estimateScalingFactors}{ChIPSeqSpikeDatasetBoost}(theObject, 
paired = FALSE, verbose = TRUE)

\S4method{estimateScalingFactors}{ChIPSeqSpikeDatasetList}(theObject, 
paired = FALSE, verbose = TRUE)

\S4method{estimateScalingFactors}{ChIPSeqSpikeDatasetListBoost}(theObject, 
paired = FALSE, verbose = TRUE)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}
    
    \item{paired}{Indicate if sequences are single- or paired-ended. Default is
                  FALSE}
    
    \item{verbose}{If False, do not output processing messages. Default is 
                   TRUE}
}

\value{
    Return an object of the same class of the input object containing computed 
    scaling factors.
}

\details{
  Estimating scaling factors is the first step to perform on a dataset. 
  A scaling factor is defined as:
  \itemize{
             \item 1/(bam_count/1000000)
  }
  bam_count being the number of reads aligned to the genome. The count is 
  determined for the endogenous and exogenous experiments.
  
  Scaling factors will the be applied to the bigwig files in the following 
  steps of the procedure. After estimating scaling factors, RPM normalization 
  should be performed. 
  
}


\section{Methods (by class)}{
\itemize{
\item \code{Experiment}: Method for signature theObject='Experiment'

\item \code{ExperimentLoaded}: Method for signature theObject='ExperimentLoaded
'

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

csds <- estimateScalingFactors(csds);
}

\seealso{
    \code{\link{spikeSummary}}
    \code{\link{scaling}}
    \code{\link{spikePipe}}
}

\author{Nicolas Descostes}
