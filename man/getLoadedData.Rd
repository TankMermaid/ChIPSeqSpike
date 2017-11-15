\docType{methods}
\name{getLoadedData}
\alias{getLoadedData}
\alias{getLoadedData,ChIPSeqSpikeDatasetBoost-method}
\alias{getLoadedData,ExperimentLoaded-method}


\title{Get the endogenous reference genome binding scores of an experiment}

\description{
    Accesses and returns the binding scores of an experiment or input DNA 
    experiment on the endogenous reference genome. (only available in boost 
    mode).
}

\usage{
getLoadedData(theObject)
            
            
\S4method{getLoadedData}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getLoadedData}{ExperimentLoaded}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDatasetBoost} or 
                    \code{ExperimentLoaded} object}
}

\value{
  A GRanges object of binding scores
}
\details{
  If the object is ChIPSeqSpikeDatasetBoost, getLoadedData returns a GRanges 
  object of binding scores of the input DNA experiment.
  
  If the object is ExperimentLoaded, getLoadedData returns a GRanges 
  object of binding scores of the experiment.
}

\examples{
file_vec <- system.file("extdata", 
                        c("bam_files/H3K79me2_0_dm3-filtered.bam", 
                        "bam_files/H3K79me2_0_hg19-filtered.bam", 
                        "bigwig_files/H3K79me2_0-filtered.bw", 
                        "bigwig_files/input_0-filtered.bw", 
                        "bam_files/input_0_hg19-filtered.bam"),
                        package="ChIPSeqSpike")

csds <- ChIPSeqSpikeDatasetBoost(endogenousBam_vec = file_vec[2],
                         exogenousBam_vec = file_vec[1],
                         bigWigFile_endogenous_vec = file_vec[3], 
                         inputBigWigFile = file_vec[4], 
                         inputBamFile = file_vec[5], 
                         expnames = "H3K79me2_0");

getLoadedData(csds);
}

\author{Nicolas Descostes}
