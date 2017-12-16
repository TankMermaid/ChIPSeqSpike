\docType{methods}
\name{scalingFactor}
\alias{scalingFactor}
\alias{scalingFactor<-}
\alias{scalingFactor<-,ChIPSeqSpikeDataset-method}
\alias{scalingFactor<-,ChIPSeqSpikeDatasetBoost-method}
\alias{scalingFactor<-,Experiment-method}
\alias{scalingFactor<-,ExperimentLoaded-method}
\alias{scalingFactor<-,ChIPSeqSpikeCore-method}


\title{Set the endogenous scaling factor associated to an experiment}

\description{
    Set the endogenous scaling factor associated to an experiment or its 
    corresponding input DNA experiment.
}

\usage{
scalingFactor(theObject) <- value
            
            
\S4method{scalingFactor}{ChIPSeqSpikeDataset}(theObject) <- value

\S4method{scalingFactor}{ChIPSeqSpikeDatasetBoost}(theObject) <- value

\S4method{scalingFactor}{Experiment}(theObject) <- value

\S4method{scalingFactor}{ExperimentLoaded}(theObject) <- value

\S4method{scalingFactor}{ChIPSeqSpikeCore}(theObject) <- value
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset}, 
                     \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment} or 
                     \code{ExperimentLoaded} object}
 
    \item{value}{A numeric representing the endogenous scaling factor}
}

\value{
  The modified object is returned
}

\details{
  A scaling factor is defined as:
  
  \itemize{
             \item 1/(bam_count/1000000)
  }
  
  If the object is ChIPSeqSpikeCore, ChIPSeqSpikeDataset or 
  ChIPSeqSpikeDatasetBoost, scalingFactor will set the endogenous scaling 
  factor of the input DNA experiment of the dataset.
  
  If the object is Experiment or ExperimentLoaded, scalingFactor will set the 
  endogenous scaling factor of the experiment.
}

\examples{
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike")
bam_path <- system.file("extdata/bam_files", package="ChIPSeqSpike")
bigwig_path <- system.file("extdata/bigwig_files", package="ChIPSeqSpike")
gff_vec <- system.file("extdata/test_coord.gff", package="ChIPSeqSpike")
genome_name <- "hg19"
output_folder <- "test_chipseqspike"
bigwig_files <- system.file("extdata/bigwig_files", 
                            c("H3K79me2_0-filtered.bw",
                              "H3K79me2_100-filtered.bw",
                              "H3K79me2_50-filtered.bw",
                              "input_0-filtered.bw",
                              "input_100-filtered.bw",
                              "input_50-filtered.bw"), package="ChIPSeqSpike")

## Copying example files
dir.create("./test_chipseqspike")
result <- file.copy(bigwig_files, "test_chipseqspike")

csds <- spikePipe(info_file_csv, bam_path, bigwig_path, gff_vec, genome_name, 
              verbose = TRUE, outputFolder = output_folder)

getScalingFactor(csds[[1]][[1]])
scalingFactor(csds[[1]][[1]]) <- 15
getScalingFactor(csds[[1]][[1]])

unlink("test_chipseqspike/", recursive = TRUE)
}

\seealso{
    \code{\link{exogenousScalingFactor}}
    \code{\link{estimateScalingFactors}}
}

\author{Nicolas Descostes}
