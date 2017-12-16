\docType{methods}
\name{count}
\alias{count}
\alias{count<-}
\alias{count<-,ChIPSeqSpikeDataset-method}
\alias{count<-,ChIPSeqSpikeDatasetBoost-method}
\alias{count<-,Experiment-method}
\alias{count<-,ExperimentLoaded-method}
\alias{count<-,ChIPSeqSpikeCore-method}


\title{Set the number of reads associated to an experiment}

\description{
    Set the number of endogenous reads associated to an experiment. This is 
    used to compute scaling factors.
}

\usage{
count(theObject) <- value
            
            
\S4method{count}{ChIPSeqSpikeDataset}(theObject) <- value

\S4method{count}{ChIPSeqSpikeDatasetBoost}(theObject) <- value

\S4method{count}{Experiment}(theObject) <- value

\S4method{count}{ExperimentLoaded}(theObject) <- value

\S4method{count}{ChIPSeqSpikeCore}(theObject) <- value
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDataset}, 
                    \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment}, 
                    \code{ExperimentLoaded} or \code{ChIPSeqSpikeCore} object}
 
    \item{value}{A numeric representing the number of mapped reads}
}

\value{
  The modified object is returned
}

\details{
  If the object is ChIPSeqSpikeCore, ChIPSeqSpikeDataset or 
  ChIPSeqSpikeDatasetBoost, count will set the number of input DNA mapped reads
  to the endogenous reference genome.
  
  If the object is Experiment or ExperimentLoaded, count will set the number 
  of experiment mapped reads to the endogenous reference genome.
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

getCount(csds[[1]][[1]])
count(csds[[1]][[1]]) <- 10
getCount(csds[[1]][[1]])

unlink("test_chipseqspike/", recursive = TRUE)
}

\seealso{
    \code{\link{exoCount}}
    \code{\link{estimateScalingFactors}}
}

\author{Nicolas Descostes}
