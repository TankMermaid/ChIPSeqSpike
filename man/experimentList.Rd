\docType{methods}
\name{experimentList}
\alias{experimentList}
\alias{experimentList<-}
\alias{experimentList<-,ChIPSeqSpikeDataset-method}
\alias{experimentList<-,ChIPSeqSpikeDatasetBoost-method}


\title{Set the list of Experiment and ExperimentLoaded objects}

\description{
    Set the list of Experiment and ExperimentLoaded objects to the 
    ChIPSeqSpikeDataset or ChIPSeqSpikeDatasetBoost object
}

\usage{
experimentList(theObject) <- value
            
            
\S4method{experimentList}{ChIPSeqSpikeDataset}(theObject) <- value

\S4method{experimentList}{ChIPSeqSpikeDatasetBoost}(theObject) <- value
}

\arguments{
    \item{theObject}{A ChIPSeqSpikeDataset or ChIPSeqSpikeDatasetBoost 
    object}
 
    \item{value}{A list of Experiment and ExperimentLoaded objects}
}

\value{
  The modified object is returned
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

experimentList(csds[[1]]) <- getExperimentList(csds[[1]])

unlink("test_chipseqspike/", recursive = TRUE)

}

\seealso{
    \code{\link{getExperimentList}}
}

\author{Nicolas Descostes}
