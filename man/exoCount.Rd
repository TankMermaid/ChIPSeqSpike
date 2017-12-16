\docType{methods}
\name{exoCount}
\alias{exoCount}
\alias{exoCount<-}
\alias{exoCount<-,Experiment-method}
\alias{exoCount<-,ExperimentLoaded-method}


\title{Set the number of reads associated to an experiment}

\description{
    Set the number of exogenous reads associated to an experiment. This is used
    to compute scaling factors.
}

\usage{
exoCount(theObject) <- value
            
            
\S4method{exoCount}{Experiment}(theObject) <- value

\S4method{exoCount}{ExperimentLoaded}(theObject) <- value
}

\arguments{
    \item{theObject}{An Experiment or ExperimentLoaded object}
 
    \item{value}{A numeric representing the number of mapped reads to the 
                exogenous reference genome}
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

getExoCount(csds[[1]][[1]])
exoCount(csds[[1]][[1]]) <- 5
getExoCount(csds[[1]][[1]])

unlink("test_chipseqspike/", recursive = TRUE)
}

\seealso{
    \code{\link{count}}
    \code{\link{estimateScalingFactors}}
}

\author{Nicolas Descostes}
