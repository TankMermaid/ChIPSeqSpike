\docType{methods}
\name{exogenousScalingFactor}
\alias{exogenousScalingFactor}
\alias{exogenousScalingFactor<-}
\alias{exogenousScalingFactor<-,Experiment-method}
\alias{exogenousScalingFactor<-,ExperimentLoaded-method}


\title{Set the exogenous scaling factor associated to an experiment}

\description{
    Set the exogenous scaling factor associated to an experiment
}

\usage{
exogenousScalingFactor(theObject) <- value
            
            
\S4method{exogenousScalingFactor}{Experiment}(theObject) <- value

\S4method{exogenousScalingFactor}{ExperimentLoaded}(theObject) <- value
}

\arguments{
    \item{theObject}{An \code{Experiment} or \code{ExperimentLoaded} object}
 
    \item{value}{A numeric representing the exogenous scaling factor}
}

\details{
  A scaling factor is defined as:
  
  \itemize{
             \item 1/(bam_count/1000000)
  }
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

getExogenousScalingFactor(csds[[1]][[1]])
exogenousScalingFactor(csds[[1]][[1]]) <- 5
getExogenousScalingFactor(csds[[1]][[1]])

unlink("test_chipseqspike/", recursive = TRUE)
}

\seealso{
    \code{\link{scalingFactor}}
    \code{\link{estimateScalingFactors}}
}

\author{Nicolas Descostes}
