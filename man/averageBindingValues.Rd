\docType{methods}
\name{averageBindingValues}
\alias{averageBindingValues}
\alias{averageBindingValues<-}
\alias{averageBindingValues<-,ChIPSeqSpikeDataset-method}
\alias{averageBindingValues<-,ChIPSeqSpikeDatasetBoost-method}
\alias{averageBindingValues<-,ChIPSeqSpikeCore-method}


\title{Set the binding values used to perform profiles and heatmaps}

\description{
    Set a plotSetArray object to the plotSetArrayList slot of the given object.
}

\usage{
averageBindingValues(theObject) <- value
            
            
\S4method{averageBindingValues}{ChIPSeqSpikeDataset}(theObject) <- value

\S4method{averageBindingValues}{ChIPSeqSpikeDatasetBoost}(theObject) <- value

\S4method{averageBindingValues}{ChIPSeqSpikeCore}(theObject) <- value
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDataset}, 
                     \code{ChIPSeqSpikeDatasetBoost} or \code{ChIPSeqSpikeCore}
                      object}
 
    \item{value}{A \code{PlotSetArray} object}
}

\value{
  The modified object is returned
}

\examples{
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike");
bam_path <- system.file("extdata/bam_files", package="ChIPSeqSpike");
bigwig_path <- system.file("extdata/bigwig_files", package="ChIPSeqSpike");
gff_vec <- system.file("extdata/test_coord.gff", package="ChIPSeqSpike");
genome_name <- "hg19";
output_folder <- "test_chipseqspike";
bigwig_files <- system.file("extdata/bigwig_files", 
                            c("H3K79me2_0-filtered.bw",
                              "H3K79me2_100-filtered.bw",
                              "H3K79me2_50-filtered.bw",
                              "input_0-filtered.bw",
                              "input_100-filtered.bw",
                              "input_50-filtered.bw"), package="ChIPSeqSpike");

## Copying example files
dir.create("./test_chipseqspike");
result <- file.copy(bigwig_files, "test_chipseqspike");

csds <- spikePipe(info_file_csv, bam_path, bigwig_path, gff_vec, genome_name, 
              verbose = TRUE, outputFolder = output_folder);

test <- getAverageBindingValues(csds[[1]]);
averageBindingValues(csds[[1]]) <- test;

unlink("test_chipseqspike/", recursive = TRUE);
}

\seealso{
    \code{\link{plotProfile}}
    \code{\link{plotTransform}}
    \code{\link{plotHeatmaps}}
    \code{\link[seqplots]{PlotSetArray}}
}

\author{Nicolas Descostes}
