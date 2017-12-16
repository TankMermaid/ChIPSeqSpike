\docType{methods}
\name{matBindingValues}
\alias{matBindingValues}
\alias{matBindingValues<-}
\alias{matBindingValues<-,ChIPSeqSpikeDataset-method}
\alias{matBindingValues<-,ChIPSeqSpikeDatasetBoost-method}
\alias{matBindingValues<-,ChIPSeqSpikeCore-method}


\title{Set the matrices of binding values}

\description{
    Set a list of matrices of binding values that are used to plot boxplots and
    correlation plots.
}

\usage{
matBindingValues(theObject) <- value
            
            
\S4method{matBindingValues}{ChIPSeqSpikeDataset}(theObject) <- value

\S4method{matBindingValues}{ChIPSeqSpikeDatasetBoost}(theObject) <- value

\S4method{matBindingValues}{ChIPSeqSpikeCore}(theObject) <- value
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset} or 
                     \code{ChIPSeqSpikeDatasetBoost} object}
 
    \item{value}{A list of matrices}
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
              verbose = FALSE, outputFolder = output_folder)

new_list <- list(matrix(seq_len(10)), matrix(2:18))
matBindingValues(csds[[1]]) <- new_list

unlink("test_chipseqspike/", recursive = TRUE)
}

\seealso{
    \code{\link{plotCor}}
    \code{\link{boxplotSpike}}
}

\author{Nicolas Descostes}
