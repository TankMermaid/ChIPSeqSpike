\name{ChIPSeqSpikeDatasetList-class}
\docType{class}
\alias{ChIPSeqSpikeDatasetList-class}
\alias{ChIPSeqSpikeDatasetList}
\alias{[[,ChIPSeqSpikeDatasetList,ANY,ANY-method}
\alias{[[<-,ChIPSeqSpikeDatasetList,ANY,ANY-method}
\title{ChIPSeqSpikeDatasetList Reference Class}

\description{
 Class containing a list of ChIPSeqSpikeDataset objects
}

\section{Fields}{
  \describe{
    datasetList: A list of ChIPSeqSpikeDataset objects
    }
}

\section{Constructor}{
  \describe{
    ChIPSeqSpikeDatasetList(dataset_list, verbose)
  }
}

\section{Arguments}{
  \describe{
    \item{dataset_list}{A properly formatted list of information needed 
                               to create the object (see details).}

    \item{verbose}{Indicate if processing messages should be output.}
  }
}

\section{Getters}{
  \describe{
    \item{\code{getBigWigFile}}{Returns a vector of all bigWig paths 
                                corresponding to all files specified in 
                                info.csv}

    \item{\code{x[[i]]}}{Get the ChIPSeqSpikeDataset object at index i}
  }
}

\section{Setters}{
  \describe{
    \item{\code{x[[i]] <- value}}{Set \code{value} to ChIPSeqSpikeDatasetList i}
  }
}

\section{Details}{
This class enables to process datasets containing different input files. It will
creates a list of ChIPSeqSpikeDataset objects, each of them containing a 
different input.

The above indicated constructor should not be used directly. One would rather
create the object by calling the "meta-constructor" spikeDataset function which
takes a info.csv file as input. spikeDataset function formats properly the 
different information into a dataset_list which is submitted to the 
ChIPSeqSpikeDatasetList constructor.
}


\examples{
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike");
bam_path <- system.file(c("extdata/bam_files"), package="ChIPSeqSpike");
bigwig_path <- system.file(c("extdata/bigwig_files"), package="ChIPSeqSpike");

csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path);

is(csds);
csds
}

\seealso{
  \code{\link{ChIPSeqSpikeDataset-class}}
  \code{\link{spikeDataset}}
  \code{\link{spikeSummary}}
  \code{\link{getRatio}}
}

\author{Nicolas Descostes}
