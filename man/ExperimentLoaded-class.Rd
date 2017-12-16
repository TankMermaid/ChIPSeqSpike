\name{ExperimentLoaded-class}
\docType{class}
\alias{ExperimentLoaded-class}
\alias{ExperimentLoaded}
\title{ExperimentLoaded Reference Class}

\description{
    Class inheriting from Experiment and containing loaded binding values
}


\section{Fields}{
  \describe{

    loadedBigWigFile: GRanges object representing binding values

    endogenousBam: File path to the experiment BAM file aligned to the 
                  reference genome

    exogenousBam: File path to the experiment BAM file aligned to the 
                  exogenous genome

    bigWigFile: File path to the experiment BIGWIG file aligned to the 
                reference genome

    expName: Experiment name

    endogenousScalingFactor: Experiment scaling factor

    exogenousScalingFactor: Scaling factor obtained from exogenous DNA

    endoCount: Number of reads contained in the endogenous BAM file

    exoCount: Number of reads contained in the exogenous BAM file
  }
}


\section{Constructor}{
  \describe{

    ExperimentLoaded(endogenousBamFilePath, exogenousBamFilePath, 
        bigWigFilePath, name, endoScalingFactor = 0, exoScalingFactor = 0, 
        endoNb = 0, exoNb = 0, verbose = TRUE)
  }
}


\section{Arguments}{
  \describe{
    \item{endogenousBamFilePath}{Character vector of path to experiment 
                                        BAM file aligned to the reference 
                                        genome.}

    \item{exogenousBamFilePath}{Character vector of path to experiment 
                                       BAM file aligned to the exogenous 
                                       genome.}

    \item{bigWigFilePath}{Character vector of path to experiment bigWig 
                                 file aligned to the reference genome.}

    \item{name}{Character vector of the experiment name.}

    \item{endoScalingFactor}{Numeric scaling factor computed from the 
                                    number of reads aligned to the reference 
                                    genome. Default is 0.}

    \item{exoScalingFactor}{Numeric scaling factor computed from the 
                                   number of reads aligned to the exogenous 
                                    genome. Default is 0.}

    \item{endoNb}{Number of reads aligned to the reference genome.
                         Default is 0.}

    \item{exoNb}{Number of reads aligned to the exogenous genome.
                         Default is 0.}

    \item{verbose}{Indicate if processing messages should be output.
                         Default is TRUE.}
  }
}

\section{Getters}{
  \describe{
  
    \item{\code{getBam}}{Returns the endogenous BAM path}
    
    \item{\code{getExogenousBam}}{Returns the exogenous BAM path}

    \item{\code{getBigWigFile}}{Returns the endogenous bigWig path}

    \item{\code{getExpName}}{Returns a character vector of the experiment name}

    \item{\code{getScalingFactor}}{Returns the endogenous scaling factor}

    \item{\code{getExogenousScalingFactor}}{Returns the exogenous scaling 
                                           factor}

    \item{\code{getCount}}{Returns the number of reads aligned to the reference
                           genome}

    \item{\code{getExoCount}}{Returns the number of reads aligned to the 
                              exogenous genome}

    \item{\code{getLoadedData}}{Returns the GRanges object of binding values}
  }
}

\section{Setters}{
  \describe{
    \item{\code{bigWigFile}}{Modifies the endogenous bigWig file path}

    \item{\code{loadedData}}{Modifies the GRanges object of binding values}

  }
}

\examples{

file_vec <- system.file("extdata", 
                        c("bam_files/H3K79me2_0_dm3-filtered.bam", 
                        "bam_files/H3K79me2_0_hg19-filtered.bam", 
                        "bigwig_files/H3K79me2_0-filtered.bw"),
                        package="ChIPSeqSpike")

exp <- ExperimentLoaded(endogenousBamFilePath = file_vec[2], 
                        exogenousBamFilePath = file_vec[1], 
                        bigWigFilePath = file_vec[3], 
                        name = "H3K79me2_0")
}

\seealso{
  \code{\link{Experiment-class}}
  \code{\link{ChIPSeqSpikeDatasetBoost-class}}
  \code{\link{spikeSummary}}
  \code{\link{getRatio}}
}

\author{Nicolas Descostes}
