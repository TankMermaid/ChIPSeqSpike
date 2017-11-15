\name{Experiment-class}
\docType{class}
\alias{Experiment-class}
\alias{Experiment}
\alias{[[<-,Experiment,ANY,ANY-method}
\title{Experiment Reference Class}

\description{
    Object containing information about experiment. It constitutes an element 
    of a list held by ChIPSeqSpikeDataset class objects
}


\section{Fields}{
  \describe{

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

   Experiment(endogenousBamFilePath, exogenousBamFilePath, 
        bigWigFilePath, name, endoScalingFactor = 0, exoScalingFactor = 0, 
        endoNb = 0, exoNb = 0)
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

    \item{endoNb}{Number of reads aligned to the exogenous genome.
                         Default is 0.}
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
  }
}

\section{Setters}{
  \describe{
    \item{\code{scalingFactor}}{Modifies the endogenous scaling factor value}

    \item{\code{exogenousScalingFactor}}{Modifies the exogenous scaling factor 
                                         value}

    \item{\code{count}}{Modifies the endogenous count value}

    \item{\code{exoCount}}{Modifies the exogenous count value}

    \item{\code{bigWigFile}}{Modifies the endogenous bigWig file path}

    \item{\code{x[[i]] <- value}}{Set \code{value} to experiment i}
  }
}

\examples{

file_vec <- system.file("extdata", 
                        c("bam_files/H3K79me2_0_dm3-filtered.bam", 
                        "bam_files/H3K79me2_0_hg19-filtered.bam", 
                        "bigwig_files/H3K79me2_0-filtered.bw"),
                        package="ChIPSeqSpike")

exp <- Experiment(endogenousBamFilePath = file_vec[2], 
                  exogenousBamFilePath = file_vec[1], 
                  bigWigFilePath = file_vec[3], 
                  name = "H3K79me2_0")
}

\seealso{
  \code{\link{ExperimentLoaded-class}}
  \code{\link{ChIPSeqSpikeDataset-class}}
  \code{\link{spikeSummary}}
  \code{\link{getRatio}}
}


\author{Nicolas Descostes}
