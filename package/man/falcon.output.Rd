\name{falcon.output}
\alias{falcon.output}
\title{
Quality control procedure for falcon's segmentation results
}
\description{
Applies a quality control procedure to the depth of coverage matrix both 
sample-wise and exon-wise before normalization.
}
\usage{
falcon.output(readMatrix, tauhat, cn, st_bp, end_bp, nboot)
}
\arguments{
    \item{readMatrix}{
    Allelic read depth matrix
    }
    \item{tauhat}{
    Change-points returned from falcon
    }
    \item{cn}{
    Segments from falcon
    }
    \item{st_bp}{
    Start base pair for each heterozygous loci
    }
    \item{end_bp}{
    End base pair for each heterozygous loci, usually the same as st_bp
    }
    \item{nboot}{
    Number of bootstrap samplings to get standard deviations of segmented copy numbers
    }
}
\value{
    \item{output}{
    Output table for falcon after QC with segment mean and standard deviation
    }
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\keyword{ package }