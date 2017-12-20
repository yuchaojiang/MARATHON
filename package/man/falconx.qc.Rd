\name{falconx.qc}
\alias{falconx.qc}
\title{
Quality control procedure for falconx's segmentation results
}
\description{
Applies a quality control procedure to the depth of coverage matrix both 
sample-wise and exon-wise before normalization.
}
\usage{
falconx.qc(readMatrix, biasMatrix, tauhat, cn, st_bp, end_bp,
           length.thres, delta.cn.thres)
}
\arguments{
    \item{readMatrix}{
    Allelic read depth matrix
    }
    \item{biasMatrix}{
    Bias matrix at each loci returned by CODEX2
    }
    \item{tauhat}{
    Change-points returned from falconx
    }
    \item{cn}{
    Segments from falconx
    }
    \item{st_bp}{
    Start base pair for each heterozygous loci
    }
    \item{end_bp}{
    End base pair for each heterozygous loci, usually the same as st_bp
    }
    \item{length.thres}{
    Threshold on length of falconx's segments. Set by default as 10^6 bp.
    }
    \item{delta.cn.thres}{
    Threshold on copy number change between consecutive segments. Set by default as 0.3.
    }
}
\value{
    \item{tauhat}{
    Updated change-points returned from falconx after QC.
    }
    \item{cn}{
    Updated segments from falconx after QC.
    }
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\keyword{ package }