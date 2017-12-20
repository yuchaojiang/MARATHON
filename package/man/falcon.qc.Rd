\name{falcon.qc}
\alias{falcon.qc}
\title{
Quality control procedure for falcon's segmentation results
}
\description{
Applies a quality control procedure to the depth of coverage matrix both 
sample-wise and exon-wise before normalization.
}
\usage{
falcon.qc(readMatrix, tauhat, cn, st_bp, end_bp, 
          rdep, length.thres, delta.cn.thres)
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
    \item{rdep}{
    Read depth factor. Set by default as ratio of total number of reads from tumor versus normal.
    }
    \item{length.thres}{
    Threshold on length of falcon's segments. Set by default as 10^6 bp.
    }
    \item{delta.cn.thres}{
    Threshold on copy number change between consecutive segments. Set by default as 0.3.
    }
}
\value{
    \item{tauhat}{
    Updated change-points returned from falcon after QC.
    }
    \item{cn}{
    Updated segments from falcon after QC.
    }
}
\author{
Yuchao Jiang \email{yuchaoj@email.unc.edu}
}
\keyword{ package }