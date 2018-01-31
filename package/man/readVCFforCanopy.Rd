\name{readVCFforCanopy}
\alias{readVCFforCanopy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read and convert VCF file to usable format for Canopy 
}
\description{
Read and convert VCF file to usable format for Canopy 
}
\usage{
readVCFforCanopy( myVCF )
}
\arguments{
  \item{myFileName}{
    vcf represented as a data frame
  }
}
\details{
This function takes as input a data frame representing the VCF file as read using read.lines().  First 9 column names should be c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT").  Following column names should be sample names. see \href{http://www.internationalgenome.org/wiki/Analysis/vcf4.0/}{VCF format}
}

\value{
outputs a list of three objects:

R is a matrix (SNA by sample) which contains the number of alternate reads.  Rownames of R (and X) are the "ID" of the VCF file

X is a matrix (SNA by sample) which contains the number of total (referecnce + alternate) reads

vcfTargets is target info from the vcf file.  This is useful for annotating X and R if desired.  }

\author{
 Yuchao Jiang \email{yuchaoj@live.unc.edu} 
 
 Gene Urrutia \email{urutia@live.unc.edu} 
}

\examples{
install.packages("/Users/urrutig1/Documents/gitRepos/MARATHON/package", repos = NULL, type = 'source')

library(MARATHON)

vcfFile = system.file("extdata", "sample.vcf", package="MARATHON")

myVCF = read.table(vcfFile, sep = "\t", header = T)

canopyInput = readVCFforCanopy( myVCF )

\keyword{ package }

