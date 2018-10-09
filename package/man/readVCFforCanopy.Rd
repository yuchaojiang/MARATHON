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
    path to VCF file
  }
}
\details{
This function takes as input a the path to the VCF file.  First 9 column names following the header should be c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT").  Following column names should be sample names. see \href{http://www.internationalgenome.org/wiki/Analysis/vcf4.0/}{VCF format}
}

\value{
outputs a list of four objects:

R is a matrix (SNA by sample) which contains the number of alternate reads.  

X is a matrix (SNA by sample) which contains the number of total (referecnce + alternate) reads.

vcfTargets is the snp info from the vcf file.  This is useful for annotating X and R if desired.  

GenotypeMatrix is a snp by sample matrix representing alleles at that position (i.e. 0: homozygous ref, 1: heterozygous alt, 2: homozygous alt)

ReadsMatrix is a matrix where rows are positions and each sample has two columns, first the  reference allele count, and second the alternate allele count.

}

\author{
  Gene Urrutia \email{gene.urrutia@gmail.com} 
}

\examples{
library(MARATHON)

vcfFile = system.file("extdata", "sample_w_header.vcf", package="MARATHON")
canopyInput = readVCFforCanopy(vcfFile)
lapply(canopyInput, head)

}
\keyword{ package }

