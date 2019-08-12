\name{readVCFforFalcon}
\alias{readVCFforFalcon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read and convert VCF file to usable format for falcon 
}
\description{
Read and convert VCF file to usable format for falcon
}
\usage{
readVCFforFalcon( myVCF )
}
\arguments{
  \item{myFileName}{
    path to VCF file
  }
}
\details{
This function takes as input a the path to the VCF file.  In the VCF file, First 9 column names following the header should be c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT").  Following column names should be sample names. see \href{http://www.internationalgenome.org/wiki/Analysis/vcf4.0/}{VCF format}
}

\value{
outputs a dataframe suitable for falcon :

Each row contains:

Overall info per variant:

- Chromosome 

- Start_position 

- End_position

And for each sample: 

- ReadCount_Alt: alternate allele count from the AD section of the VCF file

- ReadCount_Ref: reference allele count from the AD section of the VCF file

- ReadCount_Total: total depth from the DP section of the VCF file

- Allele1: first allele (nucleotides) from the GT section of the VCF file

- Allele2: second allele (nucleotides) from the GT section of the VCF file
}

\author{
  Gene Urrutia \email{gene.urrutia@gmail.com} 
}

\examples{
library(MARATHON)
vcfFile = system.file("extdata", "sample_w_header.vcf", package="MARATHON")
coverageData = readVCFforFalcon(vcfFile)

# Compare to relapse.demo from falcon
relapse.demo[1,]

# note that readMatrix must still be generated from the columns of interest.  In the relapse.demo, this would look like below where coverageChrHet is the post QC, chromosome specific version of coverageData

# not run
# readMatrix=as.data.frame(coverageChrHet[,c('Tumor_ReadCount_Ref',
#                                                   'Tumor_ReadCount_Alt',
#                                                   'Normal_ReadCount_Ref',
#                                                   'Normal_ReadCount_Alt')])

}
\keyword{ package }

