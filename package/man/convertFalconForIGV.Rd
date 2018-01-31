\name{convertFalconForIGV}
\alias{convertFalconForIGV}
\title{
Convert Falcon CNV calls to .seg format for IGV
}
\description{
Takes as input the output of falcon.output from MARATHON package and converts to .seg format.  This is a preferred format for Integrative Genomics Viewer (IGV).  This function takes a single data frame as input, so combine CNV calls from all chromosomes beforehand if necessary.
}
\usage{
convertFalconForIGV( falconOutput, sampleName, outPath, filterFraction, log2 )
}
\arguments{
    \item{ finalcall }{
    Final callset of CNVs with genotyping results returned from 
    \code{\link{segment}} or \code{\link{segmentHMM}} }
    \item{ sampleName }{
    Sample name corresponding to falconOutput }
    \item{ outPath }{
    Directory to save to file.  If no directory is specified it is saved to the local directory.}
    \item{ log2 }{
    output as log2 cnv TRUE by default, otherwise output as raw}
}
\value{
IGV ready file(s) in .seg format is saved to specified directory.  One file is prepared each for minor, major, and total CNV.  Total is calculated as mean CNV of major and minor.  If log2 is FALSE, cvn is doubled to ensure color compatibility.  More detail at \href{https://software.broadinstitute.org/software/igv/CN}{IGV Documentation} 
}
\note{ 
This function takes a single data frame as input, so combine CNV calls from all chromosomes beforehand if necessary.

If a file containes only non-negative log2 values, color will not be displayed correctly.  More detail at \href{https://software.broadinstitute.org/software/igv/CN}{IGV Documentation}
}

\author{
 Yuchao Jiang \email{yuchaoj@email.unc.edu} 
 
 Gene Urrutia \email{urutia@email.unc.edu} 
}
\examples{
library(MARATHON)

falconOutput  #from the MARTHON notebook section 4.3.2. 
#Running FALCON for allele-specific copy number profiling

convertFalconForIGV( 
  falconOutput = falconOutput, sampleName = "sample1", outPath = getwd(), log2 = FALSE)
}
\keyword{ package }
