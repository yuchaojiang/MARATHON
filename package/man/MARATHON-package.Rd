\name{MARATHON-package}
\alias{MARATHON-package}
\alias{MARATHON}
\docType{package}
\title{
Integrative pipeline for profiling DNA copy number and inferring tumor phylogeny
}
\description{
Copy number variation is an important and abundant source of variation in the human genome, which has been associated with a number of diseases, especially cancer. Massively parallel next-generation sequencing allows copy number profiling with fine resolution. Such efforts, however, have met with mixed successes, with setbacks arising partly from the lack of reliable analytical methods to meet the diverse and unique challenges arising from the myriad experimental designs and study goals in genetic studies. In cancer genomics, detection of somatic copy number changes and profiling of allele-specific copy number (ASCN) are complicated by experimental biases and artifacts as well as normal cell contamination and cancer subclone admixture. Furthermore, careful statistical modeling is warranted to reconstruct tumor phylogeny by both somatic ASCN changes and single nucleotide variants. Here we describe a flexible computational pipeline, MARATHON (copy nuMber vARiAtion and Tumor pHylOgeNy), which integrates multiple related statistical software for copy number profiling and downstream analyses in disease genetic studies.
}
\details{
\tabular{ll}{
Package: \tab MARATHON\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2017-12-18\cr
License: \tab GPL-2\cr
}
Figure 1 shows a flowchart outlining the procedures for profiling CNV, ASCN, and reconstructing tumor phylogeny. CNVs with common and rare population frequencies can be profiled by CODEX and CODEX2, with and without negative control samples. iCNV integrates sequencing and microarray data for CNV detection. ASCNs can be profiled by FALCON and FALCON-X using allelic read counts at germline heterozygous loci. Canopy infers tumor phylogeny using somatic SNVs and ASCNs.
}
\author{
Yuchao Jiang <yuchaoj@email.unc.edu>
}
\keyword{ package }