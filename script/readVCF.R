install.packages(pkgs = "~/Documents/gitRepos/MARATHON/package",
                 repos = NULL,
                 type="source")
library(MARATHON)
library(VariantAnnotation)

vcfFile = system.file("extdata", "sample_w_heade.vcf", package="MARATHON")
canopyInput = readVCF(vcfFile)
#error file does not exist

vcfFile = system.file("extdata", "sample.vcf", package="MARATHON")
canopyInput = readVCF(vcfFile)
#error file not in correct VCF format

vcfFile = system.file("extdata", "sample_w_header.vcf", package="MARATHON")
canopyInput = readVCF(vcfFile)
lapply(canopyInput, head)


