# MARATHON

Copy number variation is an important and abundant source of variation in the human genome, which has been associated with a number of diseases, especially cancer. Massively parallel next-generation sequencing allows copy number profiling with fine resolution. Such efforts, however, have met with mixed successes, with setbacks arising partly from the lack of reliable analytical methods to meet the diverse and unique challenges arising from the myriad experimental designs and study goals in genetic studies. In cancer genomics, detection of somatic copy number changes and profiling of allele-specific copy number (ASCN) are complicated by experimental biases and artifacts as well as normal cell contamination and cancer subclone admixture. Furthermore, careful statistical modeling is warranted to reconstruct tumor phylogeny by both somatic ASCN changes and single nucleotide variants. Here we describe a flexible computational pipeline, **MARATHON** (copy nu**M**ber v**AR**i**A**tion and **T**umor p**H**yl**O**ge**N**y), which integrates multiple related statistical software for copy number profiling and downstream analyses in disease genetic studies.

## Manuscript

The manuscript for MARATHON is currently under review. Pre-print is available [here via bioRixv](https://www.biorxiv.org/content/early/2017/09/28/195230).


## Questions & Problems

If you have any questions or problems when using MARATHON, please feel free to open a new issue [here](https://github.com/yuchaojiang/MARATHON/issues). You can also email the maintainers of the corresponding packages -- the contact information is shown under Developers & Maintainers.


## Installation

#### Installation Option 1: Docker Image - Good for ease of installation
A docker image is available [here](https://hub.docker.com/r/lzeppelini/marathon/).
This image is an Rstudio GUI built on rocker/tidyverse with MARATHON as well as all of its dependent packages pre-installed.  Also pre-installed is UCSC hg19 human reference genome and "WES.1KG.WUGSC" toy data set. Note that this can take a while to download the human reference genome as well as the toy sequencing dataset. Instructions for using Docker can be found [here](https://docs.docker.com/get-started/).

#### Installation Option 2: Install to R/RStudio - Good for performance
Install all packages in the latest version of [R](https://www.r-project.org/).
```r
install.packages("devtools")
library(devtools)
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")

biocLite("GenomeInfoDb")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("WES.1KG.WUGSC")
install_github("yuchaojiang/CODEX/package")

install_github("yuchaojiang/CODEX2/package")

install.packages("fields")
install.packages("truncnorm")
install.packages("ggplot2")
install_github("zhouzilu/iCNV")

install.packages("falcon")

install.packages("falconx")

install.packages("ape")
install.packages("pheatmap")
install.packages("scatterplot3d")
install_github("yuchaojiang/Canopy/package")
```


## Pipeline overview

The possible analysis scenarios are listed in Table 1. Figure 1 gives an outline for the relationship between the software: CODEX and CODEX2 perform read depth normalization for total copy number profiling; read depth normalized by CODEX/CODEX2 is received by iCNV, which combines it with allele-specific read counts and microarray data to detect CNVs; FALCON and FALCON-X perform ASCN analysis; and Canopy receives input from FALCON/FALCON-X to perform tumor phylogeny reconstruction.

**Figure 1.** A flowchart outlining the procedures for profiling CNV, ASCN, and reconstructing tumor phylogeny. CNVs with common and rare population frequencies can be profiled by CODEX and CODEX2, with and without negative control samples. iCNV integrates sequencing and microarray data for CNV detection. ASCNs can be profiled by FALCON and FALCON-X using allelic read counts at germline heterozygous loci. Canopy infers tumor phylogeny using somatic SNVs and ASCNs.

<p align="center">
  <img src='https://github.com/yuchaojiang/MARATHON/blob/master/figure/Figure1.jpg' width='600' height='300'>
</p>

**Table 1.** Analysis scenarios and pipeline design. The last column shows the sequencing of software that should be used for each analysis scenario. *By “normal” we mean samples that are not derived from tumor tissue. Broad copy number changes are not expected to cover the genome in the normal controls.

<p align="center">
  <img src='https://github.com/yuchaojiang/MARATHON/blob/master/figure/Table1.png' width='600' height='320'>
</p>


## Running MARATHON

* **Total copy number analysis of normal**
  * *Experimental design*: WES/WGS with or without matched microarray, greater than 20 samples
  * *Notebook*: [CODEX normalization](http://htmlpreview.github.io/?https://github.com/yuchaojiang/MARATHON/blob/master/notebook/codexNormalization/rNotebook.html);  [Rmd script](https://github.com/yuchaojiang/MARATHON/blob/master/notebook/codexNormalization/rNotebook.Rmd)
  * *Script*: [CODEX normalization](https://github.com/yuchaojiang/MARATHON/blob/master/script/CODEX.R), [CODEX --> iCNV](https://github.com/yuchaojiang/MARATHON/blob/master/script/CODEX_iCNV.R)
  * *More details*: [CODEX](https://github.com/yuchaojiang/CODEX), [iCNV](https://github.com/zhouzilu/iCNV)

* **Total copy number analysis of tumor**
  * *Experimental design*: WES/WGS/targeted sequencing of tumor with greater than 20 normal controls (no need to be matched)
  * *Notebook*: [CODEX2 normalization](http://htmlpreview.github.io/?https://github.com/yuchaojiang/MARATHON/blob/master/notebook/CODEX2/CODEX2.html);  [Rmd script](https://github.com/yuchaojiang/MARATHON/blob/master/notebook/CODEX2/CODEX2.Rmd)
  * *Script*: [CODEX2 normalization and segmentation](https://github.com/yuchaojiang/MARATHON/blob/master/script/CODEX2.R)
  * *More details*: [CODEX2](https://github.com/yuchaojiang/CODEX2)

* **Tumor allele-specific copy number by WGS**
  * *Experimental design*: WGS of matched tumor normal pair
  * *Notebook*: [FALCON for ASCN detection](http://htmlpreview.github.io/?https://github.com/yuchaojiang/MARATHON/blob/master/notebook/falcon/falconNotebook.html)
  * *Script*: [Profiling germline heterozygous loci](https://github.com/yuchaojiang/MARATHON/blob/master/script/germline_het_loci.sh), [FALCON for ASCN detection](https://github.com/yuchaojiang/MARATHON/blob/master/script/FALCON.R)
  * *More details*: [FALCON](https://CRAN.R-project.org/package=falcon)

* **Tumor allele-specific copy number by WES**
  * *Experimental design*: WES/WGS of tumor with matched normal, one or more tumor samples, greater than 20 normal controls (no need to be matched)
  * *Notebook*: [CODEX2 --> FALCON-X](http://htmlpreview.github.io/?https://github.com/yuchaojiang/MARATHON/blob/master/notebook/falconX/falconXNotebook.html)
  * *Script*: [Profiling germline heterozygous loci](https://github.com/yuchaojiang/MARATHON/blob/master/script/germline_het_loci.sh), [CODEX2 --> FALCON-X](https://github.com/yuchaojiang/MARATHON/blob/master/script/CODEX2_FALCONX.R)
  * *More details*: [CODEX2](https://github.com/yuchaojiang/CODEX2), [FALCON-X](https://CRAN.R-project.org/package=falconx)

* **Tumor phylogeny analysis by WGS**
  * *Experimental design*: WGS of multiple spatially or temporally separated tumor samples with matched normal
  * *Script*: [Profiling somatic point mutations](https://github.com/yuchaojiang/MARATHON/blob/master/script/somatic_SNV.sh), [FALCON -- > Canopy](https://github.com/yuchaojiang/MARATHON/blob/master/script/FALCON_Canopy.R)
  * *More details*: [FALCON](https://CRAN.R-project.org/package=falcon), [Canopy](https://github.com/yuchaojiang/Canopy)

* **Tumor phylogeny analysis by WES**
  * *Experimental design*: WES/WGS of multiple spatially or temporally separated tumor samples with matched normal, greater than 20 normal controls (no need to be matched)
  * *Script*: [Profiling somatic point mutations](https://github.com/yuchaojiang/MARATHON/blob/master/script/somatic_SNV.sh), [CODEX2 --> FALCON-X](https://github.com/yuchaojiang/MARATHON/blob/master/script/CODEX2_FALCONX.R), [Canopy](https://github.com/yuchaojiang/MARATHON/blob/master/script/Canopy.R)
  * *More details*: [CODEX2](https://github.com/yuchaojiang/CODEX2), [FALCON-X](https://CRAN.R-project.org/package=falconx), [Canopy](https://github.com/yuchaojiang/Canopy)



## Citation

Please cite MARATHON as well as all the dependent packages that you use.

* **MARATHON**: Integrative pipeline for profiling DNA copy number and inferring tumor phylogeny ([GitHub](https://github.com/yuchaojiang/MARATHON))
  <br>
  [Jiang et al. 2017 bioRixv](https://www.biorxiv.org/content/early/2017/09/28/195230)

* **CODEX**: A Normalization and Copy Number Variation Detection Method for Whole Exome Sequencing
  ([Bioconductor](http://bioconductor.org/packages/CODEX/), [GitHub](https://github.com/yuchaojiang/CODEX))
  <br>
  [Jiang et al. 2015 Nucleic Acids Research](https://academic.oup.com/nar/article/43/6/e39/2453417/CODEX-a-normalization-and-copy-number-variation)

* **CODEX2**: Full-spectrum copy number variation detection by high-throughput DNA sequencing
  ([GitHub](https://github.com/yuchaojiang/CODEX2))
  <br>
  [Jiang et al. 2017 bioRixv](https://www.biorxiv.org/content/early/2017/10/30/211698)

* **iCNV**: Integrated copy number variation detection toolset
  ([GitHub](https://github.com/zhouzilu/iCNV))
  <br>
  [Zhou et al. 2017 bioRixv](https://www.biorxiv.org/content/early/2017/09/01/172700)

* **FALCON**: Finding Allele-Specific Copy Number in Next-Generation Sequencing Data
  ([CRAN](https://CRAN.R-project.org/package=falcon))
  <br>
  [Chen et al. 2015 Nucleic Acids Research](https://academic.oup.com/nar/article/43/4/e23/2410993/Allele-specific-copy-number-profiling-by-next)

* **FALCON-X**: Finding Allele-Specific Copy Number in Whole-Exome Sequencing Data
  ([CRAN](https://CRAN.R-project.org/package=falconx))
  <br>
  [Chen et al. 2017 Annals of Applied Statistics](https://projecteuclid.org/euclid.aoas/1500537739)

* **Canopy**: Accessing Intra-Tumor Heterogeneity and Tracking Longitudinal and Spatial Clonal Evolutionary History by Next-Generation Sequencing
  ([CRAN](https://CRAN.R-project.org/package=Canopy), [GitHub](https://github.com/yuchaojiang/Canopy))
  <br>
  [Jiang et al. 2016 PNAS](http://www.pnas.org/content/113/37/E5528.full)


## Developers & Maintainers

* [Yuchao Jiang](http://sph.unc.edu/adv_profile/yuchao-jiang-phd/) (yuchaoj at email dot unc dot edu)
  <br>
  Department of Biostatistics & Department of Genetics, UNC-Chapel Hill
  
* [Hao Chen](https://anson.ucdavis.edu/~haochen/) (hxchen at ucdavis dot edu)
  <br>
  Department of Statistics, UC Davis

* Gene Urrutia (urrutia at email dot unc dot edu)
  <br>
  Department of Biostatistics, UNC-Chapel Hill

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, UPenn

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) (nzh at wharton dot upenn dot edu)
  <br>
  Department of Statistics, UPenn
