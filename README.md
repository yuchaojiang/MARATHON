# MARATHON

Copy number variation is an important and abundant source of variation in the human genome, which has been associated with a number of diseases, especially cancer. Massively parallel next-generation sequencing allows copy number profiling with fine resolution. Such efforts, however, have met with mixed successes, with setbacks arising partly from the lack of reliable analytical methods to meet the diverse and unique challenges arising from the myriad experimental designs and study goals in genetic studies. In cancer genomics, detection of somatic copy number changes and profiling of allele-specific copy number (ASCN) are complicated by experimental biases and artifacts as well as normal cell contamination and cancer subclone admixture. Furthermore, careful statistical modeling is warranted to reconstruct tumor phylogeny by both somatic ASCN changes and single nucleotide variants. Here we describe a flexible computational pipeline, **MARATHON** (copy nu**M**ber v**AR**i**A**tion and **T**umor p**H**yl**O**ge**N**y), which integrates multiple related statistical software for copy number profiling and downstream analyses in disease genetic studies.

## Manuscript

Urrutia E, Chen H, Zhou Z, Zhang NR, Jiang Y. Integrative pipeline for profiling DNA copy number and inferring tumor phylogeny. *Bioinformatics*, bty057, 2018. ([link](https://doi.org/10.1093/bioinformatics/bty057))


## Questions & Problems

If you have any questions or problems when using MARATHON, please feel free to open a new issue [here](https://github.com/yuchaojiang/MARATHON/issues). You can also email the maintainers of the corresponding packages -- the contact information is shown under Developers & Maintainers.


## Installation

#### Installation Option 1: Docker Image - Good for ease of installation
A docker image is available [here](https://hub.docker.com/r/lzeppelini/marathon/).
This image is an Rstudio GUI built on rocker/tidyverse with MARATHON as well as all of its dependent packages and datasets pre-installed. Note that this can take a while to download the human reference genome as well as the toy sequencing dataset. Instructions for using Docker can be found [here](https://docs.docker.com/get-started/).

```bash
docker pull lzeppelini/marathon
```

#### Installation Option 2: Install to R/RStudio - Good for performance
Install all packages in the latest version of [R](https://www.r-project.org/).
```r
install.packages(c("falcon", "falconx", "devtools"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("WES.1KG.WUGSC", version = "3.8")
devtools::install_github(c("yuchaojiang/CODEX/package", "yuchaojiang/CODEX2/package", "yuchaojiang/Canopy/package", "zhouzilu/iCNV", "yuchaojiang/MARATHON/package"))
```

## Pipeline overview

The possible analysis scenarios are listed in Table 1. Figure 1 gives an outline for the relationship between the software: CODEX and CODEX2 perform read depth normalization for total copy number profiling; read depth normalized by CODEX/CODEX2 is received by iCNV, which combines it with allele-specific read counts and microarray data to detect CNVs; FALCON and FALCON-X perform ASCN analysis; and Canopy receives input from FALCON/FALCON-X to perform tumor phylogeny reconstruction.

**Figure 1.** A flowchart outlining the procedures for profiling CNV, ASCN, and reconstructing tumor phylogeny. CNVs with common and rare population frequencies can be profiled by CODEX and CODEX2, with and without negative control samples. iCNV integrates sequencing and microarray data for CNV detection. ASCNs can be profiled by FALCON and FALCON-X using allelic read counts at germline heterozygous loci. Canopy infers tumor phylogeny using somatic SNVs and ASCNs.

<p align="center">
  <img src='https://github.com/yuchaojiang/MARATHON/blob/master/figure/Figure1.jpg' width='600' height='300'>
</p>

**Table 1.** Analysis scenarios and pipeline design. The last column shows the sequence of software that should be used for each analysis scenario. * By “normal” we mean samples that are not derived from tumor tissue, which are not expected to carry chromosome-level copy number changes.

<p align="center">
  <img src='https://github.com/yuchaojiang/MARATHON/blob/master/figure/Table1.png' width='600' height='320'>
</p>


## Running MARATHON

**R notebook** with step-by-step demonstration and rich display is available [***here***](https://rawgit.com/yuchaojiang/MARATHON/master/notebook/MARATHON.html). Corresponding **Rmd script** is available [***here***](https://github.com/yuchaojiang/MARATHON/blob/master/notebook/MARATHON.Rmd).


## Citation

Please cite MARATHON as well as all the dependent packages that you use.

* **MARATHON**: [Urrutia et al. 2018 Bioinformatics](https://doi.org/10.1093/bioinformatics/bty057)
  <br>
  Integrative pipeline for profiling DNA copy number and inferring tumor phylogeny ([GitHub](https://github.com/yuchaojiang/MARATHON))

* **CODEX**: [Jiang et al. 2015 Nucleic Acids Research](https://academic.oup.com/nar/article/43/6/e39/2453417/CODEX-a-normalization-and-copy-number-variation)
  <br>
  A Normalization and Copy Number Variation Detection Method for Whole Exome Sequencing
  ([Bioconductor](http://bioconductor.org/packages/CODEX/), [GitHub](https://github.com/yuchaojiang/CODEX))

* **CODEX2**: [Jiang et al. 2018 bioRixv](https://www.biorxiv.org/content/early/2018/04/01/211698)
  <br>
  Full-spectrum copy number variation detection by high-throughput DNA sequencing
  ([GitHub](https://github.com/yuchaojiang/CODEX2))

* **iCNV**: [Zhou et al. 2017 Bioinformatics](https://doi.org/10.1093/bioinformatics/bty104)
  <br>
  Integrated copy number variation detection toolset
  ([GitHub](https://github.com/zhouzilu/iCNV))

* **FALCON**: [Chen et al. 2015 Nucleic Acids Research](https://academic.oup.com/nar/article/43/4/e23/2410993/Allele-specific-copy-number-profiling-by-next)
  <br>
  Finding Allele-Specific Copy Number in Next-Generation Sequencing Data
  ([CRAN](https://CRAN.R-project.org/package=falcon))

* **FALCON-X**: [Chen et al. 2017 Annals of Applied Statistics](https://projecteuclid.org/euclid.aoas/1500537739)
  <br>
  Finding Allele-Specific Copy Number in Whole-Exome Sequencing Data
  ([CRAN](https://CRAN.R-project.org/package=falconx))

* **Canopy**: [Jiang et al. 2016 PNAS](http://www.pnas.org/content/113/37/E5528.full)
  <br>
  Accessing Intra-Tumor Heterogeneity and Tracking Longitudinal and Spatial Clonal Evolutionary History by Next-Generation Sequencing
  ([CRAN](https://CRAN.R-project.org/package=Canopy), [GitHub](https://github.com/yuchaojiang/Canopy))

## Developers & Maintainers

* Gene Urrutia (urrutia at email dot unc dot edu)
  <br>
  Department of Biostatistics, UNC-Chapel Hill

* [Yuchao Jiang](http://sph.unc.edu/adv_profile/yuchao-jiang-phd/) (yuchaoj at email dot unc dot edu)
  <br>
  Department of Biostatistics & Department of Genetics, UNC-Chapel Hill

* [Hao Chen](https://anson.ucdavis.edu/~haochen/) (hxchen at ucdavis dot edu)
  <br>
  Department of Statistics, UC Davis

* Zilu Zhou (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, UPenn

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) (nzh at wharton dot upenn dot edu)
  <br>
  Department of Statistics, UPenn
