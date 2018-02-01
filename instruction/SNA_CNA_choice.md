## **Which CNAs and SNAs should I use?**

  How to generate a clean set of input to Canopy is important and non-trivial. While in our phylogeney reconstruction paper we were not trying to solve this tricky but separate problem of quality control of point mutation profiling and copy number estimation, an input with too many false positives will only lead to "garbage in garbage out" by Canopy. We are currently working on automating the pipeline for generating CNA and SNA input as well as offering guidance to select the ***informative*** SNAs and CNAs. By saying ***informative***, we mean that the SNAs or CNAs show distinct patterns between different samples (from the same patient since we are looking at intratumor heterogeneity). For SNAs, this means that the observed VAFs are different (see __*Figure 4B*__ in our paper) and in this case a heatmap is a good way for visualization. For CNAs, this means that the WM and Wm are different (see __*Supplementary Figure S13*__ in our paper) and we find **[IGV](http://software.broadinstitute.org/software/igv/)** a good tool for visualization and recommend focusing on large CNA regions, which helps remove false calls and speed up computation.

### CNV example
Below is a view of allele-specific copy number for 3 samples using Integrative Genomics Viewer.  We see major, minor and total copy number for each sample across the genome.   In this example, chromosomes 2, 8, 9, 11, 12 and 17b (second CNA) show differential CNV by sample, so we would use these.  Chromosomes 7, 15, and 17a (first CNA) show near identical CNV by sample, so we would not use these.
  <p align="center">
  IGV allele-specific copy number view of 3 samples
  <img src='https://github.com/yuchaojiang/MARATHON/blob/master/instruction/igvCNA.png' width='994' height='479'>
  </p>

### SNA example
Below is a heatmap of variant allele frequency by Gene and by sample.  This is a curated set which includes only VAF showing differential by sample.  These would be good candidates for input to Canopy, in particular the lower portion where the first 5 samples have low VAF, and the next 5 samples have high VAF.  The heatmap below is plotted using the [pheatmap R package](https://CRAN.R-project.org/package=pheatmap).
<p align="center">
SNA heatmap </br>
  <img src='https://github.com/yuchaojiang/MARATHON/blob/master/instruction/snaHeatmap.png' width='480' height='900'>
</p>

  Just like SNAs, there will likely to be CNAs carried by small fractions of the cells, or that reside in hard-to-call regions of the genome, which are not detected. The former scenario, in particular, includes CNAs which may be informative about rare subclones. We do not assume that the CNAs (and SNAs) given to Canopy comprise all mutations that are carried by the sample, and similarly, do not attempt to claim that Canopy detects all subclones or resolves all branches of the evolutionary tree. Our goal is only to estimate all of the subclones that have representation among the set of input CNAs and SNAs, which are, inherently, limited in resolution by the experimental protocol (sequencing platform and coverage, number of tumor slices, etc.) We believe this is the best that can be achieved under current settings.
