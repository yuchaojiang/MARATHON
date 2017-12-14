---
    title: "Falcon-X Notebook"
    author: "Gene Urrutia"  
    date: "2017-12-14"
    output:  
      html_document:  
        keep_md: true 
---  

#Load Packages

```r
library(CODEX, quietly = T) 
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## 
## Attaching package: 'CODEX'
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     normalize
```

```r
library(falconx, quietly = T) 
```

#source additional falcon functions

```r
faclonURL = "https://github.com/yuchaojiang/Canopy/blob/master/instruction"
source(paste0(faclonURL, "/falconx.qc.R?raw=True"))
```

#  Below is demo dataset consisting of 39 tumor-normal paired whole-exome sequencing, published in [Maxwell et al.](https://www.nature.com/articles/s41467-017-00388-9) (Nature Communications, 2017).  We focus on chr17, where copy-neutral loss-of-heterozygosity has been reported.

# Below are allelic reads, genotype, and genomic locations, which can be extracted from vcf files. rda files available for download [here] (https://github.com/yuchaojiang/Canopy/tree/master/instruction)


```r
chr=17 
setwd("/home/rstudio/kitematic/notebook/falconX")
load('mymatrix.demo.rda')
load('genotype.demo.rda')
load('reads.demo.rda')
```

genomic locations and SNP info across all loci from chr17

```r
head(mymatrix) 
```

```
##      CHROM POS     REF ALT
## [1,] "17"  "6115"  "G" "C"
## [2,] "17"  "6157"  "A" "G"
## [3,] "17"  "6188"  "C" "T"
## [4,] "17"  "11830" "A" "G"
## [5,] "17"  "11837" "G" "A"
## [6,] "17"  "11843" "G" "A"
```

genotype in blood (bloodGT1 and blood GT2) and tumor (tumorGT1 and tumor GT2) across all samples

```r
head(genotype[,1:12]) 
```

```
##      samp1_bloodGT1 samp1_bloodGT2 samp1_tumorGT1 samp1_tumorGT2
## [1,]              0              1              0              1
## [2,]              0              1              0              1
## [3,]              0              1              0              1
## [4,]              1              1              1              1
## [5,]              0              0              0              0
## [6,]              0              0              0              0
##      samp2_bloodGT1 samp2_bloodGT2 samp2_tumorGT1 samp2_tumorGT2
## [1,]              0              0              0              1
## [2,]              0              0              0              1
## [3,]              0              1              0              1
## [4,]              0              1              0              1
## [5,]              0              0              0              0
## [6,]              0              0              0              0
##      samp3_bloodGT1 samp3_bloodGT2 samp3_tumorGT1 samp3_tumorGT2
## [1,]              1              1              1              1
## [2,]              1              1              1              1
## [3,]              0              0              0              0
## [4,]              1              1              1              1
## [5,]              0              0              0              0
## [6,]              0              0              0              0
```

allelic reads in blood (AN and BN) and tumor (AT and BT) across all samples

```r
head(reads[,1:12]) 
```

```
##      samp1_AN samp1_BN samp1_AT samp1_BT samp2_AN samp2_BN samp2_AT
## [1,]      123       78      122      128      240        8      191
## [2,]       80       54       83       85      194       11      170
## [3,]       35       49       43       49       76       50       80
## [4,]        0       33        1       37       20       26       24
## [5,]       37        0       50        0       49        0       70
## [6,]       39        0       55        0       55        0       79
##      samp2_BT samp3_AN samp3_BN samp3_AT samp3_BT
## [1,]       58        0      158        0       72
## [2,]       59        0      105        0       41
## [3,]       58       82        0       25        0
## [4,]       37        0       57        0       16
## [5,]        0       70        0       16        0
## [6,]        0       71        0       22        0
```

#####################################################################################
# Apply CODEX2 to get total coverage bias
#####################################################################################

# Get GC content from a 50bp window centered at the SNP

```r
pos=as.numeric(mymatrix[,'POS'])
ref=IRanges(start=pos-25,end=pos+25)
gc=getgc(chr,ref)  
```

# total read depth

```r
Y=matrix(nrow=nrow(reads),ncol=ncol(reads)/2)  
for(j in 1:ncol(Y)){
  Y[,j]=reads[,2*j-1]+reads[,2*j]
}
```

# QC procedure

```r
pos.filter=(apply(Y,1,median)>=20)
pos=pos[pos.filter]
ref=ref[pos.filter]
gc=gc[pos.filter]
genotype=genotype[pos.filter,]
reads=reads[pos.filter,]
mymatrix=mymatrix[pos.filter,]
Y=Y[pos.filter,]

normObj=normalize2(Y,gc,K=1:3,normal_index=seq(1,77,2))
choiceofK(normObj$AIC,normObj$BIC,normObj$RSS,K=1:3,filename=paste('choiceofK.',chr,'.pdf',sep=''))
```

```
## png 
##   2
```

```r
cat(paste('BIC is maximized at ',which.max(normObj$BIC),'.\n',sep=''))
```

```
## BIC is maximized at 3.
```

```r
Yhat=round(normObj$Yhat[[3]],0)
```



```r
dim(mymatrix)  # raw vcf read.table, including genomic locations
```

```
## [1] 9673    4
```

```r
dim(reads)   # allelic reads
```

```
## [1] 9673  156
```

```r
dim(genotype)   # genotype
```

```
## [1] 9673  156
```

```r
dim(Yhat)  # total coverage bias returned by CODEX
```

```
## [1] 9673   78
```


#################################################################################
# Generate input for FALCON-X: allelic read depth and genotype across all loci
#################################################################################


```r
n=39 # total number of samples
for (i in 1:n){
  cat('Generating input for sample',i,'...\n')
  ids = (4*i-3):(4*i)
  ids2 = (2*i-1):(2*i)
  mydata = as.data.frame(cbind(mymatrix[,1:2], genotype[,ids], reads[,ids], Yhat[,ids2]))
  colnames(mydata) = c("chr", "pos", "bloodGT1", "bloodGT2", "tumorGT1",
                       "tumorGT2", "AN", "BN", "AT", "BT", 'sN','sT')
  ids=which(as.numeric(mydata[,3])!=as.numeric(mydata[,4]))
  newdata0 = mydata[ids,]
  index.na=apply(is.na(newdata0), 1, any)
  newdata=newdata0[index.na==FALSE,]
  
  # Remove loci with multiple alternative alleles
  mul.alt.filter=rep(TRUE,nrow(newdata))
  for(s in 1:nrow(newdata)){
    filter1=!is.element(as.numeric(newdata[s,'tumorGT1']),
                        c(as.numeric(newdata[s,'bloodGT1']),
                          as.numeric(newdata[s,'bloodGT2'])))
    filter2=!is.element(as.numeric(newdata[s,'tumorGT2']),
                        c(as.numeric(newdata[s,'bloodGT1']),
                          as.numeric(newdata[s,'bloodGT2'])))
    if(filter1 | filter2){
      mul.alt.filter[s]=FALSE
    }
  }
  newdata=newdata[mul.alt.filter,]
  
  # write text at germline heterozygous loci, which is used as input for Falcon-X
  write.table(newdata, file=paste("sample",i,"_het.txt",sep=""), quote=F, row.names=F)
}
```

```
## Generating input for sample 1 ...
## Generating input for sample 2 ...
## Generating input for sample 3 ...
## Generating input for sample 4 ...
## Generating input for sample 5 ...
## Generating input for sample 6 ...
## Generating input for sample 7 ...
## Generating input for sample 8 ...
## Generating input for sample 9 ...
## Generating input for sample 10 ...
## Generating input for sample 11 ...
## Generating input for sample 12 ...
## Generating input for sample 13 ...
## Generating input for sample 14 ...
## Generating input for sample 15 ...
## Generating input for sample 16 ...
## Generating input for sample 17 ...
## Generating input for sample 18 ...
## Generating input for sample 19 ...
## Generating input for sample 20 ...
## Generating input for sample 21 ...
## Generating input for sample 22 ...
## Generating input for sample 23 ...
## Generating input for sample 24 ...
## Generating input for sample 25 ...
## Generating input for sample 26 ...
## Generating input for sample 27 ...
## Generating input for sample 28 ...
## Generating input for sample 29 ...
## Generating input for sample 30 ...
## Generating input for sample 31 ...
## Generating input for sample 32 ...
## Generating input for sample 33 ...
## Generating input for sample 34 ...
## Generating input for sample 35 ...
## Generating input for sample 36 ...
## Generating input for sample 37 ...
## Generating input for sample 38 ...
## Generating input for sample 39 ...
```


#################################################################################
# Apply FALCON-X to generate allele-specific copy number profiles
#################################################################################

# CODEX normalize total read depth across samples
# falcon-x profiles ASCN in each sample separately

```r
k=10 # calling ASCN for the 10th sample
ascn.input=read.table(paste("sample",k,"_het.txt",sep=""),head=T)
readMatrix=ascn.input[,c('AN','BN','AT','BT')]
biasMatrix=ascn.input[,c('sN','sT')]

tauhat = getChangepoints.x(readMatrix, biasMatrix, pos=ascn.input$pos)
```

```
## Number of loci: 1675 
## Scanning region between variants 1 to 1200 for change-points ... 
##   Candidate change-point(s): 21 81 91 151 341 351 551 611 621 661 671 911 941 1171 
##   Adding 21 81 91 151 341 351 551 611 621 661 671 911 941 1171 to current change-point list.
##   Current change-point list: 21 81 91 151 341 351 551 611 621 661 671 911 941 1171 
## Scanning region between variants 1001 to 1675 for change-points ... 
##   Candidate change-point(s): 1171 1201 1261 1521 1551 
##   Adding 1171 1201 1261 1521 1551 to current change-point list.
##   Current change-point list: 21 81 91 151 341 351 551 611 621 661 671 911 941 1171 1171 1201 1261 1521 1551 
## 
## Estimated change-points of the whole sequence: 21 81 91 151 341 351 551 611 621 661 671 911 941 1171 1201 1261 1521 1551
```

```r
cn = getASCN.x(readMatrix, biasMatrix, tauhat=tauhat, pos=ascn.input$pos, threshold = 0.3)
```


# cn$tauhat would give the indices of change-points.
# cn$ascn would give the estimated allele-specific copy numbers for each segment.
# cn$Haplotype[[i]] would give the estimated haplotype for the major chromosome in segment i
# if this segment has different copy numbers on the two homologous chromosomes.

```r
view(cn, pos=ascn.input$pos)
```

![](falconXNotebook_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

# Further curate FALCON-X’s segmentation:
# Remove small regions and combine consecutive regions with similar ASCN profiles.

```r
if(length(tauhat)>0){
  length.thres=10^6  # Threshold for length of segments, in base pair.
  delta.cn.thres=0.3  # Threshold of absolute copy number difference between consecutive segments.
  falcon.qc.list = falconx.qc(readMatrix = readMatrix,
                              biasMatrix = biasMatrix,
                              tauhat = tauhat,
                              cn = cn,
                              st_bp = ascn.input$pos,
                              end_bp = ascn.input$pos,
                              length.thres = length.thres,
                              delta.cn.thres = delta.cn.thres)
  
  tauhat=falcon.qc.list$tauhat
  cn=falcon.qc.list$cn
}

view(cn,pos=ascn.input$pos)
```

![](falconXNotebook_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

