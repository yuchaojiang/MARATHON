

###################################################
### Pre-processing
###################################################
library(CODEX)
library(WES.1KG.WUGSC) # Load Toy data from the 1000 Genomes Project.
dirPath <- system.file("extdata", package = "WES.1KG.WUGSC")
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampname <- as.matrix(read.table(file.path(dirPath, "sampname")))
bedFile <- file.path(dirPath, "chr22_400_to_500.bed")
chr <- 22
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX_demo", chr)
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname; chr <- bambedObj$chr


###################################################
### Get read depth coverage
###################################################
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y; readlength <- coverageObj$readlength


###################################################
### Get GC content and mappability
###################################################
gc <- getgc(chr, ref)
mapp <- getmapp(chr, ref)


###################################################
### Quality control
###################################################
qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh = c(20, 4000), 
            length_thresh = c(20, 2000), mapp_thresh = 0.9, gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc; gc_qc <- qcObj$gc_qc
mapp_qc <- qcObj$mapp_qc; ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat
write.table(qcmat, file = paste(projectname, '_', chr, '_qcmat', '.txt', sep=''),
            sep='\t', quote=FALSE, row.names=FALSE)


###################################################
### Normalization
###################################################
normObj <- normalize(Y_qc, gc_qc, K = 1:9)
Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC
RSS <- normObj$RSS; K <- normObj$K



###################################################
### Choose the number of latent factors
###################################################
 choiceofK(AIC, BIC, RSS, K, filename = paste(projectname, "_", chr, 
     "_choiceofK", ".pdf", sep = ""))



###################################################
### Poisson-likelihood segmentation
###################################################
optK = K[which.max(BIC)]
finalcall <- segment(Y_qc, Yhat, optK = optK, K = K, sampname_qc,
                     ref_qc, chr, lmax = 200, mode = "integer")
head(finalcall)
write.table(finalcall, file = paste(projectname, '_', chr, '_', optK,
             '_CODEX_frac.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE)
save.image(file = paste(projectname, '_', chr, '_image', '.rda', sep=''),
      compress='xz')

