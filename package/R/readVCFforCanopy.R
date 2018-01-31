readVCFforCanopy = function( vcf ) {
  
  correctNames = c("CHROM", "POS", "ID", "REF", 
                   "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  
  if ( !all( names(vcf)[1:9] == correctNames ) ){
    stop(
      paste0(
        "Please ensure the input file is in VCF format with first 9 column names: ", 
        paste(correctNames, collapse = ", "), 
        ", followed by one column per sample"
      )
    )
  }
  
  sampleNames = colnames(vcf)[10:ncol(vcf)]
  format = as.matrix(vcf[,'FORMAT'])
  nSNPs = nrow( vcf )
  nSamples = ncol( vcf ) - 9 
  
  AD.index.tumor = rep(0,nSNPs)
  GT.index.tumor = rep(0,nSNPs)
  DP.index.tumor = rep(0,nSNPs)
  
  for (i in 1:nSNPs){
    temp = strsplit(format[i,],":")[[1]]
    AD.index.tumor[i] = which(temp=="AD")
    GT.index.tumor[i] = which(temp=='GT')
    DP.index.tumor[i] = which(temp=='DP')
  }
  
  genotype.tumor = reads.tumor = matrix( NA, nSNPs, nSamples * 2 )
  
  for(i in 1:nSNPs){
    for(j in 1:nSamples){
      jj = j + 9
      temp = strsplit(as.matrix(vcf[i,jj]),":")[[1]]
      tempGT = strsplit(temp[GT.index.tumor[i]],"/")[[1]]
      tempAD = strsplit(temp[AD.index.tumor[i]],",")[[1]]
      tempDP = temp[DP.index.tumor[i]]
      genotype.temp=c(as.numeric(tempGT[1]), as.numeric(tempGT[2]))
      genotype.tumor[i,(2*j-1):(2*j)] = genotype.temp
      if(sum(!is.na(genotype.temp))==2){
        if(sum(genotype.temp==0)==2){
          reads.tumor[i,(2*j-1):(2*j)]=c(as.numeric(tempDP),0)
        } else if(sum(genotype.temp==1)==2){
          reads.tumor[i,(2*j-1):(2*j)]=c(0,as.numeric(tempDP))
        } else{
          reads.tumor[i,(2*j-1):(2*j)] = as.numeric(tempAD[genotype.temp+1])
        }
      }
    }
    if (i%%1000==0) cat(i,'\t')
  }
  
  mutantIndex = which( (1:ncol(reads.tumor) %% 2) == 0 )
  R = reads.tumor[, mutantIndex]  #mutant read count
  colnames(R) = sampleNames
  rownames(R) = vcf$ID
  X = reads.tumor[, mutantIndex] + reads.tumor[, mutantIndex - 1] #total read count
  colnames(X) = sampleNames
  rownames(X) = vcf$ID
  
  out = list(vcfTargets = vcf[,1:5],
             R = R,
             X = X)
  out 
}
