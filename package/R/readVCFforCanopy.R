
readVCFforCanopy = function(vcfFile){
  
  if( !file.exists(vcfFile) ){
    print("The VCF file does not exist.  Please ensure the path is correct.")
    return(NULL)
  }
  
  vcfData = try ( VariantAnnotation::readVcf(vcfFile) )
  if ( class(vcfData) == "try-error" ){
    print ( "The VCF exists but could not be parsed.  Ensure the file is in proper VCF format.  See VariantAnnotation::readVcf for further information." )
    return( NULL )
  }
  
  # Extract total read depth
  X = geno(vcfData)$DP
  
  # Extract allele information
  vcfTargets = rowRanges(vcfData)
  
  # Extract mutant allele count from Allele Depth
  AD = geno(vcfData)$AD
  RVector = unlist(lapply(AD, function(x) x[[2]]))
  R = matrix(RVector,nrow(AD), ncol(AD))
  dimnames(R) = dimnames(AD)
  
  # Convert to Genotype Matrix
  #lapply(geno(vcfData), head)
  GT = geno(vcfData)$GT
  genotypeMatrix = matrix(0, nrow(GT), ncol(GT))
  genotypeMatrix[GT %in% c("0/0", "0|0")] = 0
  genotypeMatrix[GT %in% c("0/1", "0|1", "1/0", "1|0")] = 1
  genotypeMatrix[GT %in% c("1/1", "1|1")] = 2
  dimnames(genotypeMatrix) = dimnames(GT)
  
  # Convert to Reads Matrix
  vcfSampleNames = colnames(AD)
  readsMatrixList = lapply(vcfSampleNames, function(mySample) {
    RefVector = unlist(lapply(AD[, mySample], function(x) x[[1]]))
    AltVector = unlist(lapply(AD[, mySample], function(x) x[[2]]))
    readsMatrixSub = cbind(RefVector, AltVector)
    colnames(readsMatrixSub) = paste(mySample, c("Ref", "Alt"), sep = "_")
    return( readsMatrixSub ) 
  } )
  readsMatrix = do.call( cbind, readsMatrixList )
                           
  canopyInput = list(R = R,
                     X = geno(vcfData)$DP,
                     vcfTargets = rowRanges(vcfData),
                     GenotypeMatrix = genotypeMatrix,
                     ReadsMatrix = readsMatrix
  )
  return( canopyInput ) 
  
}


readVCFDataFrameforCanopy = function( vcf ) {
  
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
  
  AD.index = rep(0,nSNPs)
  GT.index = rep(0,nSNPs)
  DP.index = rep(0,nSNPs)
  
  for (i in 1:nSNPs){
    temp = strsplit(format[i,],":")[[1]]
    AD.index[i] = which(temp=="AD")
    GT.index[i] = which(temp=='GT')
    DP.index[i] = which(temp=='DP')
  }
  
  genotype = reads = matrix( NA, nSNPs, nSamples * 2 )
  
  for(i in 1:nSNPs){
    for(j in 1:nSamples){
      jj = j + 9
      temp = strsplit(as.matrix(vcf[i,jj]),":")[[1]]
      tempGT = strsplit(temp[GT.index[i]],"/")[[1]]
      tempAD = strsplit(temp[AD.index[i]],",")[[1]]
      tempDP = temp[DP.index[i]]
      genotype.temp=c(as.numeric(tempGT[1]), as.numeric(tempGT[2]))
      genotype[i,(2*j-1):(2*j)] = genotype.temp
      if(sum(!is.na(genotype.temp))==2){
        if(sum(genotype.temp==0)==2){
          reads[i,(2*j-1):(2*j)]=c(as.numeric(tempDP),0)
        } else if(sum(genotype.temp==1)==2){
          reads[i,(2*j-1):(2*j)]=c(0,as.numeric(tempDP))
        } else{
          reads[i,(2*j-1):(2*j)] = as.numeric(tempAD[genotype.temp+1])
        }
      }
    }
    if (i%%1000==0) cat(i,'\t')
  }
  
  mutantIndex = which( (1:ncol(reads) %% 2) == 0 )
  R = reads[, mutantIndex]  #mutant read count
  colnames(R) = sampleNames
  rownames(R) = vcf$ID
  X = reads[, mutantIndex] + reads[, mutantIndex - 1] #total read count
  colnames(X) = sampleNames
  rownames(X) = vcf$ID
  
  out = list(vcfTargets = vcf[,1:5],
             R = R,
             X = X)
  out 
}


