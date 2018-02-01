convertFalconForIGV = function( 
  falconOutput, sampleName, outPath = NULL, log2 = TRUE){
  
  for (type in c( 'minor', 'major', 'total' ) ){
    sampleNameTemp = paste( sampleName, 'falcon', type, sep='.' )
    
    # multipling by 2 ensures color compatibiltity with IGV 
    
    if ( type == 'minor' ) copyNo =  2 * falconOutput$Minor_copy 
    if ( type == 'major' ) copyNo =  2 * falconOutput$Major_copy 
    if ( type == 'total' ) copyNo =  ( falconOutput$Minor_copy + 
                                              falconOutput$Major_copy ) 
    if ( log2 ) { copyNo = log2( copyNo ) } 
     
    output = cbind( sampleNameTemp, falconOutput$chr,
                    falconOutput$st_bp, falconOutput$end_bp,
                    falconOutput$end_snp - falconOutput$st_snp + 1,
                    copyNo )
    colnames( output ) = c( 'Sample', 'Chromosome','Start','End','Num_Probes',  
                            'Segment_Mean' )  
    write.table( output, file = paste0( outPath, "/", sampleNameTemp, '.seg' ),
                 sep = '\t', quote = F,  row.names = F )
  }
  
}