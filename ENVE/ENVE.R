#!/usr/bin/env Rscript
rm(list=ls())
############################################################################
#enveHome <<- "/projects/ENVE"
enveHome <<- getwd()
scriptsPath <<-paste(enveHome,"scripts",sep= "/")
supFiles <<- paste(enveHome,"support_files",sep="/")




source(paste(scriptsPath,"engine.R",sep='/'))
source(paste(scriptsPath,"ENVE_CONF.txt",sep='/'))
source(paste(scriptsPath,"ENVE_RUN_CONF.txt",sep='/'))
############################################################################


####Check_download_install_check all required R/Bioconductor Packages#######
dwnPack("DNAcopy")
dwnPack("stringr")
dwnPack("permute")
dwnPack("fExtremes")
dwnPack("IRanges")
############################################################################
#########Create all required directories####################################
dir_create()
############################################################################
chr_lengths <- chr_proc()
############################################################################

<<<<<<< Updated upstream


=======
#anaTempVScanNNres <- '/Volumes/Salendra_Data/Analysis54/temp/NorNor/copycaller_res'
setwd(anaTempVScanNNres)
>>>>>>> Stashed changes






############################################################################
if(NormNorm)
{   
<<<<<<< Updated upstream
    anaTempVScanNNres <- Input_NormNorm_adjlogratio_files
=======
    anaTempVScanNNres <- Input_NormNorm_GC_COR_LOGRATIO_files
>>>>>>> Stashed changes
    cn_called_files = filtCDS(anaTempVScanNNres,anaTempVScanNN_OC_CDSFilt)
    CBS_seg_samp(anaTempVScanNN_OC_CDSFilt,anaTempVScanNN_CBS_GC_crtd)
    Com_samp_perchr(anaTempVScanNN_CBS_GC_crtd,anaTempVScanNN_NorNor_SegMeans_CDSFilt)
    Nor_EVD_calc()
}  
############################################################################
if(TumNorm)
{
<<<<<<< Updated upstream
  anaTempVScanNNres <- Input_TumNorm_adjlogratio_files
=======
  anaTempVScanTNres <- Input_TumNorm_GC_COR_LOGRATIO_files 
>>>>>>> Stashed changes
  cn_called_files = filtCDS(anaTempVScanTNres,anaTempVScanTN_OC_CDSFilt)
  CBS_seg_samp(anaTempVScanTN_OC_CDSFilt,anaTempVScanTN_CBS_GC_crtd)
  Com_samp_perchr(anaTempVScanTN_CBS_GC_crtd,anaTempVScanTN_TumNor_SegMeans_CDSFilt)
  TumEVD_cal()
}
############################################################################



