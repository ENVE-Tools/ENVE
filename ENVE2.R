rm(list=ls())
############################################################################
enveHome <<- "/projects/ENVE"

source("/Projects/ENVE/scripts/engine2.R")
source("/Projects/ENVE/scripts/ENVE_config.R")

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
segval_windows = NULL
segval_windows = seq(0, 1, 0.05)
min_seg = 50
#anaTempVScanNNres <- '/Volumes/Salendra_Data/Analysis54/temp/NorNor/copycaller_res'
#setwd(anaTempVScanNNres)






############################################################################
if(NorNor)
{    
    cn_called_files = filtCDS(anaTempVScanNNres,anaTempVScanNN_OC_CDSFilt)
    CBS_seg_samp(anaTempVScanNN_OC_CDSFilt,anaTempVScanNN_CBS_GC_crtd)
    Com_samp_perchr(anaTempVScanNN_CBS_GC_crtd,anaTempVScanNN_NorNor_SegMeans_CDSFilt)
    Nor_EVD_calc()
}  
############################################################################
if(TumNor)
{
  cn_called_files = filtCDS(anaTempVScanTNres,anaTempVScanTN_OC_CDSFilt)
  CBS_seg_samp(anaTempVScanTN_OC_CDSFilt,anaTempVScanTN_CBS_GC_crtd)
  Com_samp_perchr(anaTempVScanTN_CBS_GC_crtd,anaTempVScanTN_TumNor_SegMeans_CDSFilt)
  TumEVD_cal()
}
############################################################################



