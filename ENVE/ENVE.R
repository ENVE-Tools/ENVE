#!/usr/bin/env Rscript
rm(list=ls())
############################################################################
enveHome <<- "/Projects/ENVE-1.0-Beta/ENVE"
#enveHome <<- getwd()
scriptsPath <<-paste(enveHome,"scripts",sep= "/")
supFiles <<- paste(enveHome,"support_files",sep="/")
############################################################################
source(paste(scriptsPath,"engine.R",sep='/'))
source(paste(scriptsPath,"Settings.txt",sep='/'))
source(paste(scriptsPath,"ENVE_RUN_CONF.txt",sep='/'))
############################################################################


####Check_download_install_check all required R/Bioconductor Packages#######
dwnPack("DNAcopy")
dwnPack("stringr")
dwnPack("permute")
dwnPack("fExtremes")
dwnPack("IRanges")
dwnPack("ggplot2")
dwnPack("grid")
############################################################################
#########Create all required directories####################################
dir_create()
############################################################################
chr_lengths <- chr_proc()
############################################################################
segval_windows = NULL
segval_windows = seq(0, 2, 0.05)
min_seg = 50
pval_sig = 0.05
num_probes= 50


############################################################################
if(NormNorm)
{    
    reqd_files = reqd_files_func(Input_NormNorm_Files_Info)
    mode_correction(anaTempVScanNNres,anaTempVScanNN_MODE_CRCT)
    GISTIC <-FALSE
    filtCDS(anaTempVScanNN_MODE_CRCT,anaTempVScanNN_OC_CDSFilt)
    CBS_seg_samp(anaTempVScanNN_OC_CDSFilt,anaTempVScanNN_CBS_GC_crtd)
    Com_samp_perchr(anaTempVScanNN_CBS_GC_crtd,anaTempVScanNN_NorNor_SegMeans_CDSFilt)
    chr_pos_neg_sep(anaTempVScanNN_NorNor_SegMeans_CDSFilt)
    Nor_EVD_calc(anaTempVScanNN_NorNor_Pos,anaTempVScanNN_Pos_EVD_Cutoff,anaTempVScanNN_Pos_Tiff_output)
    Nor_EVD_calc(anaTempVScanNN_NorNor_Neg,anaTempVScanNN_Neg_EVD_Cutoff,anaTempVScanNN_Neg_Tiff_output)
}  
############################################################################
if(TumNorm)
{
  reqd_files = reqd_files_func(Input_TumNorm_Files_Info)
  mode_correction(anaTempVScanTNres,anaTempVScanTN_MODE_CRCT)
  cn_called_files = called_files(Input_TumNorm_Files_Info)
  GISTIC <- TRUE
  filtCDS(anaTempVScanTN_MODE_CRCT,anaTempVScanTN_OC_CDSFilt)
  CBS_seg_samp(anaTempVScanTN_OC_CDSFilt,anaTempVScanTN_CBS_GC_crtd)
  Com_samp_perchr(anaTempVScanTN_CBS_GC_crtd,anaTempVScanTN_TumNor_SegMeans_CDSFilt)
  TumEVD_cal()
}
############################################################################



