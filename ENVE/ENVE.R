#!/usr/bin/env Rscript
rm(list=ls())
#########################################################################################################
run_chk <- function()
{
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  
  #args <- c("-R","/Projects/ENVE/scripts/ENVE_CONF.txt")
  
  len_args <- length(args)
  if(len_args == 0)
  {
    writeLines(" Usage : Rscript ENVE.R [-R] [ENVE_RUN_CONF(follow the template for creating CONF FILE : /ENVE/scripts/ENVE_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
  }else{
    if(args[1]=="-H")
    {
      enveHome <<- getwd()
      read_me <<- unlist(strsplit(enveHome,split='/'))
      read_me <<- paste(read_me[-length(read_me)],collapse='/')
      writeLines(" Usage : Rscript ENVE.R [-R] [ENVE_RUN_CONF(follow the template for creating CONF FILE : /ENVE/scripts/ENVE_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
      system(paste("cat ",read_me,"/README.txt",sep=""))
      return(FALSE)
    }else if(args[1]=="-R")
    {
      if(len_args!=2)
      {
        writeLines(" Usage : Rscript ENVE.R [-R] [ENVE_RUN_CONF(follow the template for creating CONF FILE : /ENVE/scripts/ENVE_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
        return(FALSE)
      } else
      {
        if(file.exists(args[2]))
        {
          source(args[2])
          return(TRUE)
        }else{
          print("ENVE RUN CONFIG FILE MISSING")
        }
      }
    }else{
      writeLines(" Usage : Rscript ENVE.R [-R] [ENVE_RUN_CONF(follow the template for creating CONF FILE : /ENVE/scripts/ENVE_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
      return(FALSE)
    }
  }
}
################################################################################################################################
run <- run_chk()
if(run)
{
  scriptsPath <<-paste(enveHome,"scripts",sep= "/")
  supFiles <<- paste(enveHome,"support_files",sep="/")
  source(paste(scriptsPath,"engine.R",sep='/'))  
  ####Check_download_install_check all required R/Bioconductor Packages#######
  dwnPack("DNAcopy")
  dwnPack("stringr")
  dwnPack("permute")
  dwnPack("fExtremes")
  dwnPack("IRanges")
  dwnPack("ggplot2")
  dwnPack("grid")
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
    mode_correction(anaTempLRNNres,anaTempLRNN_MODE_CRCT)
    GISTIC <-FALSE
    filtCDS(anaTempLRNN_MODE_CRCT,anaTempLRNN_OC_CDSFilt)
    CBS_seg_samp(anaTempLRNN_OC_CDSFilt,anaTempLRNN_CBS_GC_crtd)
    Com_samp_perchr(anaTempLRNN_CBS_GC_crtd,anaTempLRNN_NorNor_SegMeans_CDSFilt)
    chr_pos_neg_sep(anaTempLRNN_NorNor_SegMeans_CDSFilt)
    Nor_EVD_calc(anaTempLRNN_NorNor_Pos,anaTempLRNN_Pos_EVD_Cutoff,anaTempLRNN_Pos_Tiff_output)
    Nor_EVD_calc(anaTempLRNN_NorNor_Neg,anaTempLRNN_Neg_EVD_Cutoff,anaTempLRNN_Neg_Tiff_output)
  }  
  ############################################################################
  if(TumNorm)
  {
    reqd_files = reqd_files_func(Input_TumNorm_Files_Info)
    mode_correction(anaTempLRTNres,anaTempLRTN_MODE_CRCT)
    cn_called_files = called_files(Input_TumNorm_Files_Info)
    GISTIC <- TRUE
    filtCDS(anaTempLRTN_MODE_CRCT,anaTempLRTN_OC_CDSFilt)
    CBS_seg_samp(anaTempLRTN_OC_CDSFilt,anaTempLRTN_CBS_GC_crtd)
    Com_samp_perchr(anaTempLRTN_CBS_GC_crtd,anaTempLRTN_TumNor_SegMeans_CDSFilt)
    TumEVD_cal()
  }
  ############################################################################
  
  
}



