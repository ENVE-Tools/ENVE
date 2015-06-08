#!/usr/bin/env Rscript
rm(list=ls())
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
      preENVE <<- getwd()
      writeLines(" Usage : Rscript preENVE.R [-R] [ENVE_RUN_CONF(follow the template for creating CONF FILE : /ENVE/scripts/ENVE_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
      system(paste("cat ",preENVE,"/README.txt",sep=""))
      return(FALSE)
    }else if(args[1]=="-R")
    {
      if(len_args!=2)
      {
        writeLines(" Usage : Rscript preENVE.R [-R] [ENVE_RUN_CONF(follow the template for creating CONF FILE : /ENVE/scripts/ENVE_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
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
      writeLines(" Usage : Rscript preENVE.R [-R] [ENVE_RUN_CONF(follow the template for creating CONF FILE : /ENVE/scripts/ENVE_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
      return(FALSE)
    }
  }
}


run <- run_chk()





if(run)
{
source(paste(preENVE,"Scripts","engine.R",sep='/'))
dwnPack("plyr")
dir_create()
samp_proc(NormBam)
samp_proc(TumBam)
samp_info_proc(samp_info_file)

NorNor_dataRatio_calc()
NorNorScript_gen(NorNor_samp_DataRatio)

TumNor_dataRatio_calc()
TumNorScript_gen(TumNor_samp_DataRatio)
}