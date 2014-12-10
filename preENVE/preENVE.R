rm(list=ls())
preENVE <- getwd()
preENVE <- '/Projects/ENVE-1.0-Beta/preENVE'
setwd(preENVE)

source(paste(preENVE,"Scripts","engine.R",sep='/'))
source(paste(preENVE,"Scripts","Settings.txt",sep='/'))
source(paste(preENVE,"Scripts","preENVE_PROJ_Config.txt",sep='/'))

dwnPack("plyr")

dir_create()
samp_proc(NormBam)
samp_proc(TumBam)
samp_info_proc(samp_info_file)

NorNor_dataRatio_calc()
NorNorScript_gen(NorNor_samp_DataRatio)

TumNor_dataRatio_calc()
TumNorScript_gen(TumNor_samp_DataRatio)
