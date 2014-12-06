rm(list=ls())
prENVE <- getwd()
setwd(prENVE)

source(paste(prENVE,"Scripts","engine.R",sep='/'))
source(paste(prENVE,"Scripts","prENVE_config.txt",sep='/'))
source(paste(prENVE,"Scripts","prENVE_PROJ_Config.txt",sep='/'))

dwnPack("plyr")

dir_create()
samp_proc(NormBam)
samp_proc(TumBam)
samp_info_proc(samp_info_file)

NorNor_dataRatio_calc()
NorNorScript_gen(NorNor_samp_DataRatio)

TumNor_dataRatio_calc()
TumNorScript_gen(TumNor_samp_DataRatio)
