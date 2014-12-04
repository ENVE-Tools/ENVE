rm(list=ls())
prENVE <- getwd()
prENVE <- '/Projects/prENVE/'
setwd(prENVE)



source(paste(prENVE,"Scripts","engine2.R",sep='/'))
source(paste(prENVE,"Scripts","prENVE_config.R",sep='/'))
source(paste(prENVE,"Scripts","prENVE_PROJ_Config.R",sep='/'))
#source("/Projects/prENVE/Scripts/engine.R")
#source("/Projects/prENVE/Scripts/prENVE_config.R")
#source("/Projects/prENVE/Scripts/prENVE_PROJ_Config.R")


####Check_download_install_check all required R/Bioconductor Packages#######
dwnPack("DNAcopy")
dwnPack("stringr")
dwnPack("permute")
dwnPack("fExtremes")
dwnPack("IRanges")
dwnPack("permute")
dwnPack("combinat")
dwnPack("gtools")
dwnPack("plyr")
dwnPack("compare")


dir_create()
samp_proc(NormBam)
samp_proc(TumBam)
samp_info_proc(samp_info_file)
NorNor_dataRatio_calc()
NorNorScript_gen(NorNor_samp_DataRatio)

TumNor_dataRatio_calc()
TumNorScript_gen(TumNor_samp_DataRatio)
