##########################################################################################################################################################################################
dwnPack <- function(x)
{
  if(!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE))
    {
      source("http://bioconductor.org/biocLite.R")
      biocLite(x)
      if(is.element(x, installed.packages()[,1])) stop("Package not found")
    }
  }
}
##########################################################################################################################################################################################
dir_create <- function()
{
  print(paste(toString(Sys.time()),"Creating Directories for Analysis"))
  ana <<- paste(prENVE,"Analysis",sep='/')
  #####Analysis files#########
  
  anaPath <<- paste(ana,(paste("Analysis",toString(strftime(Sys.time(),format="%H_%M_%d_%m_%Y")),sep='_')),sep ="/")
  anaInp <<- paste(anaPath,"Input",sep='/')
  anaTemp <<-paste(anaPath,"temp",sep='/')
  anaTemp_NorNor <<-paste(anaTemp,'NorNor',sep='/')
  anaTemp_NorNor_mpileup_res <<- paste(anaTemp_NorNor,'mpileup_res',sep='/')
  anaTemp_NorNor_copynumer_res <<- paste(anaTemp_NorNor,'copynumber_res',sep='/')
  anaTemp_NorNor_copycaller_res <<- paste(anaTemp_NorNor,'copycaller_res',sep='/')
  anaTemp_NorNor_ResScripts <<- paste(anaTemp_NorNor,'scripts',sep='/')
  anaTemp_NorNor_som_res <<- paste(anaTemp_NorNor,'somatic_res',sep='/')
  anaTemp_NorNor_adj_logratio <<-paste(anaTemp_NorNor,'adjusted_logratio',sep='/')
  
  anaTemp_TumNor <<- paste(anaTemp,'TumNor',sep='/')
  anaTemp_TumNor_mpileup_res <<- paste(anaTemp_TumNor,'mpileup_res',sep='/')
  anaTemp_TumNor_copynumer_res <<- paste(anaTemp_TumNor,'copynumber_res',sep='/')
  anaTemp_TumNor_copycaller_res <<- paste(anaTemp_TumNor,'copycaller_res',sep='/')
  anaTemp_TumNor_ResScripts <<- paste(anaTemp_TumNor,'scripts',sep='/')
  anaTemp_TumNor_som_res <<- paste(anaTemp_TumNor,'somatic_res',sep='/')
  anaTemp_TumNor_adj_logratio <<-paste(anaTemp_TumNor,'adjusted_logratio',sep='/')
  
  ################## Insert more dirs######################  
  dirs <- c(ana,
            anaPath,
            anaInp,
            anaTemp,
            anaTemp_NorNor,
            anaTemp_NorNor_mpileup_res,
            anaTemp_NorNor_som_res,
            anaTemp_NorNor_copynumer_res,
            anaTemp_NorNor_copycaller_res,
            anaTemp_NorNor_ResScripts,
            anaTemp_NorNor_adj_logratio,
            anaTemp_TumNor,
            anaTemp_TumNor_mpileup_res,
            anaTemp_NorNor_som_res,
            anaTemp_TumNor_copynumer_res,
            anaTemp_TumNor_copycaller_res,
            anaTemp_TumNor_ResScripts,
            anaTemp_TumNor_adj_logratio
          )
  dir2 <- as.data.frame(cbind(c('ana',
                                'anaPath',
                                'anaInp',
                                'anaTemp',
                                'anaTemp_NorNor',
                                'anaTemp_NorNor_mpileup_res',
                                'anaTemp_NorNor_som_res',
                                'anaTemp_NorNor_copynumer_res',
                                'anaTemp_NorNor_copycaller_res',
                                'anaTemp_NorNor_ResScripts',
                                'anaTemp_TumNor_adj_logratio',
                                'anaTemp_TumNor',
                                'anaTemp_TumNor_mpileup_res',
                                'anaTemp_NorNor_som_res',
                                'anaTemp_TumNor_copynumer_res',
                                'anaTemp_TumNor_copycaller_res',
                                'anaTemp_TumNor_ResScripts',
                                'anaTemp_TumNor_adj_logratio'), dirs))
  colnames(dir2) <- c('dirs','location')
  
  for(i in 1:length(dirs))
  {
    if (file.exists(dirs[i])){
    } else {
      dir.create(dirs[i])
    }
  }
  
  write.table(dir2, file=paste(anaTemp, "Directories.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  return(dir2)
  print(paste(toString(Sys.time()),"Required Directories Created"))
}

##########################################################################################################################################################################################

samp_proc <- function(x)
{
  print(paste(toString(Sys.time()),
              paste("Processing the",
                    deparse(substitute(x)),
                    sep= ' '),
              "files",
              sep=' ')
  )
  sys_cmnd2 <- paste('ls',x,'>',paste(anaInp,paste(deparse(substitute(x)),'txt',sep='.'),sep='/'),sep =" ")
  system(sys_cmnd2)
}

##########################################################################################################################################################################################

samp_info_proc <- function(z)
{
  half_number = floor(Number_of_samp/2)
  ##Import the Sample sheet and convert all factors to characters##
  smp_info <-as.data.frame(read.delim(file = z,header=T,sep ='\t' ),row.names=NULL,optional = F, stringsAsFactors = F)
  i <- sapply(smp_info, is.factor)
  smp_info[i] <- lapply(smp_info[i], as.character)
  ###########
  
  ##unifying gender in Male and Female##
  for(i in 1:length(smp_info[,1]))
  {
    if(smp_info[i,6] =="M" |smp_info[i,6] =="Male" |smp_info[i,6] =="male" | smp_info[i,6] =="m")
    {
      #print(smp_info[i,3])
      smp_info[i,6] = gsub(smp_info[i,6],"Male",smp_info[i,6])
    }
    if(smp_info[i,6] =="F" |smp_info[i,6] =="Female" |smp_info[i,6] =="female" | smp_info[i,6] =="f")
    {
      #print(smp_info[i,6])
      smp_info[i,6] = gsub(smp_info[i,6],"Female",smp_info[i,6])
    }
  }
  ###########
  
  ## chking for duplicate patient ids ##
  dup_pat_ids <- count(duplicated(smp_info[,1], incomparables = F, fromLast = F) | is.na(smp_info[,1]))[2,2]
  
  ###########
  
  ## If there are no duplicate patient Ids ## 
  
  if(dup_pat_ids ==0 | is.na(dup_pat_ids))
  {
    GNS = 0                                                         # Good Normal Samples
    BNS = 0                                                         # Bad Normal Samples
    GTS = 0                                                         # Good Tumor Samples
    BTS = 0                                                         # Bad Tumor Samples
    ## For all Normal samples Ids chk if UQ base aligned reads and Gender is not missing, if yes return error, or continue ##
    for(j in 1:length(smp_info[,1]))
    {
      if(!is.na(smp_info[j,2]))
      {
        if(is.na(smp_info[j,3]) | is.na(smp_info[j,6]))
        {
          print(paste("data missing for",smp_info[j,2],"in samplesheet",sep=' '))
          BNS = BNS + 1
        }else{
          GNS = GNS +1
        }
      }
    }
    print(GNS)
    print(BNS)
    ## If the data for Normal samples is OK, creating Normal samplesheet for calculating the dataratio
    if(BNS==0)
    {
      NorSampInfo <- smp_info[complete.cases(smp_info[,c(2,3,6)]),c(2,3,6)]
      NorSampInfo2 <- unique(NorSampInfo[,1:3])
      write.table(NorSampInfo2, file=paste(anaInp, "nor_samp_info.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t") # Returning to DataRatio calculation
      NorSamps_Avl <- as.data.frame(read.delim(file = paste(anaInp,'NormBam.txt',sep='/'),header=F,sep ='\t' ),row.names=NULL,optional = F, stringsAsFactors = F, na.strings = 'NA')
      i <- sapply(NorSamps_Avl, is.factor)
      NorSamps_Avl[i] <- lapply(NorSamps_Avl[i], as.character)
      diff1  <- data.frame(setdiff(NorSampInfo2[,1],NorSamps_Avl[,1]))
      i <- sapply(diff1, is.factor)
      diff1[i] <- lapply(diff1[i], as.character)
      if(length(diff1[,1])!=0)
      {
        for(k in 1:length(diff1[,1]))
        {
          print(paste("No Files for Normal Samples in the directory for",diff1[k,1], sep=' '))
        }
      }else{
        male_nor_samps <- NorSampInfo2[which(NorSampInfo2[,3]=='Male'),1]
        female_nor_samps <- NorSampInfo2[which(NorSampInfo2[,3]=='Female'),1]
        if(length(NorSampInfo2[,1] <= Number_of_samp))
        {
          if(length(male_nor_samps) >= length(female_nor_samps))
          {
            if(length(female_nor_samps) >= half_number)
            {
              nor_list_f <- sample(female_nor_samps,half_number)
              nor_list_m <- sample(male_nor_samps,half_number)
              #nor_list <- cbind(nor_list_f,nor_list_m)
            }else{
              nor_list_f <- sample(female_nor_samps,length(female_nor_samps))
              nor_list_m <- sample(male_nor_samps,(Number_of_samp-length(female_nor_samps)))
              #nor_list <- cbind(nor_list_f,nor_list_m)
            }
          }else{
            if(length(male_nor_samps) >= half_number)
            {
              nor_list_f <- sample(female_nor_samps,half_number)
              nor_list_m <- sample(male_nor_samps,half_number)
              #nor_list <- cbind(nor_list_f,nor_list_m)
            }else{
              nor_list_m <- sample(male_nor_samps,length(male_nor_samps))
              nor_list_f <- sample(female_nor_samps,(Number_of_samp-length(male_nor_samps)))
              #nor_list <- cbind(nor_list_f,nor_list_m)
            }
          }
          nor_list_f <- as.data.frame(nor_list_f,col.names=F)
          colnames(nor_list_f) <- 'sample'
          i <- sapply(nor_list_f, is.factor)
          nor_list_f[i] <- lapply(nor_list_f[i], as.character)
          nor_list_m <- as.data.frame(nor_list_m,col.names=F)
          colnames(nor_list_m) <- 'sample'
          i <- sapply(nor_list_m, is.factor)
          nor_list_m[i] <- lapply(nor_list_m[i], as.character)
          nor_list <- rbind(nor_list_m,nor_list_f)
          write.table(nor_list, file=paste(anaInp, "Norm_Samp_List.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
         
        }else
        {
          print("number of samples asked is less than number of samples present")
        }
        
      }
      if(BNS == 0)
      {
        for(j in 1:length(smp_info[,1]))
        {
          if(!is.na(smp_info[j,4]))
          {
            if(is.na(smp_info[j,5]) | is.na(smp_info[j,6]))
            {
              print(paste("data missing for",smp_info[j,4],sep=' '))
              BTS = BTS + 1
            }else{
              if(is.na(smp_info[j,2]))
              {
                print(paste("Matched Normal Sample Missing for",smp_info[j,4],"for patient id",smp_info[j,1],"is missing",sep=' '))
                BTS = BTS + 1
              }
              GTS = GTS +1
            }
          }
        }
        
        
      }
      if(BTS==0)
      {
        TumSampInfo <- smp_info[!is.na(smp_info[,4]),4]
        TumSampInfo2 <- unique(TumSampInfo)
        TumSamps_Avl <- as.data.frame(read.delim(file = paste(anaInp,'TumBam.txt',sep='/'),header=F,sep ='\t' ),row.names=NULL,optional = F, stringsAsFactors = F, na.strings = 'NA')
        i <- sapply(TumSamps_Avl, is.factor)
        TumSamps_Avl[i] <- lapply(TumSamps_Avl[i], as.character)
        diff2  <- data.frame(setdiff(TumSampInfo2,TumSamps_Avl[,1]))
        i <- sapply(diff2, is.factor)
        diff2[i] <- lapply(diff2[i], as.character)
        if(length(diff2[,1])!=0)
        {
          for(k in 1:length(diff2[,1]))
          {
            print(paste("No Files for tumor in the directory for sample",diff2[k,1], sep=' '))
          }
        }else{
          TumSampInfo <- smp_info[complete.cases(smp_info[,c(4,5,6)]),c(2:6)]
          TumSampInfo2 <- unique(TumSampInfo[,1:5])
          write.table(TumSampInfo2, file=paste(anaInp, "tum_samp_info.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t")
          
        }
        
      }
      print(GTS)
      print(BTS)
    }
  }else{
    print(paste("There are",dup_pat_ids,"duplicate patient ids, Kindly remove them before continuiing",sep=' '))
  }
}


NorNor_dataRatio_calc <-function()
{
  
  smp_info <-as.data.frame(read.delim(file = paste(anaInp, "nor_samp_info.txt", sep="/"),header=T,sep ='\t' ),row.names=NULL,optional = F, stringsAsFactors = F)
  i <- sapply(smp_info, is.factor)
  smp_info[i] <- lapply(smp_info[i], as.character)
  colnames(smp_info) <- c("Sample","PF_UQ_BASES_ALIGNED","Gender")
  print(paste(toString(Sys.time()),"Creating Data Ratio files for Normal-Normal pairs"))
  
  
  x_list <- as.data.frame(read.delim(file = paste(anaInp,'Norm_Samp_List.txt',sep='/'),header=F,sep ='\t' ),row.names=NULL,optional = F)
  i <- sapply(x_list, is.factor)
  x_list[i] <- lapply(x_list[i], as.character)
  x_mat <- t(combn(x_list[,1],2))
  colnames(x_mat) <- c('Sample','Sample2')
  x_mat <- merge(x_mat,smp_info, by.x_mat = 'Sample', by.smp_info='Sample',all=F)
  colnames(x_mat)<- c("Sample1","Sample","Uq_Bas_Ali_samp1","Gender_Samp1")
  x_mat <- merge(x_mat,smp_info, by.x_mat = 'Sample', by.smp_info='Sample',all=F)
  colnames(x_mat)<- c("Sample2","Sample1","Uq_Bas_Ali_samp1","Gender_Samp1","Uq_Bas_Ali_samp2","Gender_Samp2")
  x_mat <- cbind(x_mat,round(x_mat[,3]/x_mat[,5], digits=10))
  colnames(x_mat)<- c("Sample2","Sample1","Uq_Bas_Ali_samp1","Gender_Samp1","Uq_Bas_Ali_samp2","Gender_Samp2","Data_Ratio")
  x_mat <- x_mat[,c("Sample1","Sample2","Uq_Bas_Ali_samp1","Gender_Samp1","Uq_Bas_Ali_samp2","Gender_Samp2","Data_Ratio")]
  if(DEBUG)
  {
    View(x_mat)
  }
  write.table(x_mat,file=paste(anaTemp_NorNor, "NorNor_samp_DataRatio.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t")
  print(paste(toString(Sys.time()),"Data Ratio files for Normal-Normal Pairs Created and Saved to",anaTemp_NorNor,sep=" "))
}

##########################################################################################################################################################################################


NorNorScript_gen <- function(x)
{
  NorNor_Dr <- as.data.frame(read.delim(file = paste(anaTemp_NorNor,paste(deparse(substitute(x)),'txt',sep='.'),sep='/'),header=T,sep ='\t' ),row.names=NULL,optional = F)
  NorNor_mpileup <- cbind(NorNor_Dr,
                          paste(paste(samTools,'samtools',sep='/'),
                                'mpileup -f',
                                hg19_karyo,
                                paste(NormBam,NorNor_Dr[,1],sep='/'),
                                paste(NormBam,NorNor_Dr[,2],sep='/'),
                                '>',
                                paste(anaTemp_NorNor_mpileup_res,
                                      paste(paste('mpileup_res',
                                                  sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                  "Normal",
                                                  sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                  "Normal",
                                                  paste(NorNor_Dr[,4],
                                                        NorNor_Dr[,6],
                                                        sep=''),
                                                  sep='_'),
                                            'mpileup',
                                            sep='.'),
                                      sep='/'),
                                sep=' ')
  )
  colnames(NorNor_mpileup) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','mpileup_command') 
  write.table(NorNor_mpileup[,8], file=paste(anaTemp_NorNor_ResScripts, "NorNor_mpileup_commands.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t")
  if(BASH)
  {
    NorNor_mpileup_bash <- as.data.frame(NorNor_mpileup[,8],stringsAsFactors = F)
    NorNor_mpileup_bash <- sapply(NorNor_mpileup_bash, as.character)
    NorNor_mpileup_bash2 <- rbind(as.character('#!/bin/bash'), NorNor_mpileup_bash)
    #View(NorNor_mpileup_bash2)
    write.table(NorNor_mpileup_bash2, file=paste(anaTemp_NorNor_ResScripts, "NorNor_mpileup_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  if(Com_Scr)
  {
    mpileup_cmnds <<- as.data.frame(NorNor_mpileup[,8],stringsAsFactors = F, row.names= NULL, optional = F, colnames = F)
    mpileup_cmnds[,1] <<- sapply(mpileup_cmnds[,1], as.character)
  }
  print(paste(toString(Sys.time()),"samtools mpileup generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_NorNor_ResScripts,sep=" "))
  if(DEBUG)
  {
    View(NorNor_mpileup)
  }
  
  NorNor_copynumber <- cbind(NorNor_Dr,
                             paste(paste(JAVA_HOME,'java',sep='/'),
                                   '-jar',
                                   VarScan,
                                   'copynumber',
                                   paste(anaTemp_NorNor_mpileup_res, 
                                         paste(paste('mpileup_res',
                                                     sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                     "Normal",
                                                     paste(NorNor_Dr[,4],
                                                           NorNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'mpileup',
                                               sep='.'),
                                         sep='/'),
                                   paste(anaTemp_NorNor_copynumer_res,
                                         paste(paste('varscanCN',
                                                     sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                     "Normal",
                                                     paste(NorNor_Dr[,4],
                                                           NorNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               sep='.'),
                                         sep='/'),
                                   '--mpileup',
                                   '1',
                                   '--data-ratio',
                                   NorNor_Dr[,7],
                                   sep=' ')
  )
  
  colnames(NorNor_copynumber) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_copynumber_commands')
  write.table(NorNor_copynumber[,8], file=paste(anaTemp_NorNor_ResScripts, "NorNor_copynumber_commands.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t") 
  if(BASH)
  {
    NorNor_copynumber_bash <- as.data.frame(NorNor_copynumber[,8],stringsAsFactors = F)
    NorNor_copynumber_bash <- sapply(NorNor_copynumber_bash, as.character)
    NorNor_copynumber_bash2 <- rbind(as.character('#!/bin/bash'), NorNor_copynumber_bash)
    #View(NorNor_copynumber_bash2)
    write.table(NorNor_copynumber_bash2, file=paste(anaTemp_NorNor_ResScripts, "NorNor_copynumber_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  if(Com_Scr)
  {
    copynumber_cmnds <<- as.data.frame(NorNor_copynumber[,8],stringsAsFactors = F, row.names= NULL, optional = F, colnames = F)
    copynumber_cmnds[,1] <<- sapply(copynumber_cmnds[,1], as.character)
  }
  if(DEBUG)
  {
    View(NorNor_copynumber)
  }
  
  NorNor_somatic <- cbind(NorNor_Dr,
                          paste(paste(JAVA_HOME,'java',sep='/'),
                                '-jar',
                                VarScan,
                                'somatic',
                                paste(anaTemp_NorNor_mpileup_res, 
                                      paste(paste('mpileup_res',
                                                  sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                  "Normal",
                                                  sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                  "Normal",
                                                  paste(NorNor_Dr[,4],
                                                        NorNor_Dr[,6],
                                                        sep=''),
                                                  sep='_'),
                                            'mpileup',
                                            sep='.'),
                                      sep='/'),
                                paste(anaTemp_NorNor_som_res,
                                      paste(paste('varscanSOM',
                                                  sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                  "Normal",
                                                  sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                  "Normal",
                                                  paste(NorNor_Dr[,4],
                                                        NorNor_Dr[,6],
                                                        sep=''),
                                                  sep='_'),
                                            'out',
                                            sep='.'),
                                      sep='/'),
                                '--mpileup',
                                '1',
                                sep=' ')
  )
  colnames(NorNor_copynumber) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_somatic_commands')
  write.table(NorNor_somatic[,8], file=paste(anaTemp_NorNor_ResScripts, "NorNor_somatic_commands.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t")
  if(BASH)
  {
    NorNor_somatic_bash <- as.data.frame(NorNor_somatic[,8],stringsAsFactors = F)
    NorNor_somatic_bash <- sapply(NorNor_somatic_bash, as.character)
    NorNor_somatic_bash2 <- rbind(as.character('#!/bin/bash'), NorNor_somatic_bash)
    #View(NorNor_copynumber_bash2)
    write.table(NorNor_somatic_bash2, file=paste(anaTemp_NorNor_ResScripts, "NorNor_somatic_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  
  if(F)
  {
    somatic_cmnds <<- as.data.frame(NorNor_somatic[,8],stringsAsFactors = F, row.names= NULL, optional = F, colnames = F)
    somatic_cmnds[,1] <<- sapply(somatic_cmnds[,1], as.character)
  }
  if(DEBUG)
  {
    View(NorNor_somatic)
  }
  
  
  print(paste(toString(Sys.time()),"Varscan_copynumber generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_NorNor_ResScripts,sep=" "))
  if(DEBUG)
  {
    View(NorNor_copynumber)
  }  
  
  NorNor_copycaller <- cbind(NorNor_Dr,
                             paste(paste(JAVA_HOME,'java',sep='/'),
                                   '-jar',
                                   VarScan,
                                   'copyCaller',
                                   paste(anaTemp_NorNor_copynumer_res,
                                         paste(paste('varscanCN',
                                                     sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                     "Normal",
                                                     paste(NorNor_Dr[,4],
                                                           NorNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               'copynumber',
                                               sep='.'),
                                         sep='/'),
                                   '--output-file',
                                   paste(anaTemp_NorNor_copycaller_res,
                                         paste(paste('varscanCC',
                                                     sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                     "Normal",
                                                     paste(NorNor_Dr[,4],
                                                           NorNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               'called',
                                               sep='.'),
                                         sep='/'),
                                   '--output-homdel-file',
                                   paste(anaTemp_NorNor_copycaller_res,
                                         paste(paste('varscanCC',
                                                     sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                     "Normal",
                                                     paste(NorNor_Dr[,4],
                                                           NorNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               'called',
                                               'homdel',
                                               sep='.'),
                                         sep='/'),
                                   sep=' ')
  )
  colnames(NorNor_copycaller) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_CopyCaller_commands')
  write.table(NorNor_copycaller[,8], file=paste(anaTemp_NorNor_ResScripts, "NorNor_copycaller_commands.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  if(BASH)
  {
    NorNor_copycaller_bash <- as.data.frame(NorNor_copycaller[,8],stringsAsFactors = F)
    NorNor_copycaller_bash <- sapply(NorNor_copycaller_bash, as.character)
    NorNor_copycaller_bash2 <- rbind(as.character('#!/bin/bash'), NorNor_copycaller_bash)
    #View(NorNor_copynumber_bash2)
    write.table(NorNor_copycaller_bash2, file=paste(anaTemp_NorNor_ResScripts, "NorNor_copycaller_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  
  print(paste(toString(Sys.time()),"Varscan_copycaller generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_NorNor_ResScripts,sep=" "))
  if(Com_Scr)
  {
    copycaller_cmnds <<- as.data.frame(NorNor_copycaller[,8],stringsAsFactors = F, row.names= NULL, optional = F, colnames = F)
    copycaller_cmnds[,1] <<- sapply(copycaller_cmnds[,1], as.character)
  }
  
  if(DEBUG)
  {
    View(NorNor_copycaller)
  }
  
  NorNor_rm_mpileup <- cbind(NorNor_Dr,
                             paste('rm',
                                   paste(anaTemp_NorNor_mpileup_res,
                                         paste(paste('mpileup_res',
                                                     sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                     "Normal",
                                                     paste(NorNor_Dr[,4],
                                                           NorNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'mpileup',
                                               sep='.'),
                                         sep='/'),
                                   sep=' ')
  )
  colnames(NorNor_rm_mpileup) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','mpileup_rm_command') 
  write.table(NorNor_rm_mpileup[,8], file=paste(anaTemp_NorNor_ResScripts, "NorNor_rm_mpileup_commands.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  if(BASH)
  {
    NorNor_rm_mpileup_bash <- as.data.frame(NorNor_rm_mpileup[,8],stringsAsFactors = F)
    NorNor_rm_mpileup_bash <- sapply(NorNor_rm_mpileup_bash, as.character)
    NorNor_rm_mpileup_bash2 <- rbind(as.character('#!/bin/bash'), NorNor_rm_mpileup_bash)
    #View(NorNor_copynumber_bash2)
    write.table(NorNor_rm_mpileup_bash2, file=paste(anaTemp_NorNor_ResScripts, "NorNor_rm_mpileup_bash.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  
  print(paste(toString(Sys.time()),"mpileup files removal script written at",anaTemp_NorNor_ResScripts,sep=" "))
  
  if(Com_Scr)
  {
    rem_mpileup_cmnds <<- as.data.frame(NorNor_rm_mpileup[,8],stringsAsFactors = F, row.names= NULL, optional = F, colnames = F)
    rem_mpileup_cmnds[,1] <<- sapply(rem_mpileup_cmnds[,1], as.character)
  }
  
  if(DEBUG)
  {
    View(NorNor_rm_mpileup)
  }
  
  #################################################################################
  NorNor_adj_logratio <- cbind(NorNor_Dr,
                               paste('cut -f1,2,3,7',
                                     paste(anaTemp_NorNor_copycaller_res,
                                           paste(paste('varscanCC',
                                                       sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                       "Normal",
                                                       sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                       "Normal",
                                                       paste(NorNor_Dr[,4],
                                                             NorNor_Dr[,6],
                                                             sep=''),
                                                       sep='_'),
                                                 'out',
                                                 'called',
                                                 sep='.'),
                                           sep ='/'
                                     ),
                                     '>',
                                     paste(anaTemp_TumNor_adj_logratio,
                                           paste(paste('varscanCC',
                                                       sub("*.cleaned.bam","",NorNor_Dr[,1]),
                                                       "Normal",
                                                       sub("*.cleaned.bam","",NorNor_Dr[,2]),
                                                       "Normal",
                                                       paste(NorNor_Dr[,4],
                                                             NorNor_Dr[,6],
                                                             sep=''),
                                                       sep='_'),
                                                 'adj.logratio',
                                                 sep='.'),
                                           sep='/'),
                                        sep=' ')
                               )
  colnames(NorNor_adj_logratio) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_adj_logratio_commands')
  write.table(NorNor_adj_logratio[,8], file=paste(anaTemp_NorNor_ResScripts, "NorNor_adj_logratio.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  if(BASH)
  {
    NorNor_adj_logratio_bash <- as.data.frame(NorNor_adj_logratio[,8],stringsAsFactors = F)
    NorNor_adj_logratio_bash <- sapply(NorNor_adj_logratio_bash, as.character)
    NorNor_adj_logratio_bash2 <- rbind(as.character('#!/bin/bash'), NorNor_adj_logratio_bash)
    #View(NorNor_copynumber_bash2)
    write.table(NorNor_adj_logratio_bash2, file=paste(anaTemp_NorNor_ResScripts, "NorNor_adj_logratio_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  
  print(paste(toString(Sys.time()),"Adjusted logratio files generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_NorNor_ResScripts,sep=" "))
  if(Com_Scr)
  {
    adj_logratio_cmnds <<- as.data.frame(NorNor_adj_logratio[,8],stringsAsFactors = F, row.names= NULL, optional = F, colnames = F)
    adj_logratio_cmnds[,1] <<- sapply(adj_logratio_cmnds[,1], as.character)
  }
  
  if(DEBUG)
  {
    View(NorNor_adj_logratio)
  }
  #################################################################################
  
  
}
##########################################################################################################################################################################################

TumNor_dataRatio_calc <-function()
{
  smp_info <-as.data.frame(read.delim(file = paste(anaInp, "tum_samp_info.txt", sep="/"),header=T,sep ='\t' ),row.names=NULL,optional = F, stringsAsFactors = F)
  i <- sapply(smp_info, is.factor)
  smp_info[i] <- lapply(smp_info[i], as.character)
  x_mat2 <- cbind(smp_info,smp_info[,5],round(smp_info[,2]/smp_info[,4], digits=10))
  colnames(x_mat2)<- c("NORMAL_SAMPLE","NORMAL_UQ_BASES_ALIGNED","TUMOR_SAMPLE","TUMOR_UQ_BASES_ALIGNED","GENDER1","GENDER2","Data_Ratio")
  x_mat2 <- x_mat2[,c("NORMAL_SAMPLE","TUMOR_SAMPLE","NORMAL_UQ_BASES_ALIGNED","GENDER1","TUMOR_UQ_BASES_ALIGNED","GENDER2", "Data_Ratio")]
  write.table(x_mat2, file=paste(anaTemp_TumNor, "TumNor_samp_DataRatio.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t")
  print(paste(toString(Sys.time()),"Data Ratio files for Normal-Tumor Pairs Created and Saved to",anaTemp_TumNor,sep=" "))
  if(DEBUG)
  {
    View(x_mat2)
  }
}

##########################################################################################################################################################################################

TumNorScript_gen <- function(x)
{
  TumNor_Dr <- as.data.frame(read.delim(file = paste(anaTemp_TumNor,paste(deparse(substitute(x)),'txt',sep='.'),sep='/'),header=T,sep ='\t' ),row.names=NULL,optional = F)
  TumNor_mpileup <- cbind(TumNor_Dr,
                          paste(paste(samTools,'samtools',sep='/'),
                                'mpileup -f',
                                hg19_karyo,
                                paste(NormBam,TumNor_Dr[,1],sep='/'),
                                paste(TumBam,TumNor_Dr[,2],sep='/'),
                                '>',
                                paste(anaTemp_TumNor_mpileup_res,
                                      paste(paste('mpileup_res',
                                                  sub("*.bam","",TumNor_Dr[,1]),
                                                  "Normal",
                                                  sub("*.bam","",TumNor_Dr[,2]),
                                                  "Tumor",
                                                  paste(TumNor_Dr[,4],
                                                        TumNor_Dr[,6],
                                                        sep=''),
                                                  sep='_'),
                                            'mpileup',
                                            sep='.'),
                                      sep='/'),
                                sep=' ')
  )
  colnames(TumNor_mpileup) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','mpileup_commands') 
  write.table(TumNor_mpileup[,8], file=paste(anaTemp_TumNor_ResScripts, "TumNor_mpileup_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  if(BASH)
  {
    TumNor_mpileup_bash <- as.data.frame(TumNor_mpileup[,8],stringsAsFactors = F)
    TumNor_mpileup_bash <- sapply(TumNor_mpileup_bash, as.character)
    TumNor_mpileup_bash2 <- rbind(as.character('#!/bin/bash'), TumNor_mpileup_bash)
    #View(TumNor_mpileup_bash2)
    write.table(TumNor_mpileup_bash2, file=paste(anaTemp_TumNor_ResScripts, "TumNor_mpileup_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  print(paste(toString(Sys.time()),"samtools mpileup generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_TumNor_ResScripts,sep=" "))
  if(DEBUG)
  {
    View(TumNor_mpileup)
  }
  
  TumNor_copynumber <- cbind(TumNor_Dr,
                             paste(paste(JAVA_HOME,'java',sep='/'),
                                   '-jar',
                                   VarScan,
                                   'copynumber',
                                   paste(anaTemp_TumNor_mpileup_res, 
                                         paste(paste('mpileup_res',
                                                     sub("*.bam","",TumNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.bam","",TumNor_Dr[,2]),
                                                     "Tumor",
                                                     paste(TumNor_Dr[,4],
                                                           TumNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'mpileup',
                                               sep='.'),
                                         sep='/'),
                                   paste(anaTemp_TumNor_copynumer_res,
                                         paste(paste('varscanCN',
                                                     sub("*.bam","",TumNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.bam","",TumNor_Dr[,2]),
                                                     "Tumor",
                                                     paste(TumNor_Dr[,4],
                                                           TumNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               sep='.'),
                                         sep='/'),
                                   '--mpileup',
                                   '1',
                                   '--data-ratio',
                                   TumNor_Dr[,7],
                                   sep=' ')
  )
  
  colnames(TumNor_copynumber) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_copynumber_command')
  write.table(TumNor_copynumber[,8], file=paste(anaTemp_TumNor_ResScripts, "TumNor_copynumber_commands.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  if(BASH)
  {
    TumNor_copynumber_bash <- as.data.frame(TumNor_copynumber[,8],stringsAsFactors = F)
    TumNor_copynumber_bash <- sapply(TumNor_copynumber_bash, as.character)
    TumNor_copynumber_bash2 <- rbind(as.character('#!/bin/bash'), TumNor_copynumber_bash)
    #View(TumNor_mpileup_bash2)
    write.table(TumNor_copynumber_bash2, file=paste(anaTemp_TumNor_ResScripts, "TumNor_copynumber_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  
  TumNor_somatic <- cbind(TumNor_Dr,
                          paste(paste(JAVA_HOME,'java',sep='/'),
                                '-jar',
                                VarScan,
                                'somatic',
                                paste(anaTemp_TumNor_mpileup_res, 
                                      paste(paste('mpileup_res',
                                                  sub("*.bam","",TumNor_Dr[,1]),
                                                  "Normal",
                                                  sub("*.bam","",TumNor_Dr[,2]),
                                                  "Tumor",
                                                  paste(TumNor_Dr[,4],
                                                        TumNor_Dr[,6],
                                                        sep=''),
                                                  sep='_'),
                                            'mpileup',
                                            sep='.'),
                                      sep='/'),
                                paste(anaTemp_TumNor_som_res,
                                      paste(paste('varscanSOM',
                                                  sub("*.bam","",TumNor_Dr[,1]),
                                                  "Normal",
                                                  sub("*.bam","",TumNor_Dr[,2]),
                                                  "Tumor",
                                                  paste(TumNor_Dr[,4],
                                                        TumNor_Dr[,6],
                                                        sep=''),
                                                  sep='_'),
                                            'out',
                                            sep='.'),
                                      sep='/'),
                                '--mpileup',
                                '1',
                                sep=' ')
  )
  colnames(TumNor_copynumber) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_somatic_commands')
  write.table(TumNor_somatic[,8], file=paste(anaTemp_TumNor_ResScripts, "TumNor_somatic_commands.txt", sep="/"), row.names=F, col.names=T, quote=F, sep="\t")
  if(BASH)
  {
    TumNor_somatic_bash <- as.data.frame(TumNor_somatic[,8],stringsAsFactors = F)
    TumNor_somatic_bash <- sapply(TumNor_somatic_bash, as.character)
    TumNor_somatic_bash2 <- rbind(as.character('#!/bin/bash'), TumNor_somatic_bash)
    #View(TumNor_copynumber_bash2)
    write.table(TumNor_somatic_bash2, file=paste(anaTemp_TumNor_ResScripts, "TumNor_somatic_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  if(DEBUG)
  {
    View(TumNor_somatic)
  }
  
  
  
  
  print(paste(toString(Sys.time()),"Varscan_CopyNumber generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_TumNor_ResScripts,sep=" "))
  if(DEBUG)
  {
    View(TumNor_copynumber)
  }  
  
  TumNor_copycaller <- cbind(TumNor_Dr,
                             paste(paste(JAVA_HOME,'java',sep='/'),
                                   '-jar',
                                   VarScan,
                                   'copyCaller',
                                   paste(anaTemp_TumNor_copynumer_res,
                                         paste(paste('varscanCN',
                                                     sub("*.bam","",TumNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.bam","",TumNor_Dr[,2]),
                                                     "Tumor",
                                                     paste(TumNor_Dr[,4],
                                                           TumNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               'copynumber',
                                               sep='.'),
                                         sep='/'),
                                   '--output-file',
                                   paste(anaTemp_TumNor_copycaller_res,
                                         paste(paste('varscanCC',
                                                     sub("*.bam","",TumNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.bam","",TumNor_Dr[,2]),
                                                     "Tumor",
                                                     paste(TumNor_Dr[,4],
                                                           TumNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               'called',
                                               sep='.'),
                                         sep='/'),
                                   '--output-homdel-file',
                                   paste(anaTemp_TumNor_copycaller_res,
                                         paste(paste('varscanCC',
                                                     sub("*.bam","",TumNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.bam","",TumNor_Dr[,2]),
                                                     "Tumor",
                                                     paste(TumNor_Dr[,4],
                                                           TumNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'out',
                                               'called',
                                               'homdel',
                                               sep='.'),
                                         sep='/'),
                                   sep=' ')
  )
  colnames(TumNor_copycaller) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_copycaller_command')
  write.table(TumNor_copycaller[,8], file=paste(anaTemp_TumNor_ResScripts, "TumNor_copycaller_commands.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  if(BASH)
  {
    TumNor_copycaller_bash <- as.data.frame(TumNor_copycaller[,8],stringsAsFactors = F)
    TumNor_copycaller_bash <- sapply(TumNor_copycaller_bash, as.character)
    TumNor_copycaller_bash2 <- rbind(as.character('#!/bin/bash'), TumNor_copycaller_bash)
    #View(TumNor_mpileup_bash2)
    write.table(TumNor_copycaller_bash2, file=paste(anaTemp_TumNor_ResScripts, "TumNor_copycaller_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  
  print(paste(toString(Sys.time()),"Varscan_Copycaller generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_TumNor_ResScripts,sep=" "))
  if(DEBUG)
  {
    View(TumNor_copycaller)
  }
  
  
  TumNor_rm_mpileup <- cbind(TumNor_Dr,
                             paste('rm',
                                   paste(anaTemp_TumNor_mpileup_res,
                                         paste(paste('mpileup_res',
                                                     sub("*.cleaned.bam","",TumNor_Dr[,1]),
                                                     "Normal",
                                                     sub("*.cleaned.bam","",TumNor_Dr[,2]),
                                                     "Tumor",
                                                     paste(TumNor_Dr[,4],
                                                           TumNor_Dr[,6],
                                                           sep=''),
                                                     sep='_'),
                                               'mpileup',
                                               sep='.'),
                                         sep='/'),
                                   sep=' ')
  )
  colnames(TumNor_rm_mpileup) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','mpileup_rm_command') 
  write.table(TumNor_mpileup[,8], file=paste(anaTemp_TumNor_ResScripts, "TumNor_rm_mpileup_commands.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  if(BASH)
  {
    TumNor_rm_mpileup_bash <- as.data.frame(TumNor_rm_mpileup[,8],stringsAsFactors = F)
    TumNor_rm_mpileup_bash <- sapply(TumNor_rm_mpileup_bash, as.character)
    TumNor_rm_mpileup_bash2 <- rbind(as.character('#!/bin/bash'), TumNor_rm_mpileup_bash)
    #View(TumNor_copynumber_bash2)
    write.table(TumNor_rm_mpileup_bash2, file=paste(anaTemp_TumNor_ResScripts, "TumNor_rm_mpileup_bash.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
  }
  
  print(paste(toString(Sys.time()),"mpileup files removal script written at",anaTemp_TumNor_ResScripts,sep=" "))
  if(DEBUG)
  {
    View(TumNor_rm_mpileup)
  }
  
  TumNor_adj_logratio <- cbind(TumNor_Dr,
                               paste('cut -f1,2,3,7',
                                     paste(anaTemp_TumNor_copycaller_res,
                                           paste(paste('varscanCC',
                                                       sub("*.cleaned.bam","",TumNor_Dr[,1]),
                                                       "Normal",
                                                       sub("*.cleaned.bam","",TumNor_Dr[,2]),
                                                       "Tumor",
                                                       paste(TumNor_Dr[,4],
                                                             TumNor_Dr[,6],
                                                             sep=''),
                                                       sep='_'),
                                                 'out',
                                                 'called',
                                                 sep='.'),
                                           sep ='/'
                                     ),
                                     '>',
                                     paste(anaTemp_TumNor_adj_logratio,
                                           paste(paste('varscanCC',
                                                       sub("*.cleaned.bam","",TumNor_Dr[,1]),
                                                       "Normal",
                                                       sub("*.cleaned.bam","",TumNor_Dr[,2]),
                                                       "Tumor",
                                                       paste(TumNor_Dr[,4],
                                                             TumNor_Dr[,6],
                                                             sep=''),
                                                       sep='_'),
                                                 'adj.logratio',
                                                 sep='.'),
                                           sep='/'),
                                     sep=' ')
  )
colnames(TumNor_adj_logratio) <- c('Sample1','Sample2','UqBsAli_sample1','Gender1','UqBsAli_sample2','Gender2','Data_Ratio','Varscan_adj_logratio_commands')
write.table(TumNor_adj_logratio[,8], file=paste(anaTemp_TumNor_ResScripts, "TumNor_adj_logratio.txt", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
if(BASH)
{
  TumNor_adj_logratio_bash <- as.data.frame(TumNor_adj_logratio[,8],stringsAsFactors = F)
  TumNor_adj_logratio_bash <- sapply(TumNor_adj_logratio_bash, as.character)
  TumNor_adj_logratio_bash2 <- rbind(as.character('#!/bin/bash'), TumNor_adj_logratio_bash)
  #View(TumNor_copynumber_bash2)
  write.table(TumNor_adj_logratio_bash2, file=paste(anaTemp_TumNor_ResScripts, "TumNor_adj_logratio_commands.sh", sep="/"), row.names=F, col.names=F, quote=F, sep="\t")
}

print(paste(toString(Sys.time()),"Adjusted logratio files generating command file for Normal-Normal Pairs Created and Saved to",anaTemp_TumNor_ResScripts,sep=" "))
if(Com_Scr)
{
  adj_logratio_cmnds <<- as.data.frame(TumNor_adj_logratio[,8],stringsAsFactors = F, row.names= NULL, optional = F, colnames = F)
  adj_logratio_cmnds[,1] <<- sapply(adj_logratio_cmnds[,1], as.character)
}

if(DEBUG)
{
  View(TumNor_adj_logratio)
}
  
  
}


##########################################################################################################################################################################################





