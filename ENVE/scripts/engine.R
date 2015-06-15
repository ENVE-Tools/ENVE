###########~~~~~~~~~Download Required Package~~~~~~~~~~~~~~~~~########################
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

####################################################################################################



##########################~~~~~~~~~~~~~~~~~MultiPlot Function~~~~~~~~~~~~~~~~#######################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    #pushViewport(viewport(layout = grid.layout(ncol(layout), nrow(layout))))
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      #print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
      print(plots[[i]], vp = viewport(layout.pos.col = matchidx$col,layout.pos.row = matchidx$row))
    }
  }
}

###################################~~~~~~~~~~multiplot-function-ends~~~~~~~~~~~~~~#################################################


###################################~~~~~~~~~~~Silhoutte-Funtion~~~~~~~~~~~~~~~~~~~#################################################

Get_Sil_Index <-   function(freq_mat, seg_index) 
  {
    dis_vals = 1-cor(t(freq_mat), method="pearson")
    upper_segvals_dis_vals_matrix = dis_vals[(seg_index+1):(dim(dis_vals)[1]), (seg_index+1):(dim(dis_vals)[1])]
    low_tri_upper_segvals_dis_vals_matrix_index = which(lower.tri(upper_segvals_dis_vals_matrix))
    mean_upper_segvals_distances = mean(upper_segvals_dis_vals_matrix[low_tri_upper_segvals_dis_vals_matrix_index])
    lower_segvals_dis_vals_matrix = dis_vals[1:(seg_index-1), 1:(seg_index-1)]
    low_tri_lower_segvals_dis_vals_matrix_index = which(lower.tri(lower_segvals_dis_vals_matrix))
    mean_lower_segvals_distances = mean(lower_segvals_dis_vals_matrix[low_tri_lower_segvals_dis_vals_matrix_index])
    lower_segvals_index = 1:(seg_index-1)
    upper_segvals_index = (seg_index+1):(dim(dis_vals)[1])
    all_lower_upper_pairs = expand.grid(lower_segvals_index, upper_segvals_index)
    sum_lower_upper_pairs_disvals = 0
    for (i in 1:length(all_lower_upper_pairs[,1]))
      {
        sum_lower_upper_pairs_disvals = sum_lower_upper_pairs_disvals+dis_vals[all_lower_upper_pairs[i,1], all_lower_upper_pairs[i,2]]
      }
    mean_lower_upper_pairs_disvals = sum_lower_upper_pairs_disvals/length(all_lower_upper_pairs[,1])
    si_up = (mean_lower_upper_pairs_disvals-mean_upper_segvals_distances)/max(mean_upper_segvals_distances, mean_lower_upper_pairs_disvals)
    si_low = (mean_lower_upper_pairs_disvals-mean_lower_segvals_distances)/max(mean_lower_segvals_distances, mean_lower_upper_pairs_disvals)
    avg_si = mean(si_up, si_low)
    return(avg_si)
  }
######################################################################################################################################



#############################~~~~~~~~~~~~~~~~~Directory_Creator~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####################################
dir_create <- function()
{
  ana <<- paste(enveHome,"analysis",sep='/')
  #####Analysis files#########
  anaPath <<- paste(ana,(paste("Analysis",toString(strftime(Sys.time(),format="%d_%m_%Y_%H_%M")),sep='_')),sep ="/")
  #anaInp <<- paste(anaPath,"Input",sep='/')
  #anaInpBamsNorm <<-paste(anaInp,"Normal_BAMs",sep='/')
  #anaInpBamsTum <<-paste(anaInp,"Tumor_BAMs",sep='/')
  anaTemp <<-paste(anaPath,"temp",sep='/')
  anaTempLR <<-paste(anaTemp,"LogRatio_CNA",sep='/')
  anaTempLRNN <<-paste(anaTempLR,"NormalNormal",sep='/')
  anaRes <<- paste(anaPath,"Results",sep='/')
  anaEVDPVal_AS <<- paste(anaTemp,"EVDPVal_AllSamps",sep='/')
  
  
  
  anaTempLRNNres <<-paste(anaTempLRNN,"Logratio_Results",sep='/')
  anaTempLRNN_MODE_CRCT <<- paste(anaTempLRNN,"Mode_Corrected_files",sep='/')
  anaTempLRNN_CBS_GC_crtd <<-paste(anaTempLRNN,"CBS_Segments_GC_Corrected_CDSFilt",sep='/')
  anaTempLRNN_NorNor_SegMeans_CDSFilt <<-paste(anaTempLRNN,"chr_ALLNormNorm",sep='/')
  anaTempLRNN_NorNor_Pos <<-paste(anaTempLRNN,"chr_ALLNormNorm_Pos",sep='/')
  anaTempLRNN_NorNor_Neg <<-paste(anaTempLRNN,"chr_ALLNormNorm_Neg",sep='/')
  anaTempLRNN_OC_CDSFilt <<- paste(anaTempLRNN,"LogRatioCC_OC_FiltCDS",sep='/')
  anaTempLRNN_Tiff_output<<- paste(anaTempLRNN,"TIFF_output",sep='/')
  anaTempLRNN_EVD_Cutoff <<-paste(anaTempLRNN,"EVD_cutoff",sep='/')
  anaTempLRNN_Pos_Tiff_output<<- paste(anaTempLRNN,"Pos_TIFF_output",sep='/')
  anaTempLRNN_Pos_EVD_Cutoff <<-paste(anaTempLRNN,"Pos_EVD_cutoff",sep='/')
  anaTempLRNN_Neg_Tiff_output<<- paste(anaTempLRNN,"Neg_TIFF_output",sep='/')
  anaTempLRNN_Neg_EVD_Cutoff <<-paste(anaTempLRNN,"Neg_EVD_cutoff",sep='/')
  
  
  
  anaTempLRTN <<-paste(anaTempLR,"TumorNormal",sep='/')
  anaTempLRTN_MODE_CRCT <<- paste(anaTempLRTN,"Mode_Corrected_files",sep='/')
  anaTempLRTNres <<-paste(anaTempLRTN,"Logratio_Results",sep='/')
  anaTempLRTN_CBS_GC_crtd <<-paste(anaTempLRTN,"CBS_Segments_GC_Corrected_CDSFilt",sep='/')
  anaTempLRTN_TumNor_SegMeans_CDSFilt <<-paste(anaTempLRTN,"Norm_vs_Tum_SegMeans_CDSFilt",sep='/')
  anaTempLRTN_OC_CDSFilt <<- paste(anaTempLRTN,"LogRatioCC_OC_FiltCDS",sep='/')
  
  
  #anaTempLRTN_allchr_CDSFilt <<- paste(anaTempLRNN,"chr_allTNSamps_SegMeans_CDSFilt",sep='/')
  ################## Insert more dirs######################  
  dirs <- c(enveHome,
            scriptsPath,
            supFiles,
            ana,
            anaPath,
            #anaInp,
            #anaInpBamsNorm,
            #anaInpBamsTum,
            anaTemp,
            anaTempLR,
            anaTempLRNN,
            anaTempLRNNres,
            anaTempLRNN_MODE_CRCT,
            anaTempLRNN_CBS_GC_crtd,
            anaTempLRNN_NorNor_SegMeans_CDSFilt,
            anaTempLRNN_NorNor_Pos,
            anaTempLRNN_NorNor_Neg,
            anaTempLRNN_Pos_Tiff_output,
            anaTempLRNN_Pos_EVD_Cutoff,
            anaTempLRNN_Neg_Tiff_output,
            anaTempLRNN_Neg_EVD_Cutoff,
            anaTempLRNN_OC_CDSFilt,
            anaTempLRNN_Tiff_output,
            anaTempLRNN_EVD_Cutoff,
            anaTempLRTN,
            anaTempLRTNres,
            anaTempLRTN_MODE_CRCT,
            anaTempLRTN_CBS_GC_crtd,
            anaTempLRTN_TumNor_SegMeans_CDSFilt,
            anaTempLRTN_OC_CDSFilt,
            anaRes,
            anaEVDPVal_AS
  )
  
  dir2 <- as.data.frame(cbind(c('enveHome',
                                'scriptsPath',
                                'supFiles',
                                'ana',
                                'anaPath',
                                'anaTemp',
                                'anaTempLR',
                                'anaTempLRNN',
                                'anaTempLRNNres',
                                'anaTempLRNN_MODE_CRCT',
                                'anaTempLRNN_CBS_GC_crtd',
                                'anaTempLRNN_NorNor_SegMeans_CDSFilt',
                                'anaTempLRNN_NorNor_Pos',
                                'anaTempLRNN_NorNor_Neg',
                                'anaTempLRNN_Pos_Tiff_output',
                                'anaTempLRNN_Pos_EVD_Cutoff',
                                'anaTempLRNN_Neg_Tiff_output',
                                'anaTempLRNN_Neg_EVD_Cutoff',
                                'anaTempLRNN_OC_CDSFilt',
                                'anaTempLRNN_Tiff_output',
                                'anaTempLRNN_EVD_Cutoff',
                                'anaTempLRTN',
                                'anaTempLRTNres',
                                'anaTempLRTN_MODE_CRCT',
                                'anaTempLRTN_CBS_GC_crtd',
                                'anaTempLRTN_TumNor_SegMeans_CDSFilt',
                                'anaTempLRTN_OC_CDSFilt',
                                'anaRes',
                                'anaEVDPVal_AS'
  ), dirs))
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
}

################################################################################################################

chrs = as.matrix(c('chr1',
                   'chr2',
                   'chr3',
                   'chr4',
                   'chr5',
                   'chr6',
                   'chr7',
                   'chr8',
                   'chr9',
                   'chr10',
                   'chr11',
                   'chr12',
                   'chr13',
                   'chr14',
                   'chr15',
                   'chr16',
                   'chr17',
                   'chr18',
                   'chr19',
                   'chr20',
                   'chr21',
                   'chr22',
                   'chrX',
                   'chrY'))



chr_proc <- function()
{
  chr_lengths = read.table(paste(supFiles,"chr_lengths_hg19.txt",sep="/"), sep="\t")
  chr_lengths = as.matrix(chr_lengths)
  return(chr_lengths)
}

##################################################################################################################################3

reqd_files_func <- function(x)
{
  samp_info <-as.data.frame(read.delim(file = x,header=T,sep ='\t' ),row.names=NULL,optional = F, stringsAsFactors = F)
  i <- sapply(samp_info, is.factor)
  samp_info[i] <- lapply(samp_info[i], as.character)
  colnames(samp_info) <- c("Sample1","Sample2","Gender_Samp1","Gender_Samp2")
  reqd_files <- cbind(paste(samp_info[,'Sample1'],samp_info[,'Sample2'],sep="_"),paste(samp_info[,'Gender_Samp1'],samp_info[,'Gender_Samp2'],sep=""))
  colnames(reqd_files) <- c("Sample", "Gender")
  return(reqd_files)
}
#####################################################################################################################################
##################### ~~~~~~~~~~~~~~~~~~~~~~Mode-Correction Related Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#########################
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = subset(x, !is.na(x))
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

mode_correction <- function(x,y)
{
  setwd(x)
  input.files = list.files(pattern=".GC_COR_adj_logratio")
  input.files = gsub(".GC_COR_adj_logratio","",input.files)
  input.files <- as.data.frame(input.files,row.names=NULL, optional=F,stringsAsFactors = F)
  colnames(input.files) <- "Sample"
  AvailSamps <- merge(input.files,reqd_files, by.input.files='Sample',by.reqd_files='Sample',all=F)
  for( i in 1:length(AvailSamps[,1]))
  {
    file <- read.table(paste(x,paste(AvailSamps[i,1],".GC_COR_adj_logratio",sep=""),sep="/"),sep='\t',header=T)
    mode <- Mode(file[,4])
    mod_crcted_file <- cbind(file, (file[,4]-mode))
    write.table(mod_crcted_file,paste(y,paste(AvailSamps[i,1],".GC_COR_adj_logratio",sep=""),sep="/"),sep="\t",row.names=F,col.names=F,quote=F)
  }
  
}


######################################################################################################################################

###################################~~~~~~~~~~~~~FILT CDS FUNCTION~~~~~~~~~~~~~~~~~~~~~~~##############################################
filtCDS <- function(x,y)
{
  setwd(x)
  input.files = list.files(pattern=".GC_COR_adj_logratio")
  input.files = gsub(".GC_COR_adj_logratio","",input.files)
  input.files <- as.data.frame(input.files,row.names=NULL, optional=F,stringsAsFactors = F)
  colnames(input.files) <- "Sample"
  AvailSamps <- merge(input.files,reqd_files, by.input.files='Sample',by.reqd_files='Sample',all=F)
  
  ##### FILTER WINDOWS TO THOSE THAT FALL WITHIN CODING REGIONS ######
  if(Whole_Exome)
  {
  for (i in 1:length(AvailSamps[,1])) {
    sys_cmd = paste("awk 'NR>1' ", paste(AvailSamps[i,1],".GC_COR_adj_logratio",sep=""), " | sed 's/ //g' | ",IntersectBED_path," -a stdin -b ",paste(supFiles,"RefSeqGTF/RefSeq_hg19_Feb2009.gtf.CDS.nospace.bed.txt",sep='/')," -f 0.6 -u > ",paste(y,paste(AvailSamps[i,1],".FiltCDS",sep=""),sep='/'), sep="")
    system(sys_cmd)
  }
  }else{
    for(i in 1:length(AvailSamps[,1]))
    {
      file.copy(paste(x,AvailSamps[i,1],sep='/'),y, overwrite = recursive, recursive = TRUE, copy.mode = TRUE)
    }
    setwd(y)
    for(i in 1:length(AvailSamps[,1]))
    {
      file.rename(paste(AvailSamps[i,1],'.GC_COR_adj_logratio',sep='.'), paste(AvailSamps[i,1],'FiltCDS',sep='.'))
    }
  }
  
  if(GISTIC)
  {
    setwd(y)
    input.files = list.files(pattern=".FiltCDS")
    input.files = gsub(".FiltCDS","",input.files)
    input.files <- as.data.frame(input.files,row.names=NULL, optional=F,stringsAsFactors = F)
    colnames(input.files) <- "Sample"
    AvailSamps <- merge(input.files,reqd_files, by.input.files='Sample',by.reqd_files='Sample',all=F)
    
    markers_out = NULL
    for (i in 1:length(AvailSamps[,1])){
      vvi = read.table(paste(AvailSamps[i,1],".FiltCDS",sep=""), sep="\t", header=F)
      if (pmatch("chr", vvi[1,1])==1){
        outi_chr = apply(as.matrix(vvi[,1]), 1, function(x) unlist(strsplit(x, split="chr"))[2])
      } else {
        outi_chr = as.matrix(vvi[,1])
      }
      markers_out = c(markers_out, paste(outi_chr, vvi[,2], sep="_"))
    }
    markers_out = unique(markers_out)
    markers_out_split = matrix(data=NA, ncol = 2, nrow = length(markers_out))
    for (i in 1:length(markers_out)){
      markers_out_split[i,] =  unlist(strsplit(markers_out[i], split="_"))  
    }
    dd = data.frame(chr = markers_out_split[,1], location = as.numeric(markers_out_split[,2]))
    ix <- sapply(dd, is.factor)
    dd[ix] <- lapply(dd[ix], as.character)
    dd_srt = dd[order(dd$chr,dd$location),]
    dd_srt_names = paste(dd_srt$chr,dd_srt$location, sep="_")
    rownames(dd_srt) = dd_srt_names
    
    write.table(dd_srt, file=paste(anaRes,"MarkerFile_for_GISTIC.txt",sep="/"), row.names=T, col.names=F, quote=F, sep="\t", eol = "\n")
  }
  ##### END OF FILTER WINDOWS TO THOSE THAT FALL WITHIN CODING REGIONS ######
}

#################################################################################################################################
###################################################CBS-SEGMENTATION##############################################################
CBS_seg_samp <- function(x,y)
{
  setwd(x)
  input.files = list.files(pattern=".FiltCDS")
  input.files = gsub(".FiltCDS","",input.files)
  input.files <- as.data.frame(input.files,row.names=NULL, optional=F,stringsAsFactors = F)
  colnames(input.files) <- "Sample"
  AvailSamps <- merge(input.files,reqd_files, by.input.files='Sample',by.reqd_files='Sample',all=F)
  ##### CBS Segmentation Per SAMPLE ######
  
  for (i in 1:length((AvailSamps[,1]))) {
    print(i)
    filcna = paste(AvailSamps[i,1], ".FiltCDS", sep="")
    #sampcna = paste(unlist(strsplit(filcna, split="\\_|\\."))[c(1,2)], collapse="_")
    cn <- read.table(filcna,header=F)
    cn = as.matrix(cn)
    cn = str_replace_all(cn, " ", "")
    CNA.object <-CNA(genomdat = as.numeric(cn[,4]), chrom = cn[,1], maploc = as.numeric(cn[,2]), data.type = 'logratio', sampleid=AvailSamps[i,1])
    CNA.smoothed <- smooth.CNA(CNA.object)
    segs <- segment(CNA.smoothed, verbose=0, min.width=2, undo.splits = "sdundo", undo.SD=3)
    segs2 = segs$output
    write.table(segs2[,2:6], file=paste(y,paste(AvailSamps[i,1],"CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t")
  }
}

###################################################################################################################################

Com_samp_perchr <- function(x,y)
{
  setwd(x)
  input.files = list.files(pattern="_CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt")
  input.files = gsub("_CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt","",input.files)
  input.files <- as.data.frame(input.files,row.names=NULL, optional=F,stringsAsFactors = F)
  colnames(input.files) <- "Sample"
  AvailSamps <- merge(input.files,reqd_files, by.input.files='Sample',by.reqd_files='Sample',all=F)
  ##### BEGIN Combine All Samples per chromosome ######
  for (j in 1:length(chrs)) {
    chr_out = NULL
    for (i in 1:length(AvailSamps[,1])) {
      filcna = AvailSamps[i,1]
      filcna2 <<- AvailSamps[i,2]
      if (chrs[j] != "chrX" & chrs[j] != "chrY"){
        sampcna = AvailSamps[i,1]
        segs.all.chr <- read.table(paste(sampcna, "CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),header=F, sep="\t")
        ind_chr = which(segs.all.chr[,1]==chrs[j])
        samp_chr_out =  cbind(matrix(data=sampcna, nrow = length(ind_chr), ncol=1), segs.all.chr[ind_chr,])
        chr_out = rbind(chr_out, samp_chr_out)
      }else if(chrs[j] =="chrX"){
        if(filcna2 == "MaleMale" | filcna2 == "FemaleFemale"){
          sampcna = AvailSamps[i,1]
          segs.all.chr <- read.table(paste(sampcna, "CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),header=F, sep="\t")
          ind_chr = which(segs.all.chr[,1]==chrs[j])
          samp_chr_out =  cbind(matrix(data=sampcna, nrow = length(ind_chr), ncol=1), segs.all.chr[ind_chr,])
          chr_out = rbind(chr_out, samp_chr_out)
          }
      }else if(chrs[j] == "chrY"){
        if(filcna2 == "MaleMale"){
          sampcna =  AvailSamps[i,1]
          segs.all.chr <- read.table(paste(sampcna, "CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),header=F, sep="\t")
          ind_chr = which(segs.all.chr[,1]==chrs[j])
          samp_chr_out =  cbind(matrix(data=sampcna, nrow = length(ind_chr), ncol=1), segs.all.chr[ind_chr,])
          chr_out = rbind(chr_out, samp_chr_out)
        }
      }
    }
    write.table(chr_out, file=paste(y,paste(chrs[j],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t", eol="\n")
    #write.table(chr_out[which(as.numeric(chr_out[,6])>0)], file=paste(anaTempLRNN_NorNor_Pos,paste(chrs[j],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t", eol="\n")
    #write.table(chr_out[which(as.numeric(chr_out[,6])<0)], file=paste(anaTempLRNN_NorNor_Neg,paste(chrs[j],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t", eol="\n")
  }
  ##### END Combine All Samples per chromosome ######
}
##################################################################################################################################
chr_pos_neg_sep <- function(x)
{
  setwd(x)
  input.files = list.files(pattern="_AllNormNormSamps_SegMeans_CDSFilt.txt")
  input.files = gsub("_AllNormNormSamps_SegMeans_CDSFilt.txt","",input.files)
  input.files <- as.data.frame(input.files,row.names=NULL, optional=F,stringsAsFactors = F)
  colnames(input.files) <- "Sample"
  AvailSamps <- intersect(as.matrix(input.files),as.matrix(chrs))
  fil <- NULL
  for(chrind in 1:length(AvailSamps))
  {  
    file=paste(x,paste(AvailSamps[chrind],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/')
    read_test <- file.info(file)$size
    if(read_test!=0)
    {
      fil <- read.delim(file=paste(x,paste(AvailSamps[chrind],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/'),header = F,sep="\t")
      write.table(fil[which(as.numeric(fil[,6])>0),], file=paste(anaTempLRNN_NorNor_Pos,paste(AvailSamps[chrind],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t", eol="\n")
      write.table(fil[which(as.numeric(fil[,6])<0),], file=paste(anaTempLRNN_NorNor_Neg,paste(AvailSamps[chrind],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t", eol="\n")
    } 
  }
}
###################################################################################################################################

Nor_EVD_calc <- function(x,y,z)
{
  setwd(x)
  norm_norm_path <- x
  tiff_output_path = z
  EVD_calc_path = y
 
  EVD_cutoff_Chr = NULL
  pileup_windows_overlap_segments_Chr <<- NULL
  pileup_windows_overlap_counts_Chr <<- NULL
  fraction_nonzero_coverage_AllChrs <<- NULL
  entropy_coverage_AllChrs <<- NULL  
  chromosomal_plots <- NULL
  entropy_plots <- NULL
  input.files = list.files(pattern="_AllNormNormSamps_SegMeans_CDSFilt.txt")
  input.files = gsub("_AllNormNormSamps_SegMeans_CDSFilt.txt","",input.files)
  input.files <- as.data.frame(input.files,row.names=NULL, optional=F,stringsAsFactors = F)
  colnames(input.files) <- "Sample"
  AvailSamps <- intersect(as.matrix(input.files),as.matrix(chrs))
  chrs <- AvailSamps
  for(chrind in 1:(length(chrs)))
    {
      print(paste("chr : ",chrind,sep=' '))
      chr_out_n = NULL
      inds_large_segs_n = NULL    
      chr_out_n = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
      inds_large_segs_n = which(chr_out_n[,5] >=min_seg)
      chr_out_n_sizeseg = chr_out_n[inds_large_segs_n, ]
      chr_len = as.numeric(chr_lengths[which(chr_lengths[,1] == chrs[chrind]),2])
      #### create non-verlapping windows of 5kb length to cover the entire chromosome ####
      win_size = 10000
      epsilon = 1e-16
      num_wins = chr_len%/%win_size
      print(paste("num_wins",num_wins,sep=' '))
      chr_wind_ends = seq(from=0, to = chr_len, by = win_size)
      chr_wind_ends = c(chr_wind_ends[-1],chr_len) 
      chr_wind_begins = seq(from=1, to = chr_len, by = win_size)
      chr_winds <<- cbind(chr_wind_begins, chr_wind_ends)
      chr_winds_IR = IRanges(chr_wind_begins, chr_wind_ends)
      chr_winds_overlap_entropy = matrix(data=NA, nrow=1, ncol = length(segval_windows))
      chr_winds_overlap_freq_matrix = matrix(data=NA, nrow=length(segval_windows), ncol = dim(chr_winds)[1])
      chr_winds_overlap_cnt_matrix = matrix(data=NA, nrow=length(segval_windows), ncol = dim(chr_winds)[1])
      for (winind in 1:length(segval_windows)){
        segs_in_wind = chr_out_n_sizeseg[which(abs(chr_out_n_sizeseg[,6])>=segval_windows[winind]),]
        segs_in_wind_IR = IRanges(segs_in_wind[,3], segs_in_wind[,4])
        chr_winds_overlap_cnt = countOverlaps(chr_winds_IR,segs_in_wind_IR)
        chr_winds_overlap_cnt_matrix[winind,] = chr_winds_overlap_cnt
        chr_winds_overlap_freq = chr_winds_overlap_cnt/sum(chr_winds_overlap_cnt)
        chr_winds_overlap_freq_matrix[winind,] = chr_winds_overlap_freq
        chr_winds_overlap_entropy[winind] = -1*sum(chr_winds_overlap_freq*log(chr_winds_overlap_freq), na.rm=T)
        tmp123 = cbind(chrs[chrind], chr_winds, segval_windows[winind], chr_winds_overlap_cnt)
        pileup_windows_overlap_counts_Chr <<- rbind(pileup_windows_overlap_counts_Chr, tmp123)
        if (dim(segs_in_wind)[1] >0){
          segx0 = segs_in_wind[,3]
          segx1 = segs_in_wind[,4]
          segy0 = segval_windows[winind]
          segy1 = segy0      
          if (length(segs_in_wind[,1]) > 1){
            tmp1234 = cbind(as.matrix(segs_in_wind[,c(2,3,4)]),segval_windows[winind])
          } else {
            tmp1234 = t(as.matrix(c(as.matrix(segs_in_wind[,c(2,3,4)]),segval_windows[winind])))
          }
          pileup_windows_overlap_segments_Chr <- rbind(pileup_windows_overlap_segments_Chr, tmp1234)
        }
      }
  
      #identify minimum segval means that have at least 10 chromosomal windows with non-zero frequency 
      #this gives us the maximum mean segval for which we want to analyze the distribution of reads
      nonNA_segval_windows_index <<- which(apply(chr_winds_overlap_freq_matrix,1,function(x) length(x)-length(which(is.na(x))))>=10)
      nonNA_chr_winds_overlap_freq_matrix <<- chr_winds_overlap_freq_matrix[nonNA_segval_windows_index,]
      chr_winds_overlap_entropy_nonNA <<- chr_winds_overlap_entropy[nonNA_segval_windows_index]
      nonNA_segval_windows <<- segval_windows[nonNA_segval_windows_index]
    
      
      #identify fraction of chromosome that has a non-zero frequency of reads aligning to it for every meansegval cutoff 
      frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix <<- apply(nonNA_chr_winds_overlap_freq_matrix,1, function(x)  length(which(x!=0))/length(x))
      
      ###output fraction of non-zero coverage for chromosome
      fraction_nonzero_coverage_AllChrs <<- rbind(fraction_nonzero_coverage_AllChrs, c(chrs[chrind], frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix))
      entropy_coverage_AllChrs <<- rbind(fraction_nonzero_coverage_AllChrs, c(chrs[chrind], chr_winds_overlap_entropy))
      
      ### pick the seg-val cutoff where the genomic coverage drop is maximum  -- WE FOUND THIS BEFORE 
      winind_sig = which.min(diff(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix))+1
      diff_frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix = NULL
      diff_frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix <<- c(NA, -diff(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix))
    
      ### if chr fraction at maximum segval is greater than 25%, use all segvals for EVD ###
      if (frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix[length(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix)] >= 0.25){
        droppoint_segval_value = max(nonNA_segval_windows) + 0.005
      } else {
        droppoint_segval_index = which.min(diff(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix))+1
        
        ### begin estimate Silhouette index for cutoff point to ensure that segment distribution above cutoff is dissimilar to segment distribution below ###
        
        si_test1 <- try(Get_Sil_Index(nonNA_chr_winds_overlap_freq_matrix, droppoint_segval_index))
        if(class(si_test1) != "try-error")
        {
          si_1 = Get_Sil_Index(nonNA_chr_winds_overlap_freq_matrix, droppoint_segval_index)
        }
        si_test2 <- try(Get_Sil_Index(nonNA_chr_winds_overlap_freq_matrix, droppoint_segval_index-1))
        if(class(si_test2) != "try-error")
        {
          si_2 = Get_Sil_Index(nonNA_chr_winds_overlap_freq_matrix, droppoint_segval_index-1)
        }
        if(!is.na(si_2))
        {
          if (!is.na(si_1))
              {
                if (si_2 > si_1){
                droppoint_segval_index = droppoint_segval_index-1 
                si_droppoint = si_2
              } else {
                droppoint_segval_index = droppoint_segval_index 
                si_droppoint = si_1
              }    
          }else if (is.na(si_1) & si_2 > 0) {
              droppoint_segval_index = droppoint_segval_index-1 
              si_droppoint = si_2
            } else {
              droppoint_segval_index = droppoint_segval_index 
              si_droppoint = si_1
            }
        }
        
          droppoint_segval_value = nonNA_segval_windows[droppoint_segval_index]+0.005
        }

      plot_chr <- ggplot() +
                  theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(), 
                        panel.background = element_blank(),
                        panel.border = element_rect(colour="black", fill=NA),
                        axis.line = element_line(colour = "black"),
                        axis.text.x = element_text(size=30,angle=0, colour="black", vjust = 0.5),
                        axis.title.y = element_blank(), 
                        axis.text.y = element_text(size=20, colour="black", hjust =0.5),
                        axis.title.x=element_blank()
                        )+
                      geom_point(aes(x=nonNA_segval_windows, y = frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix), size = 5, colour="black")+
                      scale_x_continuous(breaks = c(seq(0,max(nonNA_segval_windows),0.2)),labels =c(as.character(seq(0,max(nonNA_segval_windows),0.2))))+
                      scale_y_continuous(breaks = c(seq(0,ceiling(max(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix)),0.2)),
                                        labels =c(as.character(seq(0,ceiling(max(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix)),0.2))))+
                      coord_cartesian(xlim = c(0,max(nonNA_segval_windows)+0.02), ylim =c(0,max(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix)+0.02))+
                  ggtitle(chrs[chrind]) + 
                  theme(plot.title = element_text(lineheight=.8, face="bold", size = 48))+
                  geom_vline(xintercept=droppoint_segval_value, linetype = "dotted", lwd = 1, colour="black")
                  
      
chromosomal_plots[[chrind]] <- plot_chr  
   
plot_entropy <- ggplot() +
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      panel.background = element_blank(),
                      panel.border = element_rect(colour="black", fill=NA),
                      axis.line = element_line(colour = "black"),
                      axis.text.x = element_text(size=30,angle=0, colour="black", vjust = 0.5),
                      axis.title.y = element_blank(), 
                      axis.text.y = element_text(size=20, colour="black", hjust =0.5),
                      axis.title.x=element_blank()
                )+
                geom_point(aes(x=nonNA_segval_windows, y = chr_winds_overlap_entropy_nonNA), size = 5, colour="black")+
                scale_x_continuous(breaks = c(seq(0,max(nonNA_segval_windows),0.2)),labels =c(as.character(seq(0,max(nonNA_segval_windows),0.2))))+
                scale_y_continuous(breaks = c(seq(0,ceiling(max(chr_winds_overlap_entropy_nonNA)),1)),
                                   labels =c(as.character(seq(0,ceiling(max(chr_winds_overlap_entropy_nonNA)),1))))+
                coord_cartesian(xlim = c(0,max(nonNA_segval_windows)+0.02), ylim =c(0,max(chr_winds_overlap_entropy_nonNA)+1))+
                ggtitle(chrs[chrind]) + 
                theme(plot.title = element_text(lineheight=.8, face="bold", size = 48))+
                geom_vline(xintercept=droppoint_segval_value, linetype = "dotted", lwd = 1, colour="black") 

  entropy_plots[[chrind]] <- plot_entropy
    
    EVD_cutoff_Chr = rbind(EVD_cutoff_Chr, c(chrs[chrind], droppoint_segval_value))
  }
  
  
tiff(filename = paste(tiff_output_path, "POS_CHROMOSOME_COVERAGE_AllChrs_Thresholds.tiff", sep=""), width = 8, height = 11, units = "in", res = 300, pointsize = 12, compression = c("none"),bg = "white")
multiplot(plotlist=chromosomal_plots, cols = 4)
dev.off()

tiff(filename = paste(tiff_output_path, "POS_ENTROPY_AllChrs_Thresholds.tiff", sep=""), width = 8, height = 11, units = "in", res = 300, pointsize = 12, compression = c("none"),bg = "white")
multiplot(plotlist=entropy_plots, cols = 4)
dev.off()



  colnames(EVD_cutoff_Chr) = c("Chromosome", "SegValCutoffForEVD")
  write.table(EVD_cutoff_Chr, file=paste(EVD_calc_path,paste("EVD_Chr_Cutoffs","txt", sep="."),sep='/'), sep="\t", quote=F, row.names=F, col.names=T, eol="\n")
  
  ###################################################################################################################################
  
  

  EVD_cutoff_Chr = read.table(paste(EVD_calc_path,paste("EVD_Chr_Cutoffs","txt", sep="."),sep='/'), header=T)
  EVD_cutoff_Chr = as.matrix(EVD_cutoff_Chr)
  
  
  EVD_params_per_chr = NULL
  
  for (chrind in 1:length(chrs))  
    {
      chr_out_n = NULL
      inds_large_segs_n = NULL
      if (chrs[chrind] != "chrX" & chrs[chrind] != "chrY")  {
        chr_out_n = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
      } else if (chrs[chrind] == "chrX"){
        chr_out_nF = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
        chr_out_nM = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
        chr_out_n = rbind(chr_out_nF, chr_out_nM)
      } else if (chrs[chrind] == "chrY") {
        chr_out_n = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
      }
      inds_large_segs_n = which(chr_out_n[,5] >=min_seg)
      chr_out_n_sizeseg = chr_out_n[inds_large_segs_n, ]
      EVD_cutoff = as.numeric(EVD_cutoff_Chr[which(EVD_cutoff_Chr[,1] == chrs[chrind]),2])
      chr_out_n_sizeseg_for_EVD = chr_out_n_sizeseg[which(abs(chr_out_n_sizeseg[,6])<=EVD_cutoff),]
      chr_len = as.numeric(chr_lengths[which(chr_lengths[,1] == chrs[chrind]),2])
      
      samp_names = unique(chr_out_n_sizeseg_for_EVD[,1])
      training_vals_abs_max = NULL
      for (i in 1:length(samp_names)){
        max_abs_val = max(abs(chr_out_n_sizeseg_for_EVD[which(chr_out_n_sizeseg_for_EVD[,1]==samp_names[i]),6]))
        if (max_abs_val > 0) training_vals_abs_max = c(training_vals_abs_max, max_abs_val)
      }
      evd_vals_abs_max=gevFit(x = training_vals_abs_max, block = 1, type = "pwm")
      save(evd_vals_abs_max, file=paste(EVD_calc_path,paste("EVDParams_", chrs[chrind], ".RData", sep=""),sep='/'))
      EVD_params_per_chr = rbind(EVD_params_per_chr, c(chrs[chrind], evd_vals_abs_max@fit$par.ests))
  }
  colnames(EVD_params_per_chr) = c("Chromosome", "AbsMax_xi", "AbsMax_mu", "AbsMax_beta")
  
  write.table(EVD_params_per_chr, file=paste(EVD_calc_path, "EVD_Params_AllChrs.txt", sep="/"), sep="\t", row.names=F, col.names=T, eol = "\n", quote=F)
  ##################################################################################################################################################
if(T)
{
  pileup_cnts_file <- pileup_windows_overlap_counts_Chr
  pileup_segs_data <- pileup_windows_overlap_segments_Chr
  pileup_segs_data = str_replace_all(as.matrix(pileup_segs_data), " ", "")
  pileup_segs_data_forplot = pileup_segs_data
  colnames(pileup_segs_data_forplot) = c("Chromosome", "Start", "End", "LogRatio")
  evd_cutoff_data <- EVD_cutoff_Chr
  evd_cutoff_data = str_replace_all(as.matrix(evd_cutoff_data), " ", "")
  
  
  #abs pos is 1-base - so need to convert from chrpos which is 0-base - so need to figure out what to add 
  #chrlen is basically giving us the chrpos of final base in chr - so for chr1, len is 249250621 but that means there are 249250621+1 bases in chr1
  #so for all chrpos in chr1, need to add 1 base to account for 0-base in chr1; similarly for chr2, need to add 1 base for chr1, length of chr1 and 1 base for chr2; and so on... 
  chr_add_for_abspos = NULL
  chr_lengths_forabspos = as.matrix(c(0, as.numeric(chr_lengths[,2])))
  for (i in 1:length(chr_lengths[,1])){
    if(chr_lengths[i,1] != "chrX" & chr_lengths[i,1] != "chrY") chr_add = sum(as.numeric(chr_lengths_forabspos[1:i,1]))+as.numeric(unlist(strsplit(chr_lengths[i,1],split="chr"))[2])
    if(chr_lengths[i,1] == "chrX") chr_add = sum(as.numeric(chr_lengths_forabspos[1:i,1]))+23
    if(chr_lengths[i,1] == "chrY") chr_add = sum(as.numeric(chr_lengths_forabspos[1:i,1]))+24
    chr_add_for_abspos = rbind(chr_add_for_abspos, c(chr_lengths[i,1], chr_add))
  }
  
  
  
  chr_cent_telo = read.delim(paste(supFiles,"hg19_chr_telomere_centromere.txt",sep='/'), sep="\t", header=T)
  chr_cent = chr_cent_telo[which(chr_cent_telo[,"type"]=="centromere"),c(1,2,3)]
  chr_cent[,1] = paste("chr", chr_cent[,1], sep="")
  
  ###to get pileup_segs_data_forplot in abspos
  indschrabspos_match = match(pileup_segs_data_forplot[,"Chromosome"], chr_add_for_abspos[,1]) 
  abs_pos_start = as.numeric(pileup_segs_data_forplot[,"Start"])+as.numeric(chr_add_for_abspos[indschrabspos_match,2])
  abs_pos_end = as.numeric(pileup_segs_data_forplot[,"End"])+as.numeric(chr_add_for_abspos[indschrabspos_match,2])
  pileup_segs_data_forplot_in_abspos = cbind(abs_pos_start, abs_pos_end, as.numeric(pileup_segs_data_forplot[,"LogRatio"]), abs_pos_end-abs_pos_start)
  colnames(pileup_segs_data_forplot_in_abspos) = c("AbsStart", "AbsEnd", "LogRatio", "Width")
  df_pileup_segs_data_forplot_in_abspos = data.frame(absstart = pileup_segs_data_forplot_in_abspos[,"AbsStart"], 
                                                     absend = pileup_segs_data_forplot_in_abspos[,"AbsEnd"], 
                                                     LogRatio = pileup_segs_data_forplot_in_abspos[,"LogRatio"])
  chr_start_ends = NULL
  for (i in 1:(length(chr_add_for_abspos[,1])-1)){
    chr_start_ends = rbind(chr_start_ends, c(i, as.numeric(chr_add_for_abspos[i,2]), as.numeric(chr_add_for_abspos[i+1,2])-2))
  }
  chr_start_ends = rbind(chr_start_ends, c(24, as.numeric(chr_add_for_abspos[24,2]), as.numeric(chr_add_for_abspos[24,2])+as.numeric(chr_lengths[24,2])-1))
  
  evd_cutoff_data[,1] = apply(as.matrix(evd_cutoff_data[,1]),1, function(x) unlist(strsplit(x,split="chr"))[2]) 
  evd_cutoff_data[which(evd_cutoff_data[,1] == "X"),1] = 23
  evd_cutoff_data[which(evd_cutoff_data[,1] == "Y"),1] = 24
  evd_cutoff_data_num = cbind(as.numeric(evd_cutoff_data[,1]), as.numeric(evd_cutoff_data[,2])) 
  evd_cutoff_data_numReord = evd_cutoff_data_num[match(chr_start_ends[,1], evd_cutoff_data_num[,1]),]
  
  logratio_thresholds_chrs = cbind(chr_start_ends, evd_cutoff_data_numReord[,2])  ### bind EVD_Chr_Cutoffs per chr
  df_logratio_thresholds_chrs = data.frame(absstart = logratio_thresholds_chrs[,2], absend = logratio_thresholds_chrs[,3], width = logratio_thresholds_chrs[,3]-logratio_thresholds_chrs[,2]+1, height = logratio_thresholds_chrs[,4])
  
  inds_odd = seq(from=1, to=23, by=2)
  inds_even = seq(from=2, to=24, by=2)
  chrname_odd = paste("",inds_odd, sep="")
  chrname_odd[12] ="X"
  chrname_even = paste("", inds_even, sep="")
  chrname_even[12] ="Y"
  df_chrs_odd = data.frame(absstart = chr_start_ends[inds_odd,2], absend = chr_start_ends[inds_odd,3], chrname = chrname_odd)
  df_chrs_even = data.frame(absstart = chr_start_ends[inds_even,2], absend = chr_start_ends[inds_even,3], chrname = chrname_even)
  
  p1 <- ggplot() + geom_rect(data=df_chrs_odd, aes(xmin = absstart, xmax = absend, ymax = 2.05, ymin=-0.05), fill="gray95") +
    geom_rect(data=df_chrs_even, aes(xmin = absstart, xmax = absend, ymax = 2.05, ymin=-0.05), fill="white") +
    geom_segment(data=df_pileup_segs_data_forplot_in_abspos, aes(x = absstart, y = LogRatio, xend = absend,  yend = LogRatio), colour="black", lwd = 0.5)+
    geom_segment(data=df_logratio_thresholds_chrs, aes(x = absstart, y = height, xend = absend,  yend=height), colour="black", lwd = 2)+
    geom_hline(aes(yintercept=0), size=0.3, color="gray")+
    theme_bw()+ scale_x_continuous(expand=c(0,0), name="")+scale_y_continuous(expand=c(0,0), breaks = seq(-1, 2, by = 0.25), name = "")+
    geom_rect(data=df_chrs_odd, aes(xmin = absstart, xmax = absend, ymax = 2.125, ymin=2.005), fill="white") +
    geom_text(data=df_chrs_odd, aes(x=absstart+(absend-absstart)/2, y=2.0625, label=chrname), size=8,  angle=90, color="black")+
    geom_rect(data=df_chrs_even, aes(xmin = absstart, xmax = absend, ymax = 2.125, ymin=2.005), fill="black") +
    geom_text(data=df_chrs_even, aes(x=absstart+(absend-absstart)/2, y=2.0625, label=chrname), size=8, angle=90, color="white")+
    theme(panel.margin = unit(0, "mm"), panel.grid = element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=16) )
    pdf(file=paste(z,"COMBINED_PLOT.pdf",sep="/"),width=19,height=8)  
    print(p1)
    dev.off()
  
}
  
if(F)
  { 
    write.table(pileup_windows_overlap_segments_Chr, file=paste(EVD_calc_path, "pileup_windows_AllChrs.txt", sep="/"), sep="\t", row.names=F, col.names=F, eol = "\n", quote=F)
    write.table(pileup_windows_overlap_counts_Chr, file=paste(EVD_calc_path, "pileup_windows_overlap_counts_AllChrs.txt", sep="/"), sep="\t", row.names=F, col.names=F, eol = "\n", quote=F)
    write.table(fraction_nonzero_coverage_AllChrs, file=paste(EVD_calc_path, "fraction_nonzero_coverage_AllChrs.txt", sep="/"), sep="\t", row.names=F, col.names=F, eol = "\n", quote=F)
    write.table(entropy_coverage_AllChrs, file=paste(EVD_calc_path, "entropy_coverage_AllChrs.txt", sep="/"), sep="\t", row.names=F, col.names=F, eol = "\n", quote=F)
  }
}

##############################################

TumEVD_cal <- function()
{
  setwd(anaTempLRNN_EVD_Cutoff)
  tumor_norm_path = anaTempLRTN_TumNor_SegMeans_CDSFilt
  tumor_cna_path = anaTempLRTN_TumNor_SegMeans_CDSFilt
  min_seg = 1
  pval_sig = 1
  
  
  for (chrind in 1:(length(chrs))) {
    chr_out_t = NULL
    inds_large_segs_t = NULL
    evd_vals_max = NULL
    evd_vals_min = NULL
    chr_out_t_sizeseg_sig_evdpval = NULL
    
    file_1 <- paste(anaTempLRNN_Pos_EVD_Cutoff, paste("EVDParams_", chrs[chrind], ".RData", sep=""),sep='/')
    chk_1 <- file.exists(file_1)
    if(chk_1)
    {
      load(file=paste(anaTempLRNN_Pos_EVD_Cutoff, paste("EVDParams_", chrs[chrind], ".RData", sep=""),sep='/'))
      evd_vals_abs_max_POS = evd_vals_abs_max
      rm(evd_vals_abs_max)
    }
    
    file_2 <- paste(anaTempLRNN_Neg_EVD_Cutoff, paste("EVDParams_", chrs[chrind], ".RData", sep=""),sep='/')
    chk_2 <- file.exists(file_2)
    if(chk_2)
    {
    load(file=paste(anaTempLRNN_Neg_EVD_Cutoff, paste("EVDParams_", chrs[chrind], ".RData", sep=""),sep='/'))
    evd_vals_abs_max_NEG = evd_vals_abs_max
    rm(evd_vals_abs_max)
    }
    
    
    file_3 <- paste(tumor_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/')
    read_test2 <- file.info(file_3)$size
    if(read_test2!=0)
    {
      chr_out_t = read.table(file_3, header=F, sep="\t")
      inds_large_segs_t = which(chr_out_t[,5] >=min_seg)
      chr_out_t_sizeseg = chr_out_t[inds_large_segs_t, ]
      
      
      #### ESTIMATE PVAL OF SEGMENT USING NORMAL-NORMAL PAIRS BASED EVD #########
      chr_out_t_sizeseg_evd_pval = cbind(chr_out_t_sizeseg, matrix(data=NA, ncol = 1, nrow = length(chr_out_t_sizeseg[,1])))
      for (j in 1:length(chr_out_t_sizeseg_evd_pval[,1])){
        seg_val_j = chr_out_t_sizeseg_evd_pval[j,6]
        if (seg_val_j > 0) chr_out_t_sizeseg_evd_pval[j,7] = pgev(abs(seg_val_j), xi=evd_vals_abs_max_POS@fit$par.ests["xi"], mu = evd_vals_abs_max_POS@fit$par.ests["mu"], beta = evd_vals_abs_max_POS@fit$par.ests["beta"], lower.tail=F)[1]
        if (seg_val_j <= 0) chr_out_t_sizeseg_evd_pval[j,7] = pgev(abs(seg_val_j), xi=evd_vals_abs_max_NEG@fit$par.ests["xi"], mu = evd_vals_abs_max_NEG@fit$par.ests["mu"], beta = evd_vals_abs_max_NEG@fit$par.ests["beta"], lower.tail=F)[1] 
        #abs_seg_val_j = abs(chr_out_t_sizeseg_evd_pval[j,6])
        #chr_out_t_sizeseg_evd_pval[j,7] = pgev(abs_seg_val_j, xi=evd_vals_abs_max@fit$par.ests["xi"], mu = evd_vals_abs_max@fit$par.ests["mu"], beta = evd_vals_abs_max@fit$par.ests["beta"], lower.tail=F)[1] 
      }
      inds_sig_pval = which(chr_out_t_sizeseg_evd_pval[,7] <=pval_sig)
      chr_out_t_sizeseg_sig_evdpval = chr_out_t_sizeseg_evd_pval[inds_sig_pval,]
      colnames(chr_out_t_sizeseg_sig_evdpval) = c("Sample_ID", "Chr", "StartPos", "EndPos", "NumMarkers", "SegVal", "ENVE_Pvalue")
      write.table(chr_out_t_sizeseg_sig_evdpval, file=paste(anaEVDPVal_AS,paste(chrs[chrind], "_MinSeg_", as.character(min_seg), "_EVDPVal_AllSamps.txt", sep=""),sep='/'), sep="\t", row.names=F, col.names=T, quote=F, eol = "\n")
    }
    
    
    
    
    

  }
  
  
  #####combine all chrs CNA into single sheet######### 
  All_Chrs_All_Samps_CNA = NULL
  for (chrind in 1:(length(chrs))) {
    fil_in=paste(anaEVDPVal_AS,paste(chrs[chrind], "_MinSeg_", as.character(min_seg), "_EVDPVal_AllSamps.txt", sep=""),sep='/')
    if(file.exists(fil_in))
    {
      read_test3 <- file.info(fil_in)$size
      if(read_test3!=0)
      {
        tmp = read.table(file=fil_in, sep="\t", header=T)
        All_Chrs_All_Samps_CNA = rbind(All_Chrs_All_Samps_CNA, as.matrix(tmp))
      }
    }
  }
  write.table(All_Chrs_All_Samps_CNA, file=paste(anaRes,paste("AllChrs", "_MinSeg_", as.character(min_seg), "_EVDPVal_AllSamps.txt", sep=""),sep='/'), sep="\t", row.names=F, col.names=T, quote=F, eol = "\n")
  
 
  ##############~~~~~~~~~~~~~~~CREATION OF GISTIC FILES~~~~~~~~~~~~~~~~~~#############################
  enve_segs = All_Chrs_All_Samps_CNA
  enve_segs = as.matrix(enve_segs)
  enve_segs_GISTIC = enve_segs[,c("Sample_ID", "Chr", "StartPos", "EndPos", "NumMarkers", "SegVal")]
  for (i in 1:length(enve_segs[,1])){
    if (as.numeric(enve_segs[i,"ENVE_Pvalue"]) > pval_sig | as.numeric(enve_segs[i,"NumMarkers"]) < num_probes)  enve_segs_GISTIC[i,"SegVal"] = "0"
  }
  colnames(enve_segs_GISTIC) = c("Sample","Chromosome",  "Start",	"End",	"Num_Probes",	"Segment_Mean")
  enve_segs_GISTIC_out = enve_segs_GISTIC
  if (pmatch("chr", enve_segs_GISTIC[1,"Chromosome"])==1){
    segs_chr = apply(as.matrix(enve_segs_GISTIC[,"Chromosome"]), 1, function(x) unlist(strsplit(x, split="chr"))[2])
    enve_segs_GISTIC_out[,"Chromosome"] = segs_chr
  }
  write.table(enve_segs_GISTIC_out, paste(anaRes,"Gistic_ENVE_file.txt",sep="/"), sep = "\t", row.names=F, col.names=T, eol = "\n", quote=F)
  
  ### create the arraylist file for GISTIC
  arraylist = as.matrix(unique(enve_segs_GISTIC_out[,"Sample"]))
  colnames(arraylist) = c("Array")
  write.table(arraylist, paste(anaRes,"arraylistfile.txt",sep='/'), sep = "\t", row.names=F, col.names=T, eol = "\n", quote=F)
  
  
  
}

##############################################


