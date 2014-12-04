##Function to download required packages#######
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

##############################################

##############################################

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
##############################################
dir_create <- function()
{
  ana <<- paste(enveHome,"analysis",sep='/')
  #####Analysis files#########
  anaPath <<- paste(ana,(paste("Analysis",toString(strftime(Sys.time(),format="%H_%M_%d_%m_%Y")),sep='_')),sep ="/")
  #anaInp <<- paste(anaPath,"Input",sep='/')
  #anaInpBamsNorm <<-paste(anaInp,"Normal_BAMs",sep='/')
  #anaInpBamsTum <<-paste(anaInp,"Tumor_BAMs",sep='/')
  anaTemp <<-paste(anaPath,"temp",sep='/')
  anaTempVScan <<-paste(anaTemp,"VarScan_CNA",sep='/')
  anaTempVScanNN <<-paste(anaTempVScan,"NormalNormal",sep='/')
  anaRes <<- paste(anaPath,"Results",sep='/')
  anaEVDPVal_AS <<- paste(anaTemp,"EVDPVal_AllSamps",sep='/')
  
  
  
  anaTempVScanNNres <<-paste(anaTempVScanNN,"VarScan_Results",sep='/')
  anaTempVScanNN_CBS_GC_crtd <<-paste(anaTempVScanNN,"CBS_Segments_GC_Corrected_CDSFilt",sep='/')
  anaTempVScanNN_NorNor_SegMeans_CDSFilt <<-paste(anaTempVScanNN,"chr_ALLNormNorm",sep='/')
  anaTempVScanNN_OC_CDSFilt <<- paste(anaTempVScanNN,"VScanCC_OC_FiltCDS",sep='/')
  anaTempVScanNN_Tiff_output<<- paste(anaTempVScanNN,"TIFF_output",sep='/')
  anaTempVScanNN_EVD_Cutoff <<-paste(anaTempVScanNN,"EVD_cutoff",sep='/')
  
  
  
  anaTempVScanTN <<-paste(anaTempVScan,"TumorNormal",sep='/')
  anaTempVScanTNres <<-paste(anaTempVScanTN,"VarScan_Results",sep='/')
  anaTempVScanTN_CBS_GC_crtd <<-paste(anaTempVScanTN,"CBS_Segments_GC_Corrected_CDSFilt",sep='/')
  anaTempVScanTN_TumNor_SegMeans_CDSFilt <<-paste(anaTempVScanTN,"Norm_vs_Tum_SegMeans_CDSFilt",sep='/')
  anaTempVScanTN_OC_CDSFilt <<- paste(anaTempVScanTN,"VScanCC_OC_FiltCDS",sep='/')
  
  
  #anaTempVScanTN_allchr_CDSFilt <<- paste(anaTempVScanNN,"chr_allTNSamps_SegMeans_CDSFilt",sep='/')
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
            anaTempVScan,
            anaTempVScanNN,
            anaTempVScanNNres,
            anaTempVScanNN_CBS_GC_crtd,
            anaTempVScanNN_NorNor_SegMeans_CDSFilt,
            anaTempVScanNN_OC_CDSFilt,
            anaTempVScanNN_Tiff_output,
            anaTempVScanNN_EVD_Cutoff,
            anaTempVScanTN,
            anaTempVScanTNres,
            anaTempVScanTN_CBS_GC_crtd,
            anaTempVScanTN_TumNor_SegMeans_CDSFilt,
            anaTempVScanTN_OC_CDSFilt,
            anaRes,
            anaEVDPVal_AS
  )
  
  
  
  
  
  
  dir2 <- as.data.frame(cbind(c('enveHome',
                                'scriptsPath',
                                'supFiles',
                                'ana',
                                'anaPath',
                                'anaTemp',
                                'anaTempVScan',
                                'anaTempVScanNN',
                                'anaTempVScanNNres',
                                'anaTempVScanNN_CBS_GC_crtd',
                                'anaTempVScanNN_NorNor_SegMeans_CDSFilt',
                                'anaTempVScanNN_OC_CDSFilt',
                                'anaTempVScanNN_Tiff_output',
                                'anaTempVScanNN_EVD_Cutoff',
                                'anaTempVScanTN',
                                'anaTempVScanTNres',
                                'anaTempVScanTN_CBS_GC_crtd',
                                'anaTempVScanTN_TumNor_SegMeans_CDSFilt',
                                'anaTempVScanTN_OC_CDSFilt',
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

##############################################
chrs = as.matrix(c('chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
                   'chr18','chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6',
                   'chr7','chr8','chr9','chrX','chrY'))



chr_proc <- function()
{
  chr_lengths = read.table("/Projects/ENVE/support_files/chr_lengths_hg19.txt", sep="\t")
  chr_lengths = as.matrix(chr_lengths)
  return(chr_lengths)
}

##############################################

filtCDS <- function(x,y)
{
  setwd(x)
  #setwd(anaTempVScanNNres)
  called.files = list.files(pattern="out.called")
  #inds_gc = grep("called.gc", called.files )
  #inds_homdel = grep("called.homdel", called.files )
  #cn_called_files = called.files[c(-inds_gc, -inds_homdel)]
  #inds_FiltCDS = grep("copycaller.FiltCDS", called.files )
  #cn_called_files = called.files[c(-inds_gc, -inds_homdel, -inds_FiltCDS)]
  #y <- anaTempVScanNN_OC_CDSFilt
  cn_called_files = called.files
  
  ##### FILTER WINDOWS TO THOSE THAT FALL WITHIN CODING REGIONS ######
  if(Whole_Exome)
  {
  for (i in 1:length(cn_called_files)) {
    sys_cmd = paste("awk 'NR>1' ", cn_called_files[i], " | sed 's/ //g' | ",IntersectBED_path," -a stdin -b ",paste(supFiles,"RefSeqGTF/RefSeq_hg19_Feb2009.gtf.CDS.nospace.bed.txt",sep='/')," -f 0.6 -u > ",paste(y,paste(cn_called_files[i],".FiltCDS",sep=""),sep='/'), sep="")
    system(sys_cmd)
  }
  }else{
    for(i in 1:length(cn_called_files))
    {
      file.copy(paste(x,cn_called_files[i],sep='/'),y, overwrite = recursive, recursive = TRUE, copy.mode = TRUE)
    }
    setwd(y)
    for(i in 1:length(cn_called_files))
    {
      file.rename(paste(cn_called_files[i],'out.called',sep='.'), paste(cn_called_files[i],'FiltCDS',sep='.'))
    }
  }
  return(cn_called_files)
  ##### END OF FILTER WINDOWS TO THOSE THAT FALL WITHIN CODING REGIONS ######
}

##############################################

CBS_seg_samp <- function(x,y)
{
  setwd(x)
  ##### CBS Segmentation Per SAMPLE ######
  
  for (i in 1:length(cn_called_files)) {
    print(i)
    filcna = paste(cn_called_files[i], ".FiltCDS", sep="")
    sampcna = paste(unlist(strsplit(filcna, split="\\_|\\."))[c(2,4,6)], collapse="_")
    cn <- read.table(filcna,header=F)
    cn = as.matrix(cn)
    cn = str_replace_all(cn, " ", "")
    CNA.object <-CNA(genomdat = as.numeric(cn[,7]), chrom = cn[,1], maploc = as.numeric(cn[,2]), data.type = 'logratio', sampleid=sampcna)
    CNA.smoothed <- smooth.CNA(CNA.object)
    segs <- segment(CNA.smoothed, verbose=0, min.width=2, undo.splits = "sdundo", undo.SD=3)
    segs2 = segs$output
    write.table(segs2[,2:6], file=paste(y,paste(sampcna, "CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t")
  }
}

##############################################

Com_samp_perchr <- function(x,y)
{
  setwd(x)
  ##### BEGIN Combine All Samples per chromosome ######
  for (j in 1:length(chrs)) {
    chr_out = NULL
    for (i in 1:length(cn_called_files)) {
      filcna = cn_called_files[i]
      filcna2 <<- unlist(strsplit(filcna, split="\\_|\\.|\\-"))[c(6)]
      if (chrs[j] != "chrX" & chrs[j] != "chrY"){
        sampcna = paste(unlist(strsplit(filcna, split="\\_|\\.|\\-"))[c(2,4,6)], collapse="_")
        #setwd("/Projects/ENVE/temp/VarScan_CNA/NormalNormal")
        segs.all.chr <- read.table(paste(sampcna, "CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),header=F, sep="\t")
        ind_chr = which(segs.all.chr[,1]==chrs[j])
        samp_chr_out =  cbind(matrix(data=sampcna, nrow = length(ind_chr), ncol=1), segs.all.chr[ind_chr,])
        chr_out = rbind(chr_out, samp_chr_out)
      }else if(chrs[j] =="chrX"){
        if(filcna2 == "MaleMale" | filcna2 == "FemaleFemale"){
          sampcna = paste(unlist(strsplit(filcna, split="\\_|\\.|\\-"))[c(2,4,6)], collapse="_")
          segs.all.chr <- read.table(paste(sampcna, "CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),header=F, sep="\t")
          ind_chr = which(segs.all.chr[,1]==chrs[j])
          samp_chr_out =  cbind(matrix(data=sampcna, nrow = length(ind_chr), ncol=1), segs.all.chr[ind_chr,])
          chr_out = rbind(chr_out, samp_chr_out)
          }
      }else if(chrs[j] == "chrY"){
        if(filcna2 == "MaleMale"){
          sampcna = paste(unlist(strsplit(filcna, split="\\_|\\.|\\-"))[c(2,4,6)], collapse="_")
          segs.all.chr <- read.table(paste(sampcna, "CBS_Segments_GC_Corrected_CDSFilt.logRatio.txt", sep="_"),header=F, sep="\t")
          ind_chr = which(segs.all.chr[,1]==chrs[j])
          samp_chr_out =  cbind(matrix(data=sampcna, nrow = length(ind_chr), ncol=1), segs.all.chr[ind_chr,])
          chr_out = rbind(chr_out, samp_chr_out)
        }
      }
    }
    write.table(chr_out, file=paste(y,paste(chrs[j],"AllNormNormSamps_SegMeans_CDSFilt.txt",sep="_"),sep='/'), row.names=F, col.names=F, quote=F, sep="\t", eol="\n")
  }
  ##### END Combine All Samples per chromosome ######
}

##############################################

Nor_EVD_calc <- function()
{
  setwd(anaTempVScanNN_NorNor_SegMeans_CDSFilt)
  norm_norm_path <- anaTempVScanNN_NorNor_SegMeans_CDSFilt
  tiff_output_path = anaTempVScanNN_Tiff_output
  EVD_calc_path = anaTempVScanNN_EVD_Cutoff
 
  EVD_cutoff_Chr = NULL

  for (chrind in 1:(length(chrs)-2)) {
    #chrind = 1
    chr_out_n = NULL
    inds_large_segs_n = NULL
    #chr_out_n = read.table(paste(norm_norm_path,chrs[chrind],"_SegMeans_CDSFilt.txt",sep=""), header=F, sep="\t")
    if (chrs[chrind] != "chrX" & chrs[chrind] != "chrY")  {
      chr_out_n = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
    } else if (chrs[chrind] == "chrX"){
      chr_out_nF = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
      chr_out_nM = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
      chr_out_n = rbind(chr_out_nF, chr_out_nM)
    } else if (chrs[chrind] == "chrY") {
      chr_out_n = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
      chr_out_n = read.table(paste(norm_norm_path,paste(chrs[chrind],"_AllNormNormSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
    }
  
    inds_large_segs_n = which(chr_out_n[,5] >=min_seg)
    chr_out_n_sizeseg = chr_out_n[inds_large_segs_n, ]
    chr_len = as.numeric(chr_lengths[which(chr_lengths[,1] == chrs[chrind]),2])
    #### create non-verlapping windows of 5kb length to cover the entire chromosome ####
    win_size = 5000
    epsilon = 1e-16
    num_wins = chr_len%/%win_size
    chr_wind_ends = seq(from=0, to = chr_len, by = win_size)
    chr_wind_ends = c(chr_wind_ends[-1],chr_len) 
    chr_wind_begins = seq(from=1, to = chr_len, by = win_size)
    chr_winds = cbind(chr_wind_begins, chr_wind_ends)
    chr_winds_IR = IRanges(chr_wind_begins, chr_wind_ends)
    #find all segments whose absolute mean segval falls within segval window and plot on chr ####
    ##########chk######
     tiff(filename = paste(tiff_output_path,paste("MaxSegvalEstForEVD_", chrs[chrind],"tiff", sep="."),sep='/'), width = 800, height = 800, pointsize = 12, compression = c("none"),bg = "white")
     par(mfrow=c(1,3))
     plot(0,0, xlim=c(0,chr_len), ylim = c(0,segval_windows[length(segval_windows)]), lwd=0.001, xlab = "", ylab = "", xaxt='n', yaxt='n', main = chrs[chrind], col = "grey90", mai = c(0.5, 0,0.5,0), mar = c(0,1,1,0), oma = c(0,1,1,0))
    #########chkend#####
    chr_winds_overlap_entropy = matrix(data=NA, nrow=1, ncol = length(segval_windows))
    chr_winds_overlap_freq_matrix = matrix(data=NA, nrow=length(segval_windows), ncol = dim(chr_winds)[1])
    chr_winds_overlap_cnt_matrix = matrix(data=NA, nrow=length(segval_windows), ncol = dim(chr_winds)[1])
  
  
    for (winind in 1:length(segval_windows)){
      #segs_in_wind = chr_out_n_sizeseg[which(abs(chr_out_n_sizeseg[,6])>=segval_windows[winind, 1]  & abs(chr_out_n_sizeseg[,6])<segval_windows[winind, 2]),]
      segs_in_wind = chr_out_n_sizeseg[which(abs(chr_out_n_sizeseg[,6])>=segval_windows[winind]),]
      segs_in_wind_IR = IRanges(segs_in_wind[,3], segs_in_wind[,4])
      chr_winds_overlap_cnt = countOverlaps(chr_winds_IR,segs_in_wind_IR)
      chr_winds_overlap_cnt_matrix[winind,] = chr_winds_overlap_cnt
      chr_winds_overlap_freq = chr_winds_overlap_cnt/sum(chr_winds_overlap_cnt)
      chr_winds_overlap_freq_matrix[winind,] = chr_winds_overlap_freq
      chr_winds_overlap_entropy[winind] = -1*sum(chr_winds_overlap_freq*log(chr_winds_overlap_freq), na.rm=T)
      #print(chrind)
      #print(segs_in_wind)
      if (dim(segs_in_wind)[1] >0){
        segx0 = segs_in_wind[,3]
        segx1 = segs_in_wind[,4]
        segy0 = segval_windows[winind]
        segy1 = segy0
        segments(segx0, segy0, x1=segx1, y1=segy1)         
      }
    }
  
    #identify minimum segval means that have at least 10 chromosomal windows with non-zero frequency 
    #this gives us the maximum mean segval for which we want to analyze the distribution of reads
    nonNA_segval_windows_index = which(apply(chr_winds_overlap_freq_matrix,1,function(x) length(x)-length(which(is.na(x))))>=10)
    nonNA_chr_winds_overlap_freq_matrix = chr_winds_overlap_freq_matrix[nonNA_segval_windows_index,]
    chr_winds_overlap_entropy_nonNA = chr_winds_overlap_entropy[nonNA_segval_windows_index]
    nonNA_segval_windows = segval_windows[nonNA_segval_windows_index]
    
    
    #identify fraction of chromosome that has a non-zero frequency of reads aligning to it for every meansegval cutoff 
    frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix = apply(nonNA_chr_winds_overlap_freq_matrix,1, function(x)  length(which(x!=0))/length(x))
    
    ### pick the seg-val cutoff where the genomic coverage drop is maximum  -- WE FOUND THIS BEFORE 
    winind_sig = which.min(diff(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix))+1
    diff_frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix = NULL
    diff_frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix = c(NA, -diff(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix))
  
    ### if chr fraction at maximum segval is greater than 25%, use all segvals for EVD ###
    if (frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix[length(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix)] >= 0.25){
      droppoint_segval_value = max(nonNA_segval_windows) + 0.005
    } else {
      droppoint_segval_index = which.min(diff(frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix))+1
      
      ### begin estimate Silhouette index for cutoff point to ensure that segment distribution above cutoff is dissimilar to segment distribution below ###
      
      
      #test_sil1 <- try(Get_Sil_Index(nonNA_chr_winds_overlap_freq_matrix, droppoint_segval_index))
      #print(test_sil1)           
      #test_sil2 <- try(Get_Sil_Index(nonNA_chr_winds_overlap_freq_matrix, droppoint_segval_index-1))
      #print(test_sil2)
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
    
    #   ### BEGIN randomize all the windows for every row above the midpoint ####
    #   RND_nonNA_chr_winds_overlap_freq_matrix = nonNA_chr_winds_overlap_freq_matrix
    #   for (i in (droppoint_segval_index+1):dim(nonNA_chr_winds_overlap_freq_matrix)[1]){
    #     x = RND_nonNA_chr_winds_overlap_freq_matrix[i,]
    #     RND_nonNA_chr_winds_overlap_freq_matrix[i,] = x[shuffle(length(x))]
    #   }
    #   ### END randomize all the windows for every row above the midpoint ####
    #   si_RND = Get_Sil_Index(RND_nonNA_chr_winds_overlap_freq_matrix, droppoint_segval_index)
    
    
    
    ####chk################
    abline(h=droppoint_segval_value, lty="solid", col="grey10", lwd = 9)
    #plot(nonNA_segval_windows, frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix, family = "sans", type="o", lwd = 2, cex.main = 1.5, cex.axis = 0.9, pch = 16, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), main = chrs[chrind])
    plot(nonNA_segval_windows, diff_frac_nonzero_freq_nonNA_chr_winds_overlap_freq_matrix, type="o", lwd = 2, cex.main = 1, cex.axis = 1.5, pch = 16, xlab = "", ylab = "", xlim = c(0,1), main = chrs[chrind])
    abline(v=droppoint_segval_value, lty="dotted", col="grey10", lwd = 2)
    plot(nonNA_segval_windows, chr_winds_overlap_entropy_nonNA, type="o", family = "sans", lwd = 2, cex.main = 1.5, cex.axis = 0.9, pch = 16, xlab = "", ylab = "", xlim = c(0,1), main = chrs[chrind])
    abline(v=droppoint_segval_value, lty="dotted", col="grey10", lwd = 2)
    #######################
    
    dev.off()
    EVD_cutoff_Chr = rbind(EVD_cutoff_Chr, c(chrs[chrind], droppoint_segval_value))
  }
  
  #dev.off()
  colnames(EVD_cutoff_Chr) = c("Chromosome", "SegValCutoffForEVD")
  write.table(EVD_cutoff_Chr, file=paste(EVD_calc_path,paste("EVD_Chr_Cutoffs","txt", sep="."),sep='/'), sep="\t", quote=F, row.names=F, col.names=T, eol="\n")
  
  ###################################################################################################################################
  
  
  EVD_calc_path = anaTempVScanNN_EVD_Cutoff
  EVD_cutoff_Chr = read.table(paste(EVD_calc_path,paste("EVD_Chr_Cutoffs","txt", sep="."),sep='/'), header=T)
  EVD_cutoff_Chr = as.matrix(EVD_cutoff_Chr)
  
  
  EVD_params_per_chr = NULL
  
  for (chrind in 1:(length(chrs))) {
    #chrind = 6
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
    training_vals_max = NULL
    training_vals_min = NULL
    for (i in 1:length(samp_names)){
      max_val = max(chr_out_n_sizeseg_for_EVD[which(chr_out_n_sizeseg_for_EVD[,1]==samp_names[i]),6])
      min_val = min(chr_out_n_sizeseg_for_EVD[which(chr_out_n_sizeseg_for_EVD[,1]==samp_names[i]),6])
      if (max_val > 0) training_vals_max = c(training_vals_max, max_val) 
      #training_vals_max = c(training_vals_max, max_val)
      if (min_val < 0) training_vals_min = c(training_vals_min, min_val) 
      #training_vals_min = c(training_vals_min, min_val) 
    }
    #print("true")
    evd_vals_max=gevFit(x = training_vals_max, block = 1, type = "pwm")
    #print("true")
    evd_vals_min=gevFit(x = training_vals_min, block = 1, type = "pwm")
    #print("true")
    save(evd_vals_max,evd_vals_min, file=paste(EVD_calc_path,paste("EVDParams_", chrs[chrind], ".RData", sep=""),sep='/'))
    EVD_params_per_chr = rbind(EVD_params_per_chr, c(chrs[chrind], evd_vals_min@fit$par.ests, evd_vals_max@fit$par.ests))
  }
  colnames(EVD_params_per_chr) = c("Chromosome", "Min_xi", "Min_mu", "Min_beta", "Max_xi", "Max_mu", "Max_beta")
  write.table(EVD_params_per_chr, file=paste(EVD_calc_path, "EVD_Params_AllChrs.txt", sep="/"), sep="\t", row.names=F, col.names=T, eol = "\n", quote=F)
  ##################################################################################################################################################
  if(F)
  {
  }
}

##############################################

TumEVD_cal <- function()
{
  setwd(anaTempVScanNN_EVD_Cutoff)
  tumor_norm_path = anaTempVScanTN_TumNor_SegMeans_CDSFilt
  #norm_norm_path = "/Projects/ENVE/temp/VarScan_CNA/NormalNormal/"
  
  EVD_calc_path = anaTempVScanNN_EVD_Cutoff
  tumor_cna_path = anaTempVScanTN_TumNor_SegMeans_CDSFilt
  min_seg = 1
  pval_sig = 1
  
  
  for (chrind in 1:(length(chrs))) {
    chr_out_t = NULL
    inds_large_segs_t = NULL
    evd_vals_max = NULL
    evd_vals_min = NULL
    chr_out_t_sizeseg_sig_evdpval = NULL
    load(file=paste(EVD_calc_path, paste("EVDParams_", chrs[chrind], ".RData", sep=""),sep='/'))
    chr_out_t = read.table(paste(tumor_norm_path,paste(chrs[chrind],"_AllSamps_SegMeans_CDSFilt.txt",sep=""),sep='/'), header=F, sep="\t")
    
    inds_large_segs_t = which(chr_out_t[,5] >=min_seg)
    chr_out_t_sizeseg = chr_out_t[inds_large_segs_t, ]
    
    
    
    
    
    #### ESTIMATE PVAL OF SEGMENT USING NORMAL-NORMAL PAIRS BASED EVD #########
    chr_out_t_sizeseg_evd_pval = cbind(chr_out_t_sizeseg, matrix(data=NA, ncol = 1, nrow = length(chr_out_t_sizeseg[,1])))
    for (j in 1:length(chr_out_t_sizeseg_evd_pval[,1])){
      seg_val_j = chr_out_t_sizeseg_evd_pval[j,6]
      if (seg_val_j > 0) chr_out_t_sizeseg_evd_pval[j,7] = pgev(seg_val_j, xi=evd_vals_max@fit$par.ests["xi"], mu = evd_vals_max@fit$par.ests["mu"], beta = evd_vals_max@fit$par.ests["beta"], lower.tail=F)[1]
      if (seg_val_j <= 0) chr_out_t_sizeseg_evd_pval[j,7] = pgev(seg_val_j, xi=evd_vals_min@fit$par.ests["xi"], mu = evd_vals_min@fit$par.ests["mu"], beta = evd_vals_min@fit$par.ests["beta"], lower.tail=T)[1]  
    }
    inds_sig_pval = which(chr_out_t_sizeseg_evd_pval[,7] <=pval_sig)
    chr_out_t_sizeseg_sig_evdpval = chr_out_t_sizeseg_evd_pval[inds_sig_pval,]
    colnames(chr_out_t_sizeseg_sig_evdpval) = c("Sample_ID", "Chr", "StartPos", "EndPos", "NumMarkers", "SegVal", "ENVE_Pvalue")
    write.table(chr_out_t_sizeseg_sig_evdpval, file=paste(anaEVDPVal_AS,paste(chrs[chrind], "_MinSeg_", as.character(min_seg), "_EVDPVal_AllSamps.txt", sep=""),sep='/'), sep="\t", row.names=F, col.names=T, quote=F, eol = "\n")
  }
  
  
  #####combine all chrs CNA into single sheet######### 
  All_Chrs_All_Samps_CNA = NULL
  for (chrind in 1:(length(chrs))) {
    fil_in=paste(anaEVDPVal_AS,paste(chrs[chrind], "_MinSeg_", as.character(min_seg), "_EVDPVal_AllSamps.txt", sep=""),sep='/')
    tmp = read.table(file=fil_in, sep="\t", header=T)
    All_Chrs_All_Samps_CNA = rbind(All_Chrs_All_Samps_CNA, as.matrix(tmp))
  }
  write.table(All_Chrs_All_Samps_CNA, file=paste(anaRes,"Results.txt",sep='/'), sep="\t", row.names=F, col.names=T, quote=F, eol = "\n")
}

##############################################


