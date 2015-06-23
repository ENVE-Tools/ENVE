##ENVE Version 1.0.0


##### For technical support, send us an email at enve_support@googlegroups.com

==============================================================================================================================

###Introduction
==============================================================================================================================

ENVE is a tool which first models inherent noise in WES data using non-tumor diploid samples, and utilizes the learned model parameters to estimate sCNAs in tumors in an unbiased way. The ENVE methodology consists of two major modules, which include capturing and modeling the inherent noise (sample- and technical-associated variability) in whole-exome sequencing data using non-tumor diploid normal samples, followed by utilizing the learned model parameters to reliably detect somatic copy-number alterations in tumors.

The key features of ENVE include 
* Empirical quantification inherent noise in WES data using random normal-normal combinations of diploid samples.
* Use of three independent measures to differentiate between segmental LogRatio variations associated with inherent noise versus germline CNVs within these normal-normal comparisons.
* Modeling extreme deviations of chromosome-specific LogRatios associated with inherent noise within the diploid normal-normal comparisons by using a generalized extreme value distribution.
* By use of these learned model parameters, to evaluate somatic copy-number alterations in tumor samples. 


In this guide we provide easy to follow instructions for installing the ENVE module and running the results from prENVE, which takes the Bam files and do a pair wise comparison to find the log ratios in the form of *.GC_COR_adj_logratio files. 

==============================================================================================================================

##System Requirements
==============================================================================================================================
###Platforms and System Requirements

###Introduction 

This page describes the various system and software required to run the ENVE V1.0.0

###Platforms and System Requirements

The module performance was tested on two platforms 
  * 2.8 GHz Intel Core i7 processor with 8 GB of RAM iMac machine
  * 32-core  2.10 GHz Intel(R) Xeon(R) CPU E5-2450 with 64 GB linux machine

The performance of the script depends on the available memory and the number of cores used. ENVE was able to process a cohort of 30 matched tumor and normal pairs, in under 10 hours on the iMac machine described above. However, the same cohort was processed by ENVE in under 2 hours on the multiple-core, high RAM Linux machine.

Hard Disk Space : The module needs around 30 MB of Hard Disk space

###Software Architecture

The script was developed on the R platform with version 3.1, and it is assumed that the user also uses the same for running the script. 

###Software Requirements

 * [**R**](http://www.r-project.org/) : Free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms.
<br> preferred Version : >= 3.1 

### preENVE

If using preENVE tool to get GC corrected log ratio, user needs to have two external software packages 

1. [**VarScan**](http://varscan.sourceforge.net/) : (Version 2.3.6 or Up) 
2. [**SAMtools**](http://samtools.sourceforge.net/)
<br><br>

The individual download and Installation guides are available on their respective websites. 

### ENVE

        
 * [**Bedtools**](http://bedtools.readthedocs.org/en/latest/) : It is a set of tools which can be used for intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF.
        <br> preferred Version : >= V.2.21.0

    The following commands will install bedtools in a local directory on an UNIX or OS X machine. 

    * $ curl http://bedtools.googlecode.com/files/BEDTools.v2.21.0.tar.gz > BEDTools.tar.gz
    * $ tar -zxvf BEDTools.tar.gz
    * $ cd BEDTools.v2.21.0
    * $ make

For more information : http://bedtools.readthedocs.org/en/latest/content/installation.html

==============================================================================================================================

* Download the ENVE-1.0.0 repository

## preENVE
* In the subdirectory Run_Conf_Files user needs to edit preENVE_PROJ_Config.txt 
 * This text file have pointers to the required external software and environment packages.   
Information required : 
    * preENVE : preENVE Home(/location of Enve-master/preENVE)
    * samTools : Directory where the SamTools is installed. 
    * VarScan : Directory where Varscan is installed
    * JAVA_HOME : Directory where Java is installed ( /usr/lib/jvm/jre<version>/bin)
    * hg19_karyo : pointer to the reference fasta file.  
  * Debug(TRUE/FALSE) : shows all the table in output
  * BASH (TRUE/FALSE) : for creating Bash scripts for all individual steps( preferable for smaller batches) 
  * HPC (TRUE/FALSE): For future development
  * Com_Scr (TRUE/FALSE) : If you need all scripts combined together in batches of 8 commands( consider use of 8 core machine).
  
   * NormBam : Directory where the BAM files for Normal Samples are saved
   * TumBam : Directory for where BAM files for Tumor Samples are saved
   * samp_info_file : Samp Info file which includes the required information about samples (follow the headers and format of the provided Sample Info Sheet).
   * Number_of_samp : Number of Norm Samples to include for EVD calculation


To Run the prENVE Module 

$ cd _directory/to/preENVE_ <br>
$ RScript preENVE.R [-R] [location to preENVE_PROJ_Config.txt]<br>
for help : <br>
$ RScript preENVE.R [-H]

The preENVE will provide with the scripts which will generate the required input files for ENVE. 

## ENVE
* In subdirectory Run_Conf_Files user needs ENVE_RUN_CONF.txt
     * This file consist of information required for running the ENVE module successfully.
        Options 
        * IntersectBED_path: this field points to the intersectBed tool present in the bedtools/bin/intersectBed, and needed to be updated.
        * NormNorm (TRUE/FALSE): Option for running the Normal module which develops the model for the ENVE. 
        * TumNorm (TRUE/FALSE): Option for running the Module which compares the Tumor Data against the Matched Normal module. 
        * Whole_Exome (TRUE/FALSE) : If you are running Whole exome samples.
        * Input_NormOutcalled_files : pointer to the directory containing the GC corrected adjusted log ratio files for Normal Normal pairs. 
        * Input_TumNorm_Outcalled_files : pointer to the directory containing the GC corrected adjusted log ratio files for Tumor Matched-Normal pairs.
        * Input_NormNorm_Files_Info, Input_TumNorm_Files_Info : Files consist of Sample IDs in pairs and the gender information, which provides the ENVE module with the required information.(See the sample Input_NormNorm_Files_info file for the format).
   

To Run the ENVE Module 

$ cd _directory/to/ENVE_ <br>
$ Rscript ENVE.R [-R] [ENVE_RUN_CONF.txt] <br>
for help : <br>
$ RScript ENVE.R [-H] <br>

The results file will be stored in the analysis folder with the analysis name with a timestamp. All the required and intermediate files can be found in the temp folder.


==============================================================================================================================
###Understanding the Results
==============================================================================================================================
###Important Files

### preENVE
1. combined_Script : It is the set of bash commands needed to get the GC corrected log ratio for a particular pair of samples, this script includes commands to create mpileup using samtools and copynumber information using VarScan.

If user does not need the preENVE module to create GC corrected log ratio, user needs to provide a file which includes chromosome Information with start position and end with the log ratio for the pairs. A sample GC_COR_adj_logratio file is present at ./ENVE_HOME/ENVE/samples/SAMPLE.GC_COR_adj_logratio

2. NorNor_CalledFiles.txt , TumNor_CalledFiles.txt : these file consist of the sample ids and gender information. which is required by ENVE module to understand the pair names and gender information. A sample NorNor_CalledFiles.txt file is present at ./ENVE_HOME/ENVE/samples/NorNor_CalledFiles.txt


###ENVE
  1. Tiff Files
     * location : analysis/temp/NorNor/LogRatio_CNA/NormalNormal/Tiff_output
     * these files shows the Absolute LogRatio Thresholds for each chromosome and the segments which are above the threshold. 
  2. EVD_cutoff 
      * location : analysis/temp/NorNor/LogRatio_CNA/NormalNormal/EVD_cutoff/
      *  this file includes the Absolute LogRatio Thresholds for each chromosome

  3. Results.txt
      * location : analysis/Results/
      * This file contains the final ENVE p values of all the segments after removing the noise by the use of Normal Diploid comparisons.

==============================================================================================================================
