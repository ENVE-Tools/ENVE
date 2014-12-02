ENVE Version 1.0
====

====Introduction====

ENVE is a tool which first models inherent noise in WES data using non-tumor diploid samples, and utilizes the learned model parameters to estimate sCNAs in tumors in an unbiased way. The ENVE methodology, in general, consists of two major modules, which include capturing and modeling the inherent noise (sample- and technical-associated variability) in whole-exome sequencing data using non-tumor diploid normal samples, followed by utilizing the learned model parameters to reliably detect somatic copy-number alterations in tumors.

The key features of ENVE include 
* Empirical quantification inherent noise in WES data using random normal-normal combinations of diploid samples.
* Use of three independent measures to differentiate between segmental LogRatio variations associated with inherent noise versus germline CNVs within these normal-normal comparisons.
* Modeling extreme deviations of chromosome-specific LogRatios associated with inherent noise within the diploid normal-normal comparisons by using a generalized extreme value distribution.
* By use of these learned model parameters, to evaluate somatic copy-number alterations in tumor samples. 


In this guide we provide easy to follow instructions for installing the ENVE module and running the results from prENVE. which usually takes the Bam files and do pair wise comparison to find the log ratios in the form of *.out.called files. 

==============================================================================================================================

====System Requirements====

==Platforms and System Requirements==

The module performance was tested on two platforms

* 2.8 GHz Intel Core i7 processor with 8 GB of RAM iMac machine
* 32-core 2.10 GHz Intel(R) Xeon(R) CPU E5-2450 with 64 GB linux machine

The performance of the script depends on the available memory and number of cores used. on iMac machine 30 samples( 435 Normal Normal pairs) were processed in around 10 hrs. The performance on the multiple core high RAM machine it was less than couple of hours.

Hard Disk Space : The module needs around 30 MB of Hard Disk space

==Software Architecture==

The script was developed on the R platform with version 3.1, and it is assumed that the user also use the same for running the script.

==Software Requirements==

R : Free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms. 
preferred Version : >= 3.1
Bedtools : It is a set of tools which can be used for intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF. 
preferred Version : >= V.2.21.0

The following commands will install bedtools in a local directory on an UNIX or OS X machine. Note that the “” refers to the latest posted version number on

* $ curl http://bedtools.googlecode.com/files/BEDTools..tar.gz > BEDTools.tar.gz
* $ tar -zxvf BEDTools.tar.gz
* $ cd BEDTools-
* $ make

For more information : http://bedtools.readthedocs.org/en/latest/content/installation.html

==============================================================================================================================

====Installation And Running Tool====

Download the ENVE-1.0 repository
In subdirectory scripts one need to edit two files

ENVE_CONF.txt : this file contains a field which points to the intersectBed tool present in the bedtools/bin/intersectBed, and needed to be updated.
ENVE_RUN_CONF.TXT : This file consist of information required for running the batch. Paths to the Normal Normal Paired data (*.out.called) files needed to be provided, The text file also asks user if one is running Whole exome data and If wants to run both Normal-Normal ENVE module and Tumor-Normal module together.
To Run the ENVE Module

* $ cd directory/to/ENVE 
* $ R CMD BATCH ENVE.R

The results file be stored in analysis folder with the analysis name with timestamp. all the required and intermediate files can be found in temp folder.

==============================================================================================================================


====Understanding the Results====

==Important Files==

Tiff Files
location : analysis/temp/NorNor/VarScan_CNA/NormalNormal/Tiff_output
these files shows the Absolute LogRatio Thresholds for each chromosome and the segments which are above the threshold.
EVD_cutoff

location : analysis/temp/NorNor/VarScan_CNA/NormalNormal/EVD_cutoff/
this file includes the Absolute LogRatio Thresholds for each chromosome
Results.txt

location : analysis/Results/
This file contains the final ENVE p values of all the segments after removing the noise by the use of Normal Diploid comparisons.


For further documentation or technical support while installing or using the package, kindly mail to sxs1528@case.edu 

