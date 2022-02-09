############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  Results from "QTL_analysis_figure_3.R"
###
### Author:   Niccolo` Bassetti - niccolo.bassetti@protonmail.com
###
### Principal Investigator:   Nina Fatouros - nina.fatouros@wur.nl



# Clean R workspace
rm(list=ls())
gc()

# set up new R workspace
sessionInfo() #check packages available

.libPaths()   # see Path to R library
library()     # inspect packages available in Path


# Install R packages required for this analysis
update.packages()   # update all packages already installed

# install multiple packages at same time
x=c("qtl", "ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse")

lapply(x, install.packages, character.only = TRUE)


# Set up working directory
getwd()                 # see current working directory
setwd("paste_here_you_path_to_file")	  # set working directory where input datasets are located. 




####       Load results of QTL analysis      ####
#_______________________________________________#


## NOTE: results/objects generated in the following analysis and reported in the manuscript are store in the .RData file below.
#  Those results can be loaded to quickly generate plots and/or compare results from a new analysis. 
load("pieris.celldeath.brapa.QTLmapping.RData")






####       Supplementary figure S6      ####              
#__________________________________________#

library("qtl")

### plot results from Twoscan() - the analysis of two-QTLs models

# For details on the analysis see "QTL_analysis_figure_3.R", line 104-143 

# Plot
plot(out2, upper="fv1"  # fv1 = full model (epistatic additive and interaction effects) 
     ,lower="av1"       # av1 = additive model (without epistatic interaction effects)
     , ylab=c("")
     , xlab=c("")
     , xaxt='n'
     , yaxt='n')
title(#main=paste("A",i, sep=""),
  #, sub="sub-title"
  xlab="Chromosome"
  , ylab="Chromosome"
  , font.main=2, font.lab=2, font.sub=2
  
  , cex.main=4, cex.main=4, cex.main=6)
#axis(1, seq(0,max.cM.scale,20), font=2, cex=6)
#axis(2, seq(0,max.lod.score,2), font=2, cex=6)
dev.off()

png(file=paste("Two_scan_QTLs",".png", sep=""), width=2500,height=2000,res=300) #

dev.off()
