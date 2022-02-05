############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "dataset_TableS11.txt"
###
### Author:   Niccolo` Bassetti - niccolo.bassetti@wur.nl 
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
x=c("FSA","rcompanion")

lapply(x, install.packages, character.only = TRUE)

# Set up working directory
getwd()                 # see current working directory
setwd("path_to_file")	  # set working directory where input datasets are located. 

setwd("C:/Users/nicco/OneDrive - Wageningen University & Research/5_Articles/Chapter_4_QTL_Brassica_rapa/BMC - revision/repeat_analysis")

# Load datasets
mydata <- read.table(file = "dataset_TableS11.txt",sep = "\t", header = TRUE) 
names(mydata)
str(mydata)


# Examine dataset - remove NAs
sapply(mydata, function(x)(sum(is.na(x))))		# search NAs


# Transform factor written as "integer"
sapply(FUN = class, mydata)
mydata$Rep <- as.factor(mydata$Rep)
mydata$Plant_ID <- as.factor(mydata$Plant_ID)
mydata$Genotype <- as.factor(mydata$Genotype)
mydata$Crop_type <- as.factor(mydata$Crop_type)
mydata$Wash_HR_MAX_plant <- as.factor(mydata$Wash_HR_MAX_plant)

sapply(FUN = class, mydata)
str(mydata)



############################
### Statistical analysis ###
###


### Select phenotypic trait to be analysed
names(mydata)

trait = "Wash_HR_MAX_plant" # select phenotypic trait to analyse!!!

mydata$Response = mydata[,which(names(mydata)==trait)]
head(mydata)


### Non-parametric test
#   final results reported in manuscript
kruskal.test(mydata$Response~mydata$Genotype)


### Pairwise mulriple comparison: Dunn`s test
#   (reported in Supplementary Table S2)

library("FSA")

mydata$Response <- as.integer(mydata$Response)
mydata$Genotype = as.factor(gsub("-", "_", mydata$Genotype)) # replace dash character from Genotype names. Otherwise there will be an error later. 

PT = dunnTest(mydata$Response~mydata$Genotype,
              data=mydata,
              method="bh")
PT
PT = PT$res # prepare dataset to extract P-values from multiple comparison

library("rcompanion")
post.hoc = cldList(comparison = PT$Comparison, # gives weird results, not correspondin with Wilcoxon
                   p.value    = PT$P.adj,
                   threshold  = 0.05) # require "rcompanion"
post.hoc = as.data.frame(post.hoc, row.names = TRUE)

# export results
write.table(post.hoc, sep = "\t", file = "Germplasm_screen_posthoc_P.adj.txt")



################  
### Plot data
###
### Supplementary Figure S1
###
### Data were plotted in Excel 365 using "Stacked Column" chart
