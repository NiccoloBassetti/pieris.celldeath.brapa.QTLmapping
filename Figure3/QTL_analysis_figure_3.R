############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "rqtl_HR_BLUE.csv"
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

# This step can be skipped if interested in repeating all steps of teh analysis

load("pieris.celldeath.brapa.QTLmapping.RData")




####     QTL analysis       ####
#______________________________#


# load data
require("qtl")
require("ggplot2")
data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BLUE.csv", crosstype = "riself")

## NOTE: rqtl_HR_BLUE.csv 
#  this file can be assembled using the BLUE phenotypic data of "dataset_Table_S13.txt" (column "Egg_HR_size_MAX_BLUE")
#  the genotypic data included in "dataset_Table_S14.txt".
#  r/QTL requires a specific format for the input files (including white spaces!). please check r/QTL documentation for details

summary(data)
plot(data)  # plot summary data
plotMap(data, chr=2, show.marker.names = TRUE) # CHANGE #chr to focus on a specific chromosome

# pdf("genetic_map.pdf") # to save genetic map in .pdf
data = jittermap(data)
plotMap(data, chr=2, show.marker.names = TRUE) # CHANGE #chr to focus on a specific chromosome
dev.off()




#### 1) Single QTL analysis ####
#______________________________#

# method: interval mapping with "hk regression"
data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BLUE.csv", crosstype = "riself")
data = jittermap(data)
data_2 = calc.genoprob(data, step=1, error.prob=0.001)    # step = 1 USED

out.hk <- scanone(data_2, method="hk")  # STEP=1 USED
operm.hk <- scanone(data_2, method="hk", n.perm=10000) # permutation test for "HK regression" method
summary(out.hk, threshold=4.2)  # gives LOD peaks above specified threshold
summary(operm.hk, alpha = c(0.1, 0.05, 0.01, 0.001))             # gives a distribution of LOD thresholds for significance testing

summary(out.hk, perms = operm.hk, alpha=0.05) # include permutation results
summary(out.hk, perms = operm.hk, pvalues = TRUE) # include permutation results


  # plotting QTL results for inspection
  par(mfrow=c(4,3))	
  for (i in 1:10) {
    plot(out.hk, chr = i, ylim = c(0,round(max(summary(out.hk)[,3]))+1))  
    abline(h=summary(operm.hk)[1], col=c("red"))
    }
  dev.off()
    
write.csv(out.hk, file = "qtl_HR_MAX_LOD_scores.CSV") # export single QTL analysis to plot in a different environment



 
#### 2) Two-QTL scan (two dimensional) ####
#_________________________________________#

data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BLUE.csv", crosstype = "riself")
data = jittermap(data)
data_2 = calc.genoprob(data, step=1, error.prob=0.001)    # step = 1 USED

out2 = scantwo(data_2, method = "hk", verbose = FALSE)
summary(out2)

opermout2 = scantwo(data_2, method = "hk", n.perm = 1000, verbose = FALSE)
summary(opermout2)  # LOD threshold obtained after 1000 permutations

summary(out2,thresholds = c(7.29, 6, 4.45, 5.3, 3.92)) # input threshold from opermout2
summary(out2, perms=opermout2, alphas=c(0.05,0.05,0.05,0.05,0.05), pvalues=TRUE) # same as above but selcting threshold at specifc error rate. we chose alpha=0.05

## NOTE: Significant interactions can be highligted in the plot below. 


  # plotting QTL results for inspection
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




#### 3) MQM analysis - Multi QTL analysis ####
#____________________________________________#

require("snow")
require("qtl")

rm(list=ls())



### MQMscan() - Unsupervised backward selection ###
#  Note that r/QTL offers also stepwise/forward model selection for MQM analysis.
#  We obtained similar results and here we show only MQMscan()

data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BLUE.csv", crosstype = "riself")
data = jittermap(data)
data_2 = calc.genoprob(data, step=1, error.prob=0.001)    # step = 1 USED

maug_min1 = mqmaugment(data, minprob=1.0)
mqm_min1 <- mqmscan(maug_min1, step.size = 1.0)

real_markers <- mqmextractmarkers(mqm_min1) # remove pseudomarkers for final results
real_markers
  

  # calculation of QTL significance
  results = mqmpermutation(maug_min1, scanfunction=mqmscan  #  change mqmscan() file
                           #, n.cluster=2
                           , n.perm=1000
                           #, batchsize=10
  )
  resultsrqtl = mqmprocesspermutation(results)
  summary(resultsrqtl, alpha = c(0.1, 0.05, 0.01, 0.001))
  
  # final QTL results
  summary(mqm_min1, perms = resultsrqtl, alpha=0.05) # MQMscan() identified two loci on A02 (84 cM) and A06 (64 cM)

  
    
  # plot results (compared to scanone())
  scanone = out.hk
  lod_scanone = operm.hk
  mqm = mqm_min1
  lod_mqm = resultsrqtl
  
  par(mfrow=c(4,3))	
  #for (i in c(2,3,6)) {
  for (i in 1:10) {  
  plot(scanone,mqm, col = c("blue", "red"),
       chr = i,
       ylim = c(0,round(max(summary(mqm)[,3]))+5),
       ylab=c("LOD score"))
  abline(h=summary(lod_mqm)[1], col=c("red")) # ENTER significance for MQM
  abline(h=summary(lod_scanone)[1], col=c("blue"))
  }
  dev.off()


  
### Set markers as cofactors to discover more QTLs

# select peak markers from identified QTLs
summary(mqm_min1, perms = resultsrqtl, alpha=0.05) # cofactors from MQM scan()
#summary(out.hk, perms = operm.hk, alpha=0.05) # cofactors from scanone

chr = c(2,6)    # select chromosome

marker = c()
marker = c(marker, find.marker(data, chr = chr[1], pos = c(86))) # QTL on A02
marker = c(marker, find.marker(data, chr = chr[2], pos = c(64))) # QTL on A06

marker


# First set 1 cofactor
multitoset = find.markerindex(maug_min1, marker[1]) # marker[1] selects QTL on A02
setcofactors <- mqmsetcofactors(maug_min1, cofactors=multitoset)
mqm_co1 <- mqmscan(maug_min1, setcofactors, step.size = 1) # MQM analysis with one QTL as cofactor

  # calculation of LOD threshold significance
  results_co1 = mqmpermutation(maug_min1, scanfunction=mqmscan  #  change mqmscan() file
                               , cofactors=setcofactors # used with three cofactors coming from twoscan()
                               #, n.cluster=2
                               , n.perm=1000
                               #, batchsize=10
  )
  resultsrqtl_co1 = mqmprocesspermutation(results_co1)
  summary(resultsrqtl_co1, alpha = c(0.1, 0.05, 0.01, 0.001))
  
  # final QTL results
  summary(mqm_co1, perms = resultsrqtl_co1, alpha=0.05) # QTL on A06 is found again. Further, a third QTL is discovered on chr A06 (129 cM)

  
# Then set 3 cofactors
summary(mqm_min1, perms = resultsrqtl, alpha=0.05) # cofactors from MQM


# select peak markers from identified QTLs  
chr = c(2,3,6)

marker = c()
marker = c(marker, find.marker(data, chr = chr[1], pos = c(86)))
marker = c(marker, find.marker(data, chr = chr[2], pos = c(129)))
marker = c(marker, find.marker(data, chr = chr[3], pos = c(64)))
marker
# find.markerindex(maug_min1, marker)

multitoset_2 = find.markerindex(maug_min1, marker) # marker includes all three QTLs previously identified
setcofactors_2 = mqmsetcofactors(maug_min1, cofactors=multitoset_2)
mqm_co2 <- mqmscan(maug_min1, setcofactors_2, step.size = 1.0) # MQM analysis with one QTL as cofactor

  # calculation of LOD threshold significance
  results_co2 = mqmpermutation(maug_min1, scanfunction=mqmscan  #  change mqmscan() file
                               , cofactors=setcofactors_2 # used with three cofactors coming from twoscan()
                               #, n.cluster=2
                               , n.perm=1000
                               #, batchsize=10
  )
  
  resultsrqtl_co2 = mqmprocesspermutation(results_co2)
  summary(resultsrqtl_co2, alpha = c(0.1, 0.05, 0.01, 0.001))
  
  # final QTL results (compare all MQMscan() performed so far)
  summary(mqm_co2, perms = resultsrqtl_co2, alpha=0.05) # MQMscan() + 3 cofactors still finds 3 QTLs. No further improvement of the model.
  summary(mqm_co1, perms = resultsrqtl_co1, alpha=0.05) # MQMscan() + 1 cofactor 
  summary(mqm_min1, perms = resultsrqtl, alpha=0.05)    # MQMscan() without cofactor identied only 2 QTLs

  # plot results for comparison
  par(mfrow = c(2,1))
  plot(mqmgetmodel(mqm_co1))
  plot(mqm_co1, mqm_min1, col = c("blue", "red"),
       lim = c(0,round(max(summary(mqm_min1)[,3]))+1),
       ylab=c("LOD score"))
  abline(h=summary(resultsrqtl)[1], col=c("red")) # This line plots the signifcance threshold at 5% error rate
  abline(h=summary(resultsrqtl_co1)[1], col=c("blue")) # This line plots the signifcance threshold at 5% error rate
  dev.off()
  
  
### OVERVIEW - comparison of different MQM methods ###
summary(mqm_co2, perms = resultsrqtl_co2, alpha=0.05)
summary(mqm_co1, perms = resultsrqtl_co1, alpha=0.05) # These are final results reported in the manuscript
summary(mqm_min1, perms = resultsrqtl, alpha=0.05)

# export name of real markers to plot QTL profiles
markers = mqmextractmarkers(mqm_min1)
markers = mqmextractmarkers(mqm_co1)    # These are final results reported in the manuscript
markers = mqmextractmarkers(mqm_co2)
write.table(markers, file="mqm_lod_score.txt", sep="\t")


  

### study variance explained by QTLs ###

data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BLUE.csv", crosstype = "riself")
data = jittermap(data)
data_2 = calc.genoprob(data, step=1, error.prob=0.001)    # step = 1 USED


# Final QTLs (Pbc1, Pbc2, Pbc3) identified in this study
summary(mqm_co1, perms = resultsrqtl_co1, alpha=0.05)

# 1) check QTL effects + phenotypic variance
chr = c(2,3,6)
pos = c(86,129,64)

marker = c()
marker = c(marker, find.marker(data, chr = chr[1], pos = pos[1]))
marker = c(marker, find.marker(data, chr = chr[2], pos = pos[2]))
marker = c(marker, find.marker(data, chr = chr[3], pos = pos[3]))
marker

pos_marker = c()
for (i in 1:length(marker)){
  pos_marker = c(pos_marker, mqm_co2[row.names(mqm_co2)==marker[i],2])
}
pos_marker


pos_marker = pos

qtl <- makeqtl(data_2, chr=chr, pos=pos   # data_2 has
               #qtl <- makeqtl(data_2, chr=c(2,3,6),pos=c(87,21,65) 
               ,what="prob" # use in combination with HK method for fitting model
)
qtl
out.fq <- fitqtl(data_2, qtl=qtl, pheno.col=1, formula=y~Q1+Q2+Q3
                 ,method = "hk"
                 ,get.ests = TRUE
)
summary(out.fq)


# Condifence interval (1LOD)
# QTL Pbc1 on chromosome A02
lodint(rqtl, qtl.index = 1, expandtomarkers=TRUE)   # change to index to explore different QTLs
# QTL Pbc1 on chromosome A02
lodint(rqtl, qtl.index = 2, expandtomarkers=TRUE)   # change to index to explore different QTLs
# QTL Pbc1 on chromosome A02
lodint(rqtl, qtl.index = 3, expandtomarkers=TRUE)   # change to index to explore different QTLs

## NOTE: +/- 1.5 LOD confidence interval reported in manuscript were selected manually on the file, NOT with the formula abobe.  

# Export marker/LOD scores

# use mqmextractmarkers(mqm_min1) to remove pseudomarkers for final results

write.table(out.hk, file="marker_out.hk.txt", col.names = TRUE,sep="\t")  
write.table(mqm_min1, file="marker_mqm_min1.txt", col.names = TRUE,sep="\t")
write.table(mqm_co1, file="marker_mqm_co1.txt", col.names = TRUE,sep="\t")




### study effects of QTLs ###

## See the script Plot_figure_3.R, lines 125-226


                                     
