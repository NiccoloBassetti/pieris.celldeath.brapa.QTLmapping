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




####       Supplementary figure S5      ####              
#__________________________________________#


library("qtl")

### plot results plot results from scanone() and mqm scan()
#   comparison between single QTL and Multi-QTL (MQM) analysis


# For details on the analysis see "QTL_analysis_figure_3.R", lines 74-97 and 148-295 




# 1) plot genome wide LOD scores (all chromosomes)

scanone = out.hk
threshold.scanone = operm.hk    # results single QTL analysis

mqm = mqm_min1
threshold.mqm = resultsrqtl     # results MQM scan (backward model selection)

#mqm = mqm_co1                   # these results are similar to MQM scan and were not reported in the manuscript
#threshold.mqm = resultsrqtl_co1 # results MQM with (forward model selection)


# save plot to file
png(file=paste("QTL_overview",".png", sep=""), width=6000,height=3000,res=300) #, units="in", width=5, height=5, res=300)
par(mfrow=c(4,3))
par(mar=c(5,5,5,5))

chr = unique(scanone[,1])
for (i in 1:length(chr)) {  
  #png(file=paste("A0",i,".png", sep=""), width=3000,height=2000,res=300) #, units="in", width=5, height=5, res=300)
  max.cM.scale = round(max(scanone[scanone[,1]==i,2])/20)*20
  max.lod.score = round(max(summary(mqm)[,3]))+1
  
  plot(scanone
       , mqm
       #, mqm_co1
       , col = c("blue"
                 , "red"
                 #, "green"
       )
       , chr = i
       , type = "l"
       #, bandcol="gray70"
       , ylim = c(0,max.lod.score)
       , ylab=c("")
       , xlab=c("")
       , xaxt='n'
       , yaxt='n'
  )
  
  title(main=paste("A",i, sep="")
        #, sub="sub-title"
        , xlab="Map position (cM)"
        , ylab="LOD score"
        , font.main=2, font.lab=2, font.sub=2
        , cex.main=2.5, cex.lab=1.5)
  #mtext("LOD score", side=2, line=2.2, cex=1)
  axis(1, seq(0,max.cM.scale,20), font=2, cex.axis=1.5)
  axis(2, seq(0,max.lod.score,2), font=2, cex.axis=1.5)
  
  abline(h=summary(threshold.scanone)[1], col=c("blue"), lwd=3, lty=3)
  abline(h=summary(threshold.mqm)[1], col=c("red"), lwd=3, lty=3) 
  #abline(h=summary(resultsrqtl_co1)[1], col=c("red"), lwd=3, lty=3)
  
  if (i==round(as.numeric(chr[length(chr)])/3)){  # all chromosomes
    x.axis.lim = max(mqm[mqm[,1]==i,2])-55
    y.axis.lim = round(max(summary(mqm)[,3]))+1
    legend(x.axis.lim, y.axis.lim, legend=c("scanone", "MQM"),
           col=c("blue", "red"), lwd=3, lty=1, cex=1.5, box.lty=0)
  }else{}
}
dev.off()



#### EXTRA figure ####
#____________________#

# 2) plot only 3 chromosomes with QTLs
png(file=paste("QTL_main",".png", sep=""), width=3000,height=1000,res=300) #, units="in", width=5, height=5, res=300)

par (mfrow=c(1,3))
chr = c(2,3,6)    # select only the three QTLs identified on chromosomes A02, A03 and A06
for (i in chr) {
  max.cM.scale = round(max(out.hk[out.hk[,1]==i,2])/20)*20
  max.lod.score = round(max(summary(mqm_co1)[,3]))+1
  
  plot(mqm_co1        # results from MQM with cofactors (forward model selection), non reported in manuscript
       , out.hk       # results from single QTL analysis
       , mqm_min1     # results from MQM analysis
       , col = c("green"
                 , "blue"
                 , "red"
       )
       , chr = i
       , type = "l"
       #, bandcol="gray70"
       , ylim = c(0,max.lod.score)
       , ylab=c("")
       , xlab=c("")
       , xaxt='n'
       , yaxt='n'
  )
  title(main=paste("A",i, sep="")
        #, sub="sub-title"
        , xlab="Map position (cM)"
        , ylab="LOD score"
        , font.main=2, font.lab=2, font.sub=2
        , cex.main=3, cex.main=2.5, cex.main=2.5)  
  axis(1, seq(0,max.cM.scale,20), font=2, cex=3)
  axis(2, seq(0,max.lod.score,2), font=2, cex=3)
  abline(h=summary(resultsrqtl)[1], col=c("red"), lwd=3, lty=3) 
  abline(h=summary(operm.hk)[1], col=c("blue"), lwd=3, lty=3)
  abline(h=summary(resultsrqtl_co1)[1], col=c("green"), lwd=3, lty=3)
  
  if (i==chr[length(chr)]){              # 3 qtls
  x.axis.lim = max(out.hk[out.hk[,1]==i,2])-120
  y.axis.lim = round(max(summary(mqm_co1)[,3]))+1
  legend(x.axis.lim, y.axis.lim, 
         legend=c("scanone", 
                  "MQM",
                  "MQM_cofactors"),
         col=c("blue", 
               "red",
               "green"), 
         lwd=3, lty=3, cex=1.5, box.lty=0)
  }else{
  
  }
}
dev.off()