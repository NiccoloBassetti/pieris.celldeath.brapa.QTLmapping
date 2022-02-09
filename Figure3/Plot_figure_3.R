############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset for Figure 3A:  "rqtl_HR_BLUE.csv", results from QTL analysis store in "pieris.celldeath.brapa.QTLmapping.RData"
### Dataset for Figure 3B:  "figure_3B.csv"
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
setwd("paste_here_your_path_to_file")	  # set working directory where input datasets are located. 



####       Load results of QTL analysis      ####
#_______________________________________________#


## NOTE: results/objects generated in the following analysis and reported in the manuscript are store in the .RData file below.
#  Those results can be loaded to quickly generate plots and/or compare results from a new analysis. 
load("pieris.celldeath.brapa.QTLmapping.RData")





####       Figure 3 - panel A      #### 
#____________________________________#

library("qtl")
library("ggplot2")
library("ggrepel")


# For details on the analysis see "QTL_analysis_figure_3.R", lines 148-295 


## Prepare dataset to extract names and location of peak markers of QTLs Pbc1-3
data = read.cross(format = "csv", dir = ".",  file = "rqtl_HR_BLUE.csv", crosstype = "riself")

chr = c(2,3,6)        # these are the chromosomes on which QTLs Pbc1-3 were found
pos = c(86,129,64)  # these are the location of peak makers of Pbc1-3

marker = c()
marker = c(marker, find.marker(data, chr = chr[1], pos = pos[1]))
marker = c(marker, find.marker(data, chr = chr[2], pos = pos[2]))
marker = c(marker, find.marker(data, chr = chr[3], pos = pos[3]))
marker


## Prepare the LOD scores profiles (contained in QTL results "pieris.celldeath.brapa.QTLmapping.RData"

#lod_threshold_hk = summary(operm.hk, alpha = c(0.05))          # use this to plot single QTL analysis
lod_threshold_mqm = summary(resultsrqtl, alpha = c(0.05))       # use this to plot MQM analysis
#lod_threshold_mqm = summary(resultsrqtl_co1, alpha = c(0.05))   # use this to plot MQM analysis with forward selection (not reported in manuscript))
# marker_qtl = find.marker(data, chr = chr[i], pos = pos[i]) 

# compare scanone() with MQMscan()
#scanone_dataset = out.hk         # use this to plot single QTL analysis
mqm_dataset = mqm_co1

#single.qtl = data.frame(scanone_dataset, method = "Single QTL")      # use this to plot single QTL analysis
#head(single.qtl)
multi.qtl = data.frame(mqm_dataset[c(1,2,3)], method = "Multi QTL")
head(multi.qtl)
names(multi.qtl) = names(single.qtl)
head(multi.qtl)
#input = rbind(single.qtl,multi.qtl)

chr = c("A02","A03","A06")
multi.qtl[,1] = paste("A0", multi.qtl[,1], sep="") # change chromosome name
multi.qtl = multi.qtl[-grep("c", row.names(multi.qtl)),] # remove pseudomarkers
head(multi.qtl)


## Prepare labels with markers` names to be inlcuded within the plot

# labels_scanone <- as.data.frame(summary(scanone, perms = operm.hk, alpha = 0.05, pvalues = TRUE)) # use this to plot single QTL analysis

labels_df <-  as.data.frame(multi.qtl[row.names(multi.qtl)==marker[1],c(1,2,3)])
labels_df <-  rbind(labels_df, as.data.frame(multi.qtl[row.names(multi.qtl)==marker[2],c(1,2,3)]))
labels_df <-  rbind(labels_df, as.data.frame(multi.qtl[row.names(multi.qtl)==marker[3],c(1,2,3)]))
labels_df


## NOTE: The qtl_plot() function needs first to be generated. First run lines 239-342

plot = qtl_plot(#input, 
  multi.qtl,
  chrs = chr, 
  lod = lod_threshold_mqm[1], 
  ncol = 3,  
  rug = TRUE, 
  labels = labels_df)

# save plot
ggsave(filename=paste("QTL_main",".png", sep=""), device = "png", plot=plot, width = 13.5, height = 5, units = "in", dpi = 300)
dev.off()




####       Figure 3 - panel B      ####
#_____________________________________#

library("ggplot2")
library("svglite")
library("tidyverse")
library("RColorBrewer")


## Prepare dataset
#  NOTE: the dataset "figure_3B.csv" for this figure can be assembled from Rqtl input dataset "rqtl_HR_BLUE.csv" 
#  Phenotypic data are BLUEs for HR-like of each RIL (1st column, square root transformed)
#  Genotypic data are alleles of each peak markers of the QTLs Pbc1 (BrID11121), Pbc2 (BrID90099) and Pbc3 (BrID90095)
mydata <- read.table("figure_3B.csv", sep = ",", header = TRUE)

head(mydata)

# Search NAs
sapply(mydata, function(x)(sum(is.na(x))))		# search NAs
mydata <- na.omit(mydata)		# revome 

head(mydata)


mydata  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
# mydata <- mydata %>% filter(`staining_type`=="TB")
# mydata <- mydata %>% dplyr::select(Species, staining, HR, staining) %>% gather(HR, staining, key = class, value = value)
mydata$Genotype <- as.character(mydata$Genotype)
head(mydata)


## choose colours
cbPalette=c("Paired") # choose palette colours
n = 12 # select nymber of colours

brewer.pal(n = n, name = cbPalette)
display.brewer.pal(n = n, name = cbPalette)


# Violin plot
size=15

par(mfrow=c(1,1))		# open graphic device (row, column) 
plot = ggplot(mydata) + # use geom_boxplot or geom_violin
  geom_violin(data=mydata,
              aes(x=Genotype, y=HR_max, fill=Genotype),
              colour="black", # black outline for boxes
              outlier.colour=NA,
              width=0.6) +
  geom_boxplot(aes(x=Genotype, y=HR_max),
               colour="black", # black outline for boxes
               outlier.colour=NA,
               width=0.1) +
  geom_point(data=mydata[mydata$Genotype=="R-o-18",], 
             mapping=aes(x=Genotype,y=HR_max, fill=Genotype), 
             colour=c("#1F78B4"), # black outline for points
             size=2,
             shape=19,
             alpha=0.5,
             position=position_jitterdodge(jitter.width=0.40,jitter.height=0,dodge.width=1,seed=NA)) +
  geom_point(data=mydata[mydata$Genotype=="L58",], 
             mapping=aes(x=Genotype,y=HR_max, fill=Genotype), 
             colour=c("#E31A1C"), # black outline for points
             size=2,
             shape=19,
             alpha=0.5,
             position=position_jitterdodge(jitter.width=0.40,jitter.height=0,dodge.width=1,seed=NA)) +
  stat_summary(aes(x=Genotype, y=HR_max),
               fun.y = mean, 
               geom = "point", 
               col = "black",
               shape = 18, size = 4) + 
  facet_wrap(~ chr, ncol = 3, scales = "free_x") +  
  labs(x="Allele", y="HR-like (mm2)", size=size) +
  scale_fill_manual(values=rep("white",3)) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  # minimal plotting theme
  #theme_minimal(base_size = 25) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  theme(legend.position="none") +
  # increase strip title size
  theme(strip.text = element_text(face = "bold", size = 15)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size=size, color="black"),
        axis.text.y=element_text(size=size, color="black"),
        axis.title=element_text(size=size, color="black")) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  theme(legend.position="none") + # change to set position, "none" remove legend
  #theme(legend.position = c(0.85, 0.90)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size-10), face="plain")) + #+ # change "face" value for style text
  #scale_x_discrete(limits=unique(mydata$Genotype)) + # reorder scale x and legend
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), limits = c(0, 3.5))
#scale_y_continuous(breaks=c(200,seq(500,3000, by=500)),
#breaks=c(seq(0,200, by=50)),
#                  limits=c(200,3000)
#)+
#theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)

ggsave("Fig.3_QTL_effect.png", device = "png", plot=plot, width = 13.5, height = 5, units = "in", dpi = 300)

dev.off()




#### APPENDIX ####
#### qtl_plot() funtion
# modified from Shirin Elsinghorst - https://shiring.github.io/ggplot2/2017/02/12/qtl_plots?

qtl_plot <- function(input,              # data frame input from scanone
                     mult.pheno = FALSE, # multiple phenotypes?
                     model = "normal",   # model used in scanone
                     chrs = NA,          # chromosomes to display
                     lod = NA,           # LOD threshold
                     rug = FALSE,        # plot marker positions as rug?
                     ncol = NA,          # number of columns for facetting
                     labels = NA         # optional dataframe to plot QTL labels
) {
  
  # if we have multiple phenotypes and/or a 2part model, gather input
  if (mult.pheno & model == "2part") {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (mult.pheno) {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (model == "2part") {
    input <- gather(input, method, lod, lod.p.mu:lod.mu)
  }
  
  # if not all chromosomes should be displayed, subset input
  if (!is.na(chrs)[1]) {
    input <- input[as.character(input$chr) %in% chrs, ]
  }
  
  # if there is more than one LOD column, gather input
  if (!any(colnames(input) == "lod")) {
    input$lod <- input[, grep("lod", colnames(input))]
  }
  
  # if no number of columns for facetting is defined, plot all in one row
  if (is.na(ncol)) {
    ncol <- length(unique(input$chr))
  }
  
  # if labels are set and there is no name column, set from rownames
  if (!is.na(labels)[1]) {
    if (is.null(labels$name)) {
      labels$name <- rownames(labels)
    }
  }
  
  # plot input data frame position and LOD score
  plot <- ggplot(input, aes(x = pos, y = lod)) + {
    
    # if LOD threshold is given, plot as horizontal line
    if (!is.na(lod)[1] & length(lod) == 1) geom_hline(yintercept = lod, linetype = "dashed")
  } + {
    
    if (!is.na(lod)[1] & length(lod) > 1) geom_hline(data = lod, aes(yintercept = lod, linetype = group))
  } + {
    
    # plot rug on bottom, if TRUE
    if (rug) geom_rug(size = 0.1, sides = "b")
  } + {
    
    # if input has column method but not group, plot line and color by method
    if (!is.null(input$method) & is.null(input$group)) geom_line(aes(color = method), size = 2, alpha = 0.6)
  } + {
    
    # if input has column group but not method, plot line and color by group
    if (!is.null(input$group) & is.null(input$method)) geom_line(aes(color = group), size = 2, alpha = 0.6)
  } + {
    
    # if input has columns method and group, plot line and color by method & linetype by group
    if (!is.null(input$group) & !is.null(input$method)) geom_line(aes(color = method, linetype = group), size = 2, alpha = 0.6)
  } + {
    
    # set linetype, if input has columns method and group
    if (!is.null(input$group) & !is.null(input$method)) scale_linetype_manual(values = c("solid", "twodash", "dotted"))
  } + {
    
    # if input has neither columns method nor group, plot black line
    if (is.null(input$group) & is.null(input$method)) geom_line(size = 2, alpha = 0.6)
  } + {
    
    # if QTL positions are given in labels df, plot as point...
    if (!is.na(labels)[1]) geom_point(data = labels, aes(x = pos, y = lod))
  } + {
    
    # ... and plot name as text with ggrepel to avoid overlapping
    if (!is.na(labels)[1]) geom_text_repel(data = labels, aes(x = pos, y = lod, label = name), size=5, nudge_y = 0.5)
  } + 
    # facet by chromosome
    facet_wrap(~ chr, ncol = ncol, scales = "free_x") +
    # minimal plotting theme
    # theme_minimal(base_size = 20) +
    theme_bw(base_size = 20) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
    theme(legend.position="none") +
    # increase strip title size
    theme(strip.text = element_text(face = "bold", size = 20)) +
    # set y-axis values
    scale_y_continuous(breaks=c(0, 2, 4, 6), limits = c(0, 6)) +
    # use RcolorBrewer palette
    scale_color_brewer(palette = "Dark2") +
    #scale_fill_manual(values=c("green")) + 
    # Change plot labels
    labs(x = "Genetic map (cM)",
         y = "LOD",
         color = "",
         linetype = "")
  
  print(plot)
}
