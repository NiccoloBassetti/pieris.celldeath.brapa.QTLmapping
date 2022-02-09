############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "dataset_figure_S3.txt"
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
x=c("ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse")

lapply(x, install.packages, character.only = TRUE)

# Set up working directory
getwd()                 # see current working directory
setwd("paste_here_you_path_to_file")	  # set working directory where input datasets are located. 




#### Supplementary. figure S3 ####
# WEKA segmentation vs manual segmentation # 
#_________________________ ________________#



# Load datasets
mydata <- read.table(file = "dataset_figure_S3.txt",sep = "\t", header = TRUE)
names(mydata)
str(mydata)


# Transform factor written as "integer"
sapply(FUN = class, mydata)
mydata$RIL_genotype <- as.factor(mydata$RIL_genotype)
mydata$Method <- as.factor(mydata$Method)

str(mydata)
head(mydata)


mydata  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
mydata$RIL_genotype <- as.character(mydata$RIL_genotype)
head(mydata)

# choose colours
cbPalette=c("Paired") # choose palette colours
n = 12 # select nymber of colours

brewer.pal(n = n, name = cbPalette)
display.brewer.pal(n = n, name = cbPalette)

size=15


# Violin plot
par(mfrow=c(1,1))		# open graphic device (row, column) 
plot = ggplot(mydata) + # use geom_boxplot or geom_violin
  geom_boxplot(aes(x=RIL_genotype, y=Response, fill=Method),
               colour="black", # black outline for boxes
               outlier.colour=NA,
               width=1) +
  geom_point(data=mydata, 
             mapping=aes(x=RIL_genotype, y=Response, fill=Method), 
             colour=c("black"), # black outline for points
             size=2,
             shape=19,
             alpha=1,
             position=position_jitterdodge(jitter.width=0.01,jitter.height=0,dodge.width=1,seed=NA)) +
  coord_flip() +
  labs(x="Picture", y="Cell death size (mm2)", size=size) +
  scale_fill_manual(values=c("#1F78B4", "#E31A1C")) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  # minimal plotting theme
  # minimal plotting theme
  #theme_minimal(base_size = 25) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  # increase strip title size
  theme(strip.text = element_text(face = "bold", size = 15)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size=size, color="black"),
        axis.text.y=element_text(size=size, color="black"),
        axis.title=element_text(size=size, color="black")) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  #theme(legend.position="none") + # change to set position, "none" remove legend
  theme(legend.position = c(0.80, 0.20)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size-2), face="plain")) + #+ # change "face" value for style text
  #scale_x_discrete(limits=unique(mydata$Genotype)) + # reorder scale x and legend
  #scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), limits = c(0, 3.5))
  #scale_y_continuous(breaks=c(200,seq(500,3000, by=500)),
  #breaks=c(seq(0,200, by=50)),
  #                  limits=c(200,3000)
  #)+
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)
ggsave("Fig.S3_WEKA_reproducibility.png", device = "png", plot=plot, width = 13.5, height = 5, units = "in", dpi = 300)
#ggsave(file="Figure_PR1_butterflies_TOP.svg", plot=plot, width=6, height=6, units="in")
dev.off()

