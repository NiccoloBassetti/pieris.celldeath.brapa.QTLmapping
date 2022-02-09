############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "dataset_Table_S13.txt"
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

library(lme4); library(lmerTest); library(pbkrtest); library(MASS)
x=c("ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse")

lapply(x, install.packages, character.only = TRUE)

# Set up working directory
getwd()                 # see current working directory
setwd("paste_here_your_path_to_file")	  # set working directory where input datasets are located. 



#### Statistical analysis ####
#____________________________#

mydata <- read.table(file = "dataset_Table_S13.txt",sep = "\t", header = TRUE) 

str(mydata)
head(mydata)

# Transform categorical variables into "factors"
sapply(FUN = class, mydata)
mydata$Genotype <- factor(mydata$Genotype)
mydata$Egg_HR_size_MAX_plant_Exp_1 <- as.numeric(mydata$Egg_HR_size_MAX_plant_Exp_1)
mydata$Egg_HR_size_MAX_plant_Exp_2 <- as.numeric(mydata$Egg_HR_size_MAX_plant_Exp_2)
mydata$Egg_HR_size_MAX_plant_Exp_3 <- as.numeric(mydata$Egg_HR_size_MAX_plant_Exp_3)
mydata$Egg_HR_size_MAX_BLUE <- as.numeric(mydata$Egg_HR_size_MAX_BLUE)

head(mydata)



### Calculate BLUEs for multienvironment phenotypic measurements

# prepare dataset
require(tidyverse)
data = as_tibble(mydata)[,-5] # remove column with BLUE values
names(data)[c(2,3,4)] = c("A", "B", "C") # rename QTL experiments for easiness
names(data)
data = data %>% gather("A", "B", "C", key = "Experiment", value = "HR_MAX")
#data <- na.omit(data)

head(data)
tail(data)


# linear mixed model
require(lme4)
mm1 <- lmer(HR_MAX ~ Genotype + (1|Experiment), data = data)

anova(mm1)
summary(mm1)

BLUE <- fixef(mm1) # extract BLUE values (estimated fixed effects from mixed model)


# Prepare BLUEs values for export 
head(BLUE)
tail(BLUE)
str(BLUE)
length(BLUE)
names(BLUE) = substring(names(BLUE),9) # remove "Genotype" label from strings
for (i in 2: length(BLUE)) {
  BLUE[i] = BLUE[i] + BLUE[1]
}
names(BLUE)[1] = "L58"
head(BLUE)

# Note: These are the final BLUE phenotypic values used for QTL analysis. 
# They are attached to "dataset_TableS_13.txt" as column "Egg_HR_size_MAX_BLUE".



### heritability BLUE values - calculated on square root transformed data

# prepare dataset
require(tidyverse)
#data = mydata[c(1:162),c(1,5)] # c(1:162) include only one BLUE value per RIL/Genotype

head(data) # check that same dataset as above in stored in object "data"

data$HR_MAX_sqrt = sqrt(data$HR_MAX)

dim(data)
head(data)
tail(data)


# square_root
mm_h2 <- lmer(HR_MAX_sqrt ~ (1|Genotype) + (1|Experiment), data = data)
anova(mm_h2)
summary(mm_h2)

var_genotype = as.data.frame(VarCorr(mm_h2))$vcov[1]
var_experiment = as.data.frame(VarCorr(mm_h2))$vcov[2]
var_environment = as.data.frame(VarCorr(mm_h2))$vcov[1]

broad_h2 = var_genotype/(var_genotype+var_experiment+var_environment)
broad_h2 # 0.44



#### Fig.2 - Phenotypic data distribution   ####
#______________________________________________#

### Boxplot data

library(ggplot2)
library(tidyverse)
library(RColorBrewer)


## NOTE: "Phenotypic_data.csv" is assembled ONLY for visualization purpose.
#
#  Data includes the BLUEs values from column "Egg_HR_size_MAX_BLUE" of dataset "dataset_TableS13.txt"
#  and the raw data of the two parents (L58, RILs) across all Experiments 1,2,3.


mydata <- read.table("figure_2.CSV", sep = ",", header = TRUE) 
head(mydata)


mydata  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble
# mydata <- mydata %>% filter(`staining_type`=="TB")
# mydata <- mydata %>% dplyr::select(Species, staining, HR, staining) %>% gather(HR, staining, key = class, value = value)
mydata$Genotype <- as.character(mydata$Genotype)
head(mydata)

# Prepare color palette
require(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE) # display all colorblind friendly palettes - https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#use-a-colorblind-friendly-palette
cbPalette=c("Paired") # choose palette colours
n = 12 # select nymber of colours

brewer.pal(n = n, name = cbPalette)
display.brewer.pal(n = n, name = cbPalette)

brewer.pal(n = n, name = "Dark2")
display.brewer.pal(n = n, name = "Dark2")


# Boxplot
par(mfrow=c(1,1))		# open graphic device (row, column) 

size=15

plot = ggplot(mydata) + # use geom_boxplot or geom_violin
  geom_violin(data=mydata[mydata$Genotype=="RIL",],
              aes(x=Genotype, y=HR_max, fill=Genotype),
              colour="black", # black outline for boxes
              outlier.colour=NA,
              width=0.8) +
  geom_boxplot(data=mydata,
               aes(x=reorder(Genotype,HR_max), y=HR_max), #reorder(Genotype, Egg_1_HR_size)
               colour="black", # black outline for boxes
               outlier.colour=NA,
               width=0.1) +
  #coord_flip() +
  geom_point(data=mydata[mydata$Point!="mean_BLUE",][mydata[mydata$Point!="mean_BLUE",]$Genotype=="L58",], 
             mapping=aes(x=Genotype,y=HR_max, fill=Genotype), 
             colour=c("#E31A1C"), # black outline for points
             size=2,
             shape=19,
             alpha=0.5,
             position=position_jitterdodge(jitter.width=0.05,jitter.height=0,dodge.width=1,seed=NA)) +
  geom_point(data=mydata[mydata$Point!="mean_BLUE",][mydata[mydata$Point!="mean_BLUE",]$Genotype=="RIL",], 
             mapping=aes(x=Genotype,y=HR_max, fill=Genotype), 
             colour=c("#33A02C"), # black outline for points
             size=2,
             shape=19,
             alpha=0.5,
             position=position_jitterdodge(jitter.width=0.40,jitter.height=0,dodge.width=1,seed=NA)) +
  geom_point(data=mydata[mydata$Point!="mean_BLUE",][mydata[mydata$Point!="mean_BLUE",]$Genotype=="R-o-18",], 
             mapping=aes(x=Genotype,y=HR_max, fill=Genotype), 
             colour=c("#1F78B4"), # black outline for points
             size=2,
             shape=19,
             alpha=0.5,
             position=position_jitterdodge(jitter.width=0.05,jitter.height=0,dodge.width=1,seed=NA)) +
  geom_point(data=mydata[mydata$Point=="mean_BLUE",], 
             mapping=aes(x=Genotype,y=HR_max, fill=Genotype), 
             colour=c("black"), # black outline for points
             size=4,
             shape=18,
             #alpha=0.5,
             position=position_jitterdodge(jitter.width=0.05,jitter.height=0,dodge.width=1,seed=NA)) +
  #stat_summary(data = mydata[,c(1,2)],
  #            aes(x=Genotype, y=HR_max),
  #           fun.y = mean, 
  #          geom = "point", 
  #         col = "white",
  #        shape = 18, size = 4) +
  labs(x="Allele", y="HR-like (mm2)", size=size) +
  scale_fill_manual(values=c("white", "white", "white")) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  theme_bw() + # remove gray background 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size=size, color="black"),
        axis.text.y=element_text(size=size, color="black"),
        axis.title=element_text(size=size, color="black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  theme(legend.position="none") + # change to set position, "none" remove legend
  #theme(legend.position = c(0.85, 0.90)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size-10), face="plain")) + # change "face" value for style text
  scale_x_discrete(limits=unique(mydata$Genotype)) + # reorder scale x and legend
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), limits = c(0, 3)) +
  #scale_y_continuous(breaks=c(200,seq(500,3000, by=500)),
  #breaks=c(seq(0,200, by=50)),
  #                  limits=c(200,3000)
  #)+
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)
ggsave("fig.2_distribution_HR_BLUE.png", plot=plot, width = 6, height = 6, units="in", dpi = 300) # CHANGE
#ggsave(file="Figure_PR1_butterflies_TOP.svg", plot=plot, width=6, height=6, units="in")
dev.off()
