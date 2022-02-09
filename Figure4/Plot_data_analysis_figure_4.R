############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "dataset_Table_S15.txt"
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



## Load dataset

mydata <- read.table("dataset_Table_S15.txt", sep = "\t", header = TRUE)

head(mydata)
str(mydata)


# Search NAs
sapply(mydata, function(x)(sum(is.na(x))))		# search NAs
mydata <- na.omit(mydata)		# revome 

head(mydata)


# Add column with data transformation
trait = "Egg_HR_size_MAX" # select phenotypic trait to analyse

mydata$HR_max = as.numeric(mydata$HR_max)
mydata$Response = mydata[,which(names(mydata)==trait)] # create column with raw phenotypic data
mydata$sqrt.Response = sqrt(mydata[,which(names(mydata)==trait)]) # create column wiht square root-transformed data
head(mydata)

# Rename columns
for (i in c(2,3,4)){
  names(mydata)[i]=substr(names(mydata)[i],0,4) # shorten name of columns names
}

head(mydata)
tail(mydata)


# Transform categorical variables into "factors"
sapply(FUN = class, mydata)
mydata$Genotype <- factor(mydata$Genotype)
mydata$Pbc1 <- factor(mydata$Pbc1)
mydata$Pbc2 <- factor(mydata$Pbc2)
mydata$Pbc3 <- factor(mydata$Pbc3)
str(mydata)


#### Statistical analysis ####
#____________________________#

data = mydata

# diagnostic plots
par(mfrow=c(2,2)) 

lm1 <- lm(Response ~ Genotype, data)
plot(lm1)

lm2 <- lm(sqrt.Response ~ Genotype, data)
plot(lm2) # square root transformed data appear better suited for normality assumptions


# CHECK ASSUMPTIONS (normality and homogeneity of variance)
shapiro.test(data$Response) # check normality of residuals (small P-value = NOT NORMAL)
fligner.test(data$Response~data$Genotype)	# test for equal variances

shapiro.test(data$sqrt.Response) # check normality of residuals (small P-value = NOT NORMAL)
fligner.test(data$sqrt.Response~data$Genotype)

require(car)
outlierTest(lm1) # test for outliers
# lm2<-update(lm1, subset(data, mydata$Genotype !=57)) # remove outliers/not working



### LINEAR MODEL 1

data = mydata
#data = na.omit(data) # revome NAs for plotting
head(data)
dim(data)

# test effect of RIL genotype on phenotypic differences
lm1 <- lm(sqrt.Response ~ Genotype, data)
anova(lm1) # (reported in manuscript)


# effect of specific QTL haplotype on phenotypic differences
lm2 <- lm(sqrt.Response ~ Pbc1+Pbc2+Pbc3, data)
lm3 <- lm(sqrt.Response ~ Pbc1*Pbc2+Pbc1*Pbc3, data)

lm4 <- lm(sqrt.Response ~ Pbc1+Pbc3, data)

anova(lm2) # Pbc2 is not significant 
anova(lm3) # interactions terms are not signifcant

anova(lm4) # more parsimonius model (reported in manuscript)

# these tests can be used to compare models 
AIC(lm1) # goodness-of-fit measure - smaller values are better
BIC(lm1) # goodness-of-fit measure - smaller values are better
coef(lm1) # coefficients of the model
confint = as.data.frame(confint(lm1)) # confidence interval
confint 


### PAIRWISE COMPARISONS

library(lsmeans)
library(multcomp)

# Pbc1 vs Pbc3
lm1 <- lm(sqrt.Response ~ Pbc1+Pbc3, data)

lsmeans(lm1, ~ Pbc1+Pbc3, adjust = "Tukey")
marginal = lsmeans(lm1, ~ Pbc1+Pbc3, adjust = "Tukey")

cld(marginal,
    alpha=0.01,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey")       ### Tukey-adjusted comparisons 


# Pbc1 vs Pbc2
lm1 <- lm(sqrt.Response ~ Pbc1+Pbc2, data)

lsmeans(lm1, ~ Pbc1+Pbc2, adjust = "Tukey")
marginal = lsmeans(lm1, ~ Pbc1+Pbc2, adjust = "Tukey")

cld(marginal,
    alpha=0.01,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey")       ### Tukey-adjusted comparisons 




#### Fig.4 - Re-evaluation of twelve selected RILs ####
#_____________________________________________________#

# select packages
x=c("ggplot2", "agricolae", "RColorBrewer", "tidyverse", "svglite")   # load packages required for final plots
lapply(x, library, character.only = TRUE)


# Replace allele name for visualization purposes
data  <-  as_tibble(mydata)
options(pillar.sigfig = 5) # visualize decimals in tibble

data = data  %>% mutate(Pbc1 = str_replace(Pbc1, "B", "R"))
data = data  %>% mutate(Pbc2 = str_replace(Pbc2, "B", "R"))
data = data  %>% mutate(Pbc3 = str_replace(Pbc3, "B", "R"))
data = data  %>% mutate(Pbc1 = str_replace(Pbc1, "A", "L"))
data = data  %>% mutate(Pbc2 = str_replace(Pbc2, "A", "L"))
data = data  %>% mutate(Pbc3 = str_replace(Pbc3, "A", "L"))

head(data)

# choose colours
cbPalette=c("Paired") # choose palette colours
n = 12 # select nymber of colours

brewer.pal(n = n, name = cbPalette)
display.brewer.pal(n = n, name = cbPalette)



size=15


# reorganize order levels for visualization purposes
data$Pbc1=factor(data$Pbc1, 
                   levels = c("R","L"), ordered = TRUE)
data$Pbc2=factor(data$Pbc2, 
                   levels = c("L","R"), ordered = TRUE)
data$Pbc3=factor(data$Pbc3, 
                   levels = c("L","R"), ordered = TRUE)
str(data)



#### Plot Pbc1 vs Pbc2

plot = ggplot(data) + # use geom_boxplot or geom_violin
  geom_boxplot(aes(x=Pbc1, y=Response, fill=Pbc2),
               colour="black", # black outline for boxes
               outlier.colour=NA,
               width=1) +
  geom_point(data=data, 
             mapping=aes(x=Pbc1, y=Response, fill=Pbc2), 
             colour=c("black"), # black outline for points
             size=2,
             shape=19,
             alpha=1,
             position=position_jitterdodge(jitter.width=0.20,jitter.height=0,dodge.width=1,seed=NA)) +
  stat_summary(aes(x=Pbc1, y=Response, fill=Pbc2),
               fun.y = mean,
               position=position_dodge(width=0.80),
               geom = "point", 
               col = "white",
               shape = 18, size = 4) +
  labs(x="Pbc1 allele", y="HR-like (mm2)", size=size)  +
  scale_fill_manual(name = c("Pbc2 allele"), values=c("#1F78B4", "#E31A1C")) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
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
  #theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size-2), face="plain")) + #+ # change "face" value for style text
  #scale_x_discrete(limits=unique(data$Genotype)) + # reorder scale x and legend
  #scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), limits = c(0, 3.5))
  #scale_y_continuous(breaks=c(200,seq(500,3000, by=500)),
  #breaks=c(seq(0,200, by=50)),
  #                  limits=c(200,3000)
  #)+
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)

# save figure
ggsave("figure_4A_Pbc1vsPbc2.png", device = "png", plot=plot, width = 13.5, height = 5, units = "in", dpi = 300)
dev.off()



#### Plot Pbc1 vs Pbc3


plot = ggplot(data) + # use geom_boxplot or geom_violin
  geom_boxplot(aes(x=Pbc1, y=Response, fill=Pbc3),
               colour="black", # black outline for boxes
               outlier.colour=NA,
               width=1) +
  geom_point(data=data, 
             mapping=aes(x=Pbc1, y=Response, fill=Pbc3), 
             colour=c("black"), # black outline for points
             size=2,
             shape=19,
             alpha=1,
             position=position_jitterdodge(jitter.width=0.20,jitter.height=0,dodge.width=1,seed=NA)) +
  stat_summary(aes(x=Pbc1, y=Response, fill=Pbc3),
               fun.y = mean,
               position=position_dodge(width=0.80),
               geom = "point", 
               col = "white",
               shape = 18, size = 4) +
  labs(x="Pbc1 allele", y="HR-like (mm2)", size=size)  +
  scale_fill_manual(name = c("Pbc3 allele"), values=c("#1F78B4", "#E31A1C")) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
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
  #theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", size=(size-2), face="plain")) + #+ # change "face" value for style text
  #scale_x_discrete(limits=unique(data$Genotype)) + # reorder scale x and legend
  #scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), limits = c(0, 3.5))
  #scale_y_continuous(breaks=c(200,seq(500,3000, by=500)),
  #breaks=c(seq(0,200, by=50)),
  #                  limits=c(200,3000)
  #)+
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
plot(plot)

# save figure
ggsave("figure_4B_Pbc1vsPbc3.png", device = "png", plot=plot, width = 13.5, height = 5, units = "in", dpi = 300)
dev.off()

