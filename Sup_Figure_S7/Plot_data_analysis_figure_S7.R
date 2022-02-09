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
setwd("paste_here_you_path_to_file")	  # set working directory where input datasets are locate


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


# these tests can be used to compare models 
AIC(lm1) # goodness-of-fit measure - smaller values are better
BIC(lm1) # goodness-of-fit measure - smaller values are better
coef(lm1) # coefficients of the model
confint = as.data.frame(confint(lm1)) # confidence interval
confint 


### PAIRWISE COMPARISONS

library(lsmeans)
library(multcomp)


lm1 <- lm(sqrt.Response ~ Genotype, data)

lsmeans(lm1, ~ Genotype, adjust = "Tukey")
marginal = lsmeans(lm1, ~ Genotype, adjust = "Tukey")

cld(marginal,
    alpha=0.01,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey")       ### Tukey-adjusted comparisons 



#### Supplementary fig.S7 - Re-evaluation of twelve selected RILs ####
#____________________________________________________________________#


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


# Prepare letters from multiple comparison
data.summarized = data %>% group_by(Genotype) %>% summarize(Response=max(Response))
hsd=HSD.test(aov(sqrt.Response~Genotype,data=data), "Genotype", alpha = 0.01, group=T)
hsd=hsd$groups
names(hsd)=c("Response", "groups")
hsd =hsd[order(match(rownames(hsd), data.summarized$Genotype)), , drop = FALSE]
hsd


# Boxplot for Egg_HR
par(mfrow=c(1,1))		# (row, column)
size = 17 # size labels 

ggplot(data = data) +
  geom_boxplot(aes(x=reorder(Genotype, Egg_HR_size_MAX), y=Response, fill=Genotype),
               outlier.colour=NA) +
  geom_point(mapping=aes(x=Genotype,y=Response,fill=Genotype),
             colour="black",
             size=1.5,
             shape=19,
             alpha=0.5,
             position=position_jitter(width=0.2,height=0)) +
  labs(x="", y="HR-like size (mm2)", size=size) +
  scale_fill_manual(values=rep("#33A02C",length(unique(data$Genotype)))) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  geom_text(data=data.summarized
            ,aes(x=Genotype,
                 y=0.2+Response,label=hsd$groups)
            ,size=size/3
            ,vjust=0)+
  theme_bw(base_size = 20) + # remove gray background 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  theme(legend.position="none") +
  theme(axis.text.x = element_text(colour="black", angle = 60, hjust = 1, size=size),
        axis.text.y = element_text(colour = "black", size=size),
        axis.title=element_text(colour = "black", size=size)) +
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5), limits = c(0, 2.5)) + # use 3.5 for egg, 3 for egg wash
  #scale_y_continuous(breaks=c(0, 2.5, 5, 7.5, 10, 12.5))
  # scale_y_continuous(breaks=c(0,2.5,5,7.5), limits=c(0,7.5)) # NOT USED
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot
ggsave("supp_figure_S7.png")
dev.off()

## NOTE: Final version of Sup. figure S7 was edited on PowerPoint to the QTL haplotypes of each RILs names. 
#  QTL haplotypes can be extracted from the file "dataset_Table_S15.txt"
