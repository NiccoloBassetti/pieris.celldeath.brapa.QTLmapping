############################################################
### Script for statistical analysis and plotting of data ###
### 
### Dataset:  "dataset_TableS12.txt"
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
x=c("ggplot2", "agricolae", "RColorBrewer", "car", "multcomp", "lme4", "lsmeans", "tidyverse")

lapply(x, library, character.only = TRUE)
lapply(x, library, character.only = TRUE)

# Set up working directory
getwd()                 # see current working directory
setwd("path_to_file")	  # set working directory where input datasets are located. 


# Load datasets
mydata <- read.table(file = "dataset_TableS12.txt",sep = "\t", header = TRUE) 
names(mydata)
str(mydata)


# Transform categorical variables into "factors"
sapply(FUN = class, mydata)
mydata$Block <- as.factor(mydata$Block)
mydata$Row <- as.factor(mydata$Row)
mydata$Col <- as.factor(mydata$Col)
mydata$Plots <- as.factor(mydata$Plots)
mydata$Gen_code <- as.factor(mydata$Gen_code)
mydata$Genotype = as.factor(mydata$Genotype)
mydata$Plant_ID = as.factor(mydata$Plant_ID)

sapply(FUN = class, mydata)
str(mydata)

# Examine dataset - remove NAs
dim(mydata)   # check size dataset
sapply(mydata, function(x)(sum(is.na(x))))		# search NAs

mydata <- na.omit(mydata)		# revome NAs
sapply(mydata, function(x)(sum(is.na(x))))		# check if removal NAs worked

dim(mydata)   # compare with original size dataset



###############################
### Correlation between traits
###
data = mydata
data = data[!(data$Genotype=="B.nigra"),] # exclude B.nigra from data analysis because positive control
dim(data)

# correlation matrix
traits = names(data)[c(8,9)]
traits
data_cor = data[traits]
cor_paerson = round(cor(data_cor, method = c("pearson")),3)
cor_paerson

# P-value
cor.test( ~ Egg_HR_size_MAX_plant + Wash_HR_size_MAX_plant,
          data=data,
          method = "pearson",
          conf.level = 0.95)



############################
### Statistical analysis ###
###


### Select phenotypic trait to be analysed
names(mydata)[c(8,9)]

trait = "Egg_HR_size_MAX_plant" # repeat analysis with this trait
#trait = "Wash_HR_size_MAX_plant" # repeat analysis with this trait

mydata$Response = mydata[,which(names(mydata)==trait)]    # add a column "Response" with values of "trait" to be analysed
mydata$sqrt.Response = sqrt(mydata$Response)              # add a column with squared root transformation of "Response"
head(mydata)


### Generate summary statistics of non-transformed data
#   (Supplementary Table S3)

summary(mydata)
N = table(mydata$Genotype)
mean = tapply(mydata$Response, mydata$Genotype, mean)	 # mean variable for each genotype
sd = tapply(mydata$Response, mydata$Genotype, sd)
se = sd/sqrt(N)
tapply(mydata$Response, mydata$Genotype, var)
min = tapply(mydata$Response, mydata$Genotype, min)
max = tapply(mydata$Response, mydata$Genotype, max)

summary_stat = as.data.frame(round(cbind(N,min,max,mean,sd,se),2))
summary_stat = summary_stat[order(summary_stat$mean),]
summary_stat

# export summary statistics dataframe
file_name = paste("summary_stat_", trait, ".txt", sep = "")
write.table(summary_stat, file = file_name, sep="\t")   # this .txt was imported in Excel to edit final table layout



### LINEAR MODEL 1 - analysis on non-trasformed data
data = mydata
data = data[!(mydata$Genotype=="B.nigra"),] # to exclude B.nigra because from plots
data = na.omit(data) # revome NAs for plotting
dim(data)

#lm1 <- lm(Response ~ Block + Block:Col + Block:Row + Genotype, data)
#lm1 <- lm(Response ~ Block + Genotype, data)
lm1 <- lm(Response ~ Genotype, data)

# summary of fitted models (compare all three models above)
summary(lm1)
anova(lm1)
AIC(lm1) # goodness-of-fit measure 
BIC(lm1) # goodness-of-fit measure 
coef(lm1) # coefficients of the model
confint = as.data.frame(confint(lm1)) # confidence interval
confint 
# write.table(confint, file = "confint_model.txt", header=TRUE)

# Note: There is ittle evidence that experimental design terms (Blocks, Row and Col) improve model. 
#       More parsimonious model only including Genotype term has better fit (lower AIC and BIC).



### CHECK ASSUMPTIONS (normality and homogeneity of variance)

# Check distribution of non-transformed data
n = 100 # numbers of breaks
hist(mydata$Response, breaks=n)
# histInfo <- hist(mydata$Response, breaks=n)
# histInfo
lines(density(mydata$Response))
rug(jitter(mydata$Response))


# diagnostic plots for non-transformed data
par(mfrow=c(2,2)) 

lm1 <- lm(Response ~ Genotype, data)
plot(lm1)

require(car)
outlierTest(lm1) # test for outliers

dev.off()

# test assumptions (normality and homogeneity of variance)
shapiro.test(data$Response)   # test normality of residuals (small P-value = NOT NORMAL)
fligner.test(data$Response~data$Genotype) # test for homogeneity of residuals variance

# test assumptions on transformed data
shapiro.test(data$sqrt.Response) # check normality of residuals (small P-value = NOT NORMAL)
fligner.test(data$sqrt.Response~data$Genotype)

# Note: non-transformed data look approximately normally distributed, but residual variance not equal. 
#       Data transformation with square root function improve the assumptions of linearity.
#       Proceed with data analysis on square root-transformed data. 



### LINEAR MODEL 2 - analysis on trasformed data
#   (final results reported in publication)

#lm1 <- lm(sqrt.Response ~ Block + Block:Col + Block:Row + Genotype, data)
#lm1 <- lm(sqrt.Response ~ Block + Genotype, data)
lm1 <- lm(sqrt.Response ~ Genotype, data)

# summary of fitted models (compare all three models above)
summary(lm1)
anova(lm1)
AIC(lm1) # goodness-of-fit measure 
BIC(lm1) # goodness-of-fit measure 
coef(lm1) # coefficients of the model
confint = as.data.frame(confint(lm1)) # confidence interval
confint 
# write.table(confint, file = "confint_model.txt", header=TRUE)

# Note: Again, little evidence that exp design terms (Blocks, Row and Col) improve model. 
#       Take more parsimonious model by only including Genotype term.



### Pairwise multiple comparisons (reported in Supplementary Table S3)
lm1 <- lm(sqrt.Response ~ Genotype, data) # this is the final model

require(lsmeans)
marginal = lsmeans(lm1, ~ Genotype, adjust = "Tukey")

require(multcomp)
cld(marginal,
    alpha=0.01,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey")       ### Tukey-adjusted comparisons 



### Estimate of variance of model terms 
### (reported in Supplementary Table S3)

# LINEAR MIXED EFFECTS MODELS (lme)
require(lme4)

lmer1 = lmer(sqrt.Response ~ (1|Genotype), data)
summary(lmer1)

as.data.frame(VarCorr(lmer1))$vcov    # check that column "vcov" reports same estimated variance as in model lmer1

var_genotype = as.data.frame(VarCorr(lmer1))$vcov[1]  # extract variance of Genotype term 
var_environment = as.data.frame(VarCorr(lmer1))$vcov[2]  # extract variance of error term


# calculate coefficient of variation (Supplementary Table S3)
grand_mean_trait = mean(data$sqrt.Response)

CV_genotype = sqrt(var_genotype)/grand_mean_trait * 100
CV_genotype

CV_enviroment = sqrt(var_environment)/grand_mean_trait * 100
CV_enviroment 


# calculate broad-sense heritability

broad_h2 = var_genotype/(var_genotype+var_environment) # EGGs
broad_h2 

# export results
df = round(as.data.frame(cbind(CV_genotype, CV_enviroment, broad_h2)),2)
trait
row.names(df) = trait
df

file_name = paste("analysis_variance_", trait, ".txt", sep = "")
write.table(df, file = file_name, sep="\t")   # this .txt was imported in Excel to edit final table layout


    
################  
### Boxplot data
###
### Fig.1
  
# select pacakges
x=c("ggplot2", "agricolae", "RColorBrewer", "tidyverse")   # load packages required for final plots
lapply(x, library, character.only = TRUE)

# select dataset
data = mydata
data = data[!(mydata$Genotype=="B.nigra"),] # to exclude B.nigra because from plots
data = na.omit(data) # revome NAs for plotting
dim(data)

# choose color Palette (requires RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE) # display all colorblind friendly palettes - https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/#use-a-colorblind-friendly-palette
cbPalette=c("Paired")   # choose a palette
n = 12    # select n of colours to be visualized

brewer.pal(n = n, name = cbPalette) # display names of colours
display.brewer.pal(n = n, name = cbPalette)  # plot colours


# Prepare letters from multiple comparison (requires agricolae and tidyverse)
data.summarized = data %>% group_by(Genotype) %>% summarize(Response=max(Response))
hsd=HSD.test(aov(sqrt.Response~Genotype,data=data), "Genotype", alpha = 0.01, group=T) # calculate Tukey`s HSD that can be added on top of boxplots
hsd=hsd$groups
names(hsd)=c("Response", "groups")
hsd = hsd[order(match(rownames(hsd), data.summarized$Genotype)), , drop = FALSE]
hsd


# Boxplot for mydata$Response (non-transformed data)
par(mfrow=c(1,1))		# (row, column)
size = 17 # size labels 

ggplot(data = data) +
	geom_boxplot(aes(x=reorder(Genotype, Egg_HR_size_MAX_plant), y=Response, fill=Genotype),
		outlier.colour=NA) +
	geom_point(mapping=aes(x=Genotype,y=Response,fill=Genotype),
	  colour="black",
	  size=1.5,
	  shape=19,
	  alpha=0.5,
		position=position_jitter(width=0.2,height=0)) +
	labs(x="", y="HR-like size (mm2)", size=size) +
  scale_fill_manual(values=c("#FF7F00", "#A6CEE3","#B2DF8A", "#6A3D9A", "#FB9A99", "#33A02C",
                             "#E31A1C", "#B15928", "#1F78B4", "#FDBF6F")) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_fill_brewer(palette=cbPalette) + # set color for fills - http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  #scale_colour_brewer(palette=cbPalette) + # set color for line/points - 
  geom_text(data=data.summarized
            ,aes(x=Genotype,
                 y=0.2+Response,label=hsd$groups)
            ,size=size/3
            ,vjust=0)+
  theme_bw() + # remove gray background 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid https://felixfan.github.io/ggplot2-remove-grid-background-margin/
	theme(legend.position="none") +
  theme(axis.text.x = element_text(colour="black", angle = 60, hjust = 1, size=size),
        axis.text.y = element_text(colour = "black", size=size),
        axis.title=element_text(colour = "black", size=size)) +
  scale_y_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), limits = c(0, 3)) + # use 3.5 for egg, 3 for egg wash
  #scale_y_continuous(breaks=c(0, 2.5, 5, 7.5, 10, 12.5))
	# scale_y_continuous(breaks=c(0,2.5,5,7.5), limits=c(0,7.5)) # NOT USED
  theme(aspect.ratio = 20/20) # fix the x/y ratio of the plot

trait
plot_name = paste("fig1_", trait, ".png", sep="")
ggsave(plot_name)

dev.off()
