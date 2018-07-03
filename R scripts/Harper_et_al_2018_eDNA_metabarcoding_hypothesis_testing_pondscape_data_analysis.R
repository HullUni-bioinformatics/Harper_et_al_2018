#' ---
#' Title: "Threshold analysis for 'Understanding biodiversity at the pondscape using environmental DNA: a focus great crested newts"
#' Author: "Lynsey Rebecca Harper"
#' Date: "22nd December 2016"
#' ---
#' 
#' 
#' A total of 532 samples previously analysed by qPCR for great crested newt 
#' (GCN) were re-analysed by eDNA metabarcoding to compare method sensitivity 
#' and evaluate this method for detection of vertebrate communities in 
#' freshwater ponds.
#' 
#' Environmental metadata was also obtained for these samples from 
#' Natural England and applied to the species inventories obtained by
#' eDNA metabarcoding.
#' 
#' This analysis aims to test ecological hypotheses using eDNA metabarcoding,
#' with the GCN as a focal species due to extensive literature on its ecology.
#' We aim to identify species associations between GCN and other vertebrates, 
#' biotic and abiotic determinants of GCN pond occupancy, and determinants
#' of overall vertebrate species richness.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. 
#' Then, load functions for plotting model residuals and testing model 
#' fit.
#'

## Clear memory
rm(list=ls())

## Direct R to folder containing packages on workstation
.libPaths(c("C:\\R\\rlib", .libPaths("rlib")))

## set working directory to the location of the script
#install.packages("rstudioapi") # first time for each computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory
getwd()


## Load required packages
p <- c("ggplot2","grid","gridExtra","lme4","glmmADMB","MASS","car","scales",
       "AICcmodavg","gtools","reshape2","plyr","dplyr","arm",
       "RVAideMemoire","ResourceSelection","bbmle","RColorBrewer", 
       "MuMIn","rpart","mgcv","ncf","glmmML","digest","LMERConvenienceFunctions")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies = TRUE)

lapply(p, require, character.only = TRUE)


## Load custom functions
f <- c("CheckResidsFunction.R","OverdispersalFunction.R",
       "CheckConvergenceFunction.R", "HighstatLibV6.R",
       "MyLibrary.R", "glmmML_pres_function.R",
       "ggplotLegendFunction.R", "grid_arrange_shared_legend_Function.R")

lapply(f, source)


#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()


#'
#' ## Metabarcoding dataset
#' 
#' Next, examine the final metabarcoding dataset produced for Chapter 1
#' of my thesis and published in Harper *et al.* (2018) Needle in a haystack?
#' A comparison of eDNA metabarcoding and targeted qPCR for detection of 
#' the great crested newt (*Triturus cristatus*). *Ecology and Evolution*.
#' 

vert.dat <- read.csv("../Data/MetabarcodingData.csv", header=TRUE)

summary(vert.dat)
head(vert.dat)
names(vert.dat)
str(vert.dat)


## Remove PCR positive control, cichlid DNA, as this is no longer 
## required for analysis
## Also remove Human as this cannot be ruled out as laboratory 
## contamination
## Clupea harengus was only present in PCR controls so remove this
## species as well
vert.dat <- vert.dat[,-c(30,44,69)]


#'
#' Examine the effect of other vertebrate species on GCN occupancy.
#' 
#' First, I need to create columns containing the number of species
#' in each vertebrate group found in each eDNA sample.
#' 
#' Convert the read count data to binary presence/absence data
#' 

## Create another data frame for presence/absence detection of species in samples
vert_binary <- vert.dat

## Remove columns created for previous analysis in Harper et al. (2018) Ecol. Evol.
vert_binary <- vert_binary[,-c(2:16,77)]

## Convert sequence read counts to presence-absence
vert_binary[,2:61][vert_binary[,2:61] > 0] <- 1


## total up number of fish species detected in each sample
vert_binary$Fish <- rowSums(vert_binary[,c(3,6:7,12,15:16,19,24:25,
                                           40,45,49,52,56)])

head(vert_binary[,60:62])


## Total up number of amphibian species other than GCN detected in 
## each sample
vert_binary$Amphibian <- rowSums(vert_binary[,c(9,29:30,43,50)])
head(vert_binary[,60:63])


## Total up number of waterfowl species detected in each sample
vert_binary$Waterfowl <- rowSums(vert_binary[,c(2,4,21:22,26)])
head(vert_binary[,60:64])


## Total up number of terrestrial bird species detected in each sample
vert_binary$Woodland.Bird <- rowSums(vert_binary[,c(10,13,17,23,27,32,
                                                    39,44,46,48,54,57:58)])
head(vert_binary[,60:65])


## Total up number of mammal species detected in each sample
vert_binary$Mammal <- rowSums(vert_binary[,c(5,8,11,14,18,20,28,31,33:37,
                                             41:42,47,51,53,55,59,61)])
head(vert_binary[,65:66])



## Plot number of vertebrate species present in ponds with and without GCN

p1 <- ggplot(vert_binary, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p1 <- p1 + geom_bar(stat="count", alpha=.7, position='fill', colour='black')
p1 <- p1 + geom_text(stat='count', aes(label=..count.., ymin=0, ymax=100),
                     size=5, position=position_fill(), vjust=1.5)
p1 <- p1 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p1 <- p1 + scale_fill_manual(name="Crested newt",
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p1 <- p1 + labs(title="(a)", x="Number of other amphibian species", y="Ponds")
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))


p2 <- ggplot(vert_binary, aes(x=factor(Fish), fill=factor(Triturus_cristatus)))
p2 <- p2 + geom_bar(stat="count", alpha=.7, position='fill', colour='black')
p2 <- p2 + geom_text(stat='count', aes(label=..count.., ymin=0, ymax=100),
                     size=5, position=position_fill(), vjust=1.5)
p2 <- p2 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p2 <- p2 + scale_fill_manual(name="Crested newt",
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p2 <- p2 + labs(title="(b)", x="Number of fish species", y="Ponds")
p2 <- p2 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))


p3 <- ggplot(vert_binary, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3 <- p3 + geom_bar(stat="count", alpha=.7, position='fill', colour='black')
p3 <- p3 + geom_text(stat='count', aes(label=..count.., ymin=0, ymax=100),
                     size=5, position=position_fill(), vjust=1.5)
p3 <- p3 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p3 <- p3 + scale_fill_manual(name="Crested newt",
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p3 <- p3 + labs(title="(c)", x="Number of waterfowl species", y="Ponds")
p3 <- p3 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))


p4 <- ggplot(vert_binary, aes(x=factor(Woodland.Bird), fill=factor(Triturus_cristatus)))
p4 <- p4 + geom_bar(stat="count", alpha=.7, position='fill', colour='black')
p4 <- p4 + geom_text(stat='count', aes(label=..count.., ymin=0, ymax=100),
                     size=5, position=position_fill(), vjust=1.5)
p4 <- p4 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p4 <- p4 + scale_fill_manual(name="Crested newt",
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p4 <- p4 + labs(title="(d)", x="Number of terrestrial bird species", y="Ponds")
p4 <- p4 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=20))


p5 <- ggplot(vert_binary, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5 <- p5 + geom_bar(stat="count", alpha=.7, position='fill', colour='black')
p5 <- p5 + geom_text(stat='count', aes(label=..count.., ymin=0, ymax=100),
                     size=5, position=position_fill(), vjust=1.5)
p5 <- p5 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p5 <- p5 + scale_fill_manual(name="GCN",
                             values=c("grey45","orange"),
                             breaks=c("0","1"),
                             labels=c("Negative",
                                      "Positive"))
p5 <- p5 + labs(title="(e)", x="Number of mammal species", y="Ponds")
p5 <- p5 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "bottom",
                 text = element_text(size=20))
p5 <- p5 + guides(fill=guide_legend(reverse=T))


# SHOW ALL PLOTS
# grid.arrange() is a function within the gridEXTRA package that allows to
# arrange how ggplots are viewed, similar to par(mfrow=c())
# Within this function, I can use arrangeGrob() to specify how many plots I
# want on each row
# I want to display my factors on the top row and my continuous predictors
# on the bottom row

grid.arrange(arrangeGrob(p1,p2,p3, nrow=1, ncol=3),
             arrangeGrob(p4,p5, nrow=1, ncol=2), nrow=2)



## Statistically test relationships between number of species and 
## GCN occupancy

## Perform model with binary response variable
speciesGCN <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                    + Waterfowl + Woodland.Bird + Mammal,
                    family = binomial,
                    data = vert_binary)

## Compare to null model
nullGCN <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = vert_binary)

anova(speciesGCN, nullGCN)

## LRT p-value significant and AIC of null model much higher thus number of 
## vertebrate species are appropriate explanatory variables


summary(speciesGCN)
anova(speciesGCN)
drop1(speciesGCN, test = "Chi")  # for binomial/integer y models, 
                                 # statistics for significance of each term
display(speciesGCN)
se.ranef(speciesGCN)   # levels of random factor centred around 0


## Check for model overdispersion 
overdisp.glmer(speciesGCN)

## Overdispersion test (chi-square)

1-pchisq(588.346, df=525)    # model overdispersed


## summary() gives inflated value for model residual deviance so usual
## methods of calculating residual deviance are unreliable for GLMMs

## Use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(speciesGCN)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(speciesGCN, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(speciesGCN)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(speciesGCN)

## Plot the residuals against each x vairable
plot(sresid ~ vert_binary$Fish)
plot(sresid ~ vert_binary$Amphibian)
plot(sresid ~ vert_binary$Waterfowl)
plot(sresid ~ vert_binary$Woodland.Bird)
plot(sresid ~ vert_binary$Mammal)


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(speciesGCN)

## 9.43% explained


#'
#' An additional test of model fit is the Hosmer and Lemeshow Goodness 
#' of Fit Test, which is similar to a chi-square test and indicates the 
#' extent to which the model provides better fit than a null model with 
#' no predictors.
#' 
#' If chi-square goodness of fit is not significant, then the model has 
#' adequate fit. Similarly, if the test is significant, the model does not 
#' adequately fit the data.
#' 

hoslem.test(vert_binary$Triturus_cristatus, fitted(speciesGCN))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## Plot model output
## Obtain predicted values for the full data set: vert_binary
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(speciesGCN, newdata=vert_binary, 
                            se.fit = TRUE, print.matrix=T))

# Create new data set with fitted values and original data for plots 1 and 2
box.dat <- cbind(vert_binary, fit) 


# Plot showing relationship between great crested newt presence and 
# species number of each vertebrate group using predicted values from
# model

p1m <- ggplot(box.dat) 
p1m <- p1m + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.3, width=0.3, cex=1)
p1m <- p1m + scale_colour_gradient(low="gray45", high="orange")
p1m <- p1m + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p1m <- p1m + labs(title="", y="Probability of GCN", x = "Number of other amphibian species") 
p1m <- p1m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))


p2m <- ggplot(box.dat) 
p2m <- p2m + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.3, width=0.3, cex=1)
p2m <- p2m + scale_colour_gradient(low="gray45", high="orange")
p2m <- p2m + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p2m <- p2m + labs(title="", y="Probability of GCN", x = "Number of fish species") 
p2m <- p2m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))


p3m <- ggplot(box.dat) 
p3m <- p3m + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.3, width=0.3, cex=1)
p3m <- p3m + scale_colour_gradient(low="gray45", high="orange")
p3m <- p3m + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3m <- p3m + labs(title="", y="Probability of GCN", x = "Number of waterfowl species") 
p3m <- p3m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))


p4m <- ggplot(box.dat) 
p4m <- p4m + geom_jitter(aes(x=factor(Woodland.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.3, width=0.3, cex=1)
p4m <- p4m + scale_colour_gradient(low="gray45", high="orange")
p4m <- p4m + geom_boxplot(aes(x=factor(Woodland.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4m <- p4m + labs(title="", y="Probability of GCN", x = "Number of terrestrial bird species") 
p4m <- p4m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=20))


p5m <- ggplot(box.dat) 
p5m <- p5m + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.3, width=0.3, cex=1)
p5m <- p5m + scale_colour_gradient(low="gray45", high="orange")
p5m <- p5m + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5m <- p5m + labs(title="", y="Probability of GCN", x = "Number of mammal species") 
p5m <- p5m + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=20),
                   legend.position="none",
                   legend.key.size = unit(1, "cm"))


## extract common legend
mylegend <- g_legend(p5)

## Show all plots with common legend
p1_5 <- grid.arrange(arrangeGrob(p1,p1m,
                                 p2,p2m,
                                 p3,p3m,
                                 p4,p4m,
                                 p5 + theme(legend.position="none"),p5m, 
                                 nrow=5, ncol=2),
                     mylegend, nrow=2,heights=c(10, 1))



#'
#' Now, examine individual species in GCN and non-GCN ponds. 
#' 
#' First, need to create a new dataframe.
#' 

## total up number times an assignment occurred across GCN positive
## ponds
GCNpos <- subset(vert_binary, vert_binary$Triturus_cristatus == 1)
GCNneg <- subset(vert_binary, vert_binary$Triturus_cristatus == 0)


## Create dataframe containing only species in GCN positive ponds
## minus GCN
GCNpos_spp <- data.frame(GCNpos[,c(2:59,61)])

## Calculate number of ponds each species found in
GCNpos_spp <- rbind(GCNpos_spp, colSums(GCNpos_spp))


## Create dataframe containing only species in GCN negative ponds
## minus GCN
GCNneg_spp <- data.frame(GCNneg[,c(2:59,61)])

## Calculate number of ponds each species found in
GCNneg_spp <- rbind(GCNneg_spp, colSums(GCNneg_spp))


## Bind rows containing pond totals from each dataframe
totals <- data.frame(rbind(GCNpos_spp[150,], GCNneg_spp[384,]))
rownames(totals) <- c("P", "N")

## Transpose data frame
totals <- t(totals)

## Melt dataframe to create new factor column containing P and N
ponds_spp <- melt(totals, id=c("P","N"))
colnames(ponds_spp) <- c("Species", "GCN", "Ponds")

 
## To complete the dataframe, need to create another factor column
## indicating what vertebrate groups species belong to

ponds_spp$Group <- factor(ifelse(ponds_spp$Species %in% ponds_spp[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                          ifelse(ponds_spp$Species %in% ponds_spp[c(8,28:29,42,49),1], "Amphibian",
                          ifelse(ponds_spp$Species %in% ponds_spp[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                          ifelse(ponds_spp$Species %in% ponds_spp[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                 "NA")))))


## Remove underscores from species names
ponds_spp$Species <- gsub('_', ' ', ponds_spp$Species)


## Subset data by vertebrate groups
Fish <- subset(ponds_spp, Group == "Fish")
Amphib <- subset(ponds_spp, Group == "Amphibian")
Bird <- subset(ponds_spp, Group == "Bird")
Mammal <- subset(ponds_spp, Group == "Mammal")


## Plot proportion of amphibian species in GCN and non-GCN ponds

p6 <- ggplot(Amphib[order(Amphib$GCN, decreasing=T),],
             aes(x=Species, y=Ponds, fill=factor(GCN, levels=c("N","P"))))
p6 <- p6 + geom_bar(stat="identity", alpha=.7, position="fill", colour="black")
p6 <- p6 + geom_text(data=subset(Amphib, Ponds > 0), 
                     aes(label=Ponds), size=5, position=position_fill(), vjust=1.5)
p6 <- p6 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p6 <- p6 + scale_fill_manual(name="Great crested newt",
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive",
                                      "Negative"))
p6 <- p6 + labs(title="(a)", x="Other amphibian species", y="Ponds")
p6 <- p6 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", angle = 45, vjust = 1, hjust=1, face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p6


## Plot proportion of fish species in GCN and non-GCN ponds

p7 <- ggplot(Fish[order(Fish$GCN, decreasing=T),],
             aes(x=Species, y=Ponds, fill=factor(GCN, levels=c("N","P"))))
p7 <- p7 + geom_bar(stat="identity", alpha=.7, position="fill", colour="black")
p7 <- p7 + geom_text(data=subset(Fish, Ponds > 0), 
                     aes(label=Ponds), size=5, position=position_fill(), vjust=1.5)
p7 <- p7 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p7 <- p7 + scale_fill_manual(name="Great crested newt",
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive",
                                      "Negative"))
p7 <- p7 + labs(title="(b)", x="Fish species", y="")
p7 <- p7 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", angle = 45, vjust = 1, hjust=1, face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p7


## Plot bird species in GCN and non-GCN ponds

p8 <- ggplot(Bird[order(Bird$GCN, decreasing=T),],
             aes(x=Species, y=Ponds, fill=factor(GCN, levels=c("N","P"))))
p8 <- p8 + geom_bar(stat="identity", alpha=.7, position="fill", colour="black")
p8 <- p8 + geom_text(data=subset(Bird, Ponds > 0), 
                     aes(label=Ponds), size=5, position=position_fill(), vjust=1.5)
p8 <- p8 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p8 <- p8 + scale_fill_manual(name="Great crested newt",
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive",
                                      "Negative"))
p8 <- p8 + labs(title="(c)", x="Bird species", y="Ponds")
p8 <- p8 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", angle = 45, vjust = 1, hjust=1, face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p8


## Plot mammal species in GCN and non-GCN ponds

p9 <- ggplot(Mammal[order(Mammal$GCN, decreasing=T),],
             aes(x=Species, y=Ponds, fill=factor(GCN, levels=c("N","P"))))
p9 <- p9 + geom_bar(stat="identity", alpha=.7, position="fill", colour="black")
p9 <- p9 + geom_text(data=subset(Mammal, Ponds > 0), 
                     aes(label=Ponds), size=5, position=position_fill(), vjust=1.5)
p9 <- p9 + scale_y_continuous(labels = percent_format(), expand = c(0,0))
p9 <- p9 + scale_fill_manual(name="Great crested newt",
                             values=c("grey45","orange"),
                             breaks=c("P","N"),
                             labels=c("Positive",
                                      "Negative"))
p9 <- p9 + labs(title="(d)", x="Mammal species", y="")
p9 <- p9 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black", angle = 45, vjust = 1, hjust=1, face="italic"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p9


## Show all plots side by side
grid_arrange_shared_legend(p6,p7,p8,p9, nrow=2, ncol=2)



#'
#' ## Are these interactions significant?
#' 
#' The probabilistic model of species co-occurrence (Veech 2013) 
#' measures co-occurrence in the most straightforward way as the 
#' number of sampling sites where two species co-occur. Observed 
#' co-occurrence can be compared to the expected co-occurrence where 
#' the latter is the product of the two species probability of 
#' occurrence multiplied by the number of sampling sites.
#' 
#' This probabilistic model was developed in R to increase 
#' the availability of the model as an easy-to-use method for conducting 
#' pairwise co-occurrence analyses.
#' 
#' The package calculates the probability of selecting a site (or
#' sample) that has species #1 given that it already has species #2.
#' Subsequently, coocur() uses a hypergeometric approach
#' 

## Load package cooccur for probabalistic species coocurrence analysis 
## in R
library("cooccur")

## Adapt presence-absence dataframe to contain only species, no metadata
ponds <- vert_binary[,c(1:61)]

## Transpose data frame so that species are first column
ponds <- t(ponds)
colnames(ponds) <- ponds[1,]
ponds <- ponds[-1,]
ponds <- data.frame(ponds)
ponds <- tibble::rownames_to_column(ponds, "species")

## correct the species names
ponds$species <- gsub('_', ' ', ponds$species)

## Export new data frame as csv file
write.csv(ponds, "../Data/cooccur_data.csv", row.names=FALSE)


#'
#' The function accepts community data (e.g., species by site matrix 
#' or vice-versa) in the form of a data frame or matrix and returns 
#' a list containing pairwise species co-occurrence results.
#' 
#' The aim of the sample analysis is to determine the degree to 
#' which communities contain species that are positively, negatively, 
#' and randomly associated with one another, investigate the 
#' contribution of individual species to these patterns, and to 
#' quantify the strength of the positive and negative associations 
#' between species pairs.
#' 

## Load dataset which contains presence-absence data of all species
## detected in each pond by metabarcoding

ponds <- read.csv("../Data/cooccur_data.csv", row.names=1)


## Run analysis

cooccur.ponds <- cooccur(mat = ponds,        # the dataframe
                         type = "spp_site",  # data organized with species as rows and sites as column
                         thresh = TRUE,      # summarize the most important species associations
                         spp_names = TRUE)   # use species names in data


## NB: thresh=TRUE removes species that simply do not have sufficient 
## occurrence data.
## i.e. they are expected to share less than one site

## The cooccur() function produces an output object of class cooccur 
## containing all of the results from the co-occurrence analysis
class(cooccur.ponds)

## Calling summary() produces a readout of the total positive, 
## negative, and random species pairs classified by the algorithm.
summary(cooccur.ponds)

## 48 positive pairs
## 17 negative pairs
## 299 random pairs

#'
#' Of 1770 species pair combinations, 1406 pairs (79.44 %) were removed 
#' from the analysis because expected co-occurrence was < 1 and 364 pairs 
#' were analyzed.
#' 
#' There is sufficient statistical power to analyse all pairs as none 
#' were unclassifiable.
#'  

## To obtain the complete set of species pairs analyzed, use prob.table()
prob.table(cooccur.ponds)


## sp1 = Order of species 1 in matrix
## sp2 = Order of species 2 in matric
## sp1_inc = Number of sites (or samples) that have species 1
## sp2_inc = Number of sites that have species 2
## obs_cooccur = Observed number of sites having both species
## prob_cooccur = Probability that both species occur at a site
## exp_cooccur = Expected number of sites having both species
## p_lt = Probability that two species co-occur at a frequency less 
## than observed number of co-occurrence sites if the two species 
## were distributed randomly (independently) of one another
## p_gt = Probability of co-occurrence at a frequency greater than 
## the observed frequency
## sp1_name = supplied name of sp1
## sp2_name = supplied name of sp2


## A list of only significant species combinations can be obtained 
## using the print method
print(cooccur.ponds)


#'
#' For a given species pair, p_lt and p_gt represent the probabilities 
#' that those species could co-occur less than or greater than what 
#' is observed in the data, respectively.
#' 
#' They can be interpreted as p-values, thus indicating significance 
#' levels for negative and positive co-occurrence patterns.
#' 
#' In the ponds data, there was 17 significant negative co-occurrence
#' pattern and 48 significant positive co-occurrence patterns.
#' 

#'
#' Use the plot() method on the results object. This will produce a 
#' visualization of all of the pairwise combinations of species and 
#' their co-occurrence signs (positive or negative) using a ggplot2 
#' heatmap.
#' 
#' The plot trims out any species that do not have any significant 
#' negative or positive associations and orders the remaining species 
#' starting from those with the most negative interactions to those 
#' with the most positive interactions.
#' 

plot(cooccur.ponds)


#'
#' The pair() function can be used to specifiy a specific species, 
#' by name or number, to inspect, by default only significant results 
#' are shown, but if all = TRUE then all results will be shown
#' 

## Inspect great crested newt
pair(mod = cooccur.ponds, "Triturus cristatus")
pair(mod = cooccur.ponds, "Triturus cristatus", all = TRUE)

## Borderline significant negative cooccurrences include:
## Triturus cristatus - Anas carolinensis: p_lt = 0.09869
## Triturus cristatus - Cyprinus carpio: p_lt = 0.07035
## Triturus cristatus - Meles meles: p_lt = 0.09869
## Triturus cristatus - Cow: p_gt = 0.09714


#'
#' To understand each species individual contribution to the 
#' positive and negative species associations, need to create a 
#' pairing profile.
#' 
#' pair.attributes() produces a table of the percentage of each 
#' species total pairings that were classified as positive, negative,
#' and random (columns with prefix num are counts).
#' 

pair.attributes(cooccur.ponds)


#'
#' The primary goal is to compare the numbers of positive versus 
#' negative associations (unclassifiable pairings treated as random).
#' Same results can be visualized across all species by using the
#' function pair.profile() to create a box plot of these percentages.
#' This plot will show the percent of species pairs that were positive, 
#' negative, and random for all species.
#' 

p10 <- pair.profile(cooccur.ponds)
p10


#'
#' This communicates whether or not species tend to have mostly negative
#' or mostly positive interactions and suggests whether interactions 
#' are evenly distributed among species as opposed to being clustered
#' in a few species.
#' 


#'
#' Effect sizes can also be calculated from co-occurrence analyses; 
#' they allow for comparisons among studies and methods as well as 
#' providing a quantitative measurement of co-occurrence for use in 
#' downstream analyses.
#' 
#' Effect size = difference between expected and observed frequency 
#' of co-occurrence
#' 
#' Values are standardized by dividing these differences by the 
#' number of sampling sites in the dataset. In standardized form, 
#' these values are bounded from -1 to 1, with positive values 
#' indicating positive associations and negative values indicating 
#' negative associations.
#' 

effects.table <- cooccur(mat = ponds, type = "spp_site", thresh = FALSE,
                         spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, 
                         eff_matrix = TRUE)


## Convert dist object produced by cooccur() to a matrix
effects.table <- as.matrix(dist(effects.table)) 


## Write matrix of effect sizes as csv file to working directory
write.csv(effects.table, "../Data/EffectSize.csv")


#'
#' Cohen provided rules of thumb for interpreting these effect sizes, 
#' suggesting that an r of .1 represents a 'small' effect size, .3 
#' represents a 'medium' effect size and .5 represents a 'large' 
#' effect size. 
#' 
#' To inspect the degree to pond species pairs deviate from their 
#' expected co-occurrence levels, plot the observed values against 
#' the expected value as a visual diagnostic. 
#' 

p11 <- obs.v.exp(cooccur.ponds)
p11


#' 
#' The probability calculations are based on the number of sites and 
#' the individual frequencies of occurrence and co-occurrence for 
#' each species pair. Therefore the conditions determining 
#' statistical power change with sample size and it is valuable to 
#' examine effect sizes for species pairs regardless of statistical 
#' significance. A detailed discussion of power, Type I and II error 
#' rates, and a comparison with other methods can be found in Veech 
#' (2013).
#' 
#' 
#' There are several species pairs that exhibit fewer and greater 
#' expected co-occurrences. These pairs are largely clustered towards 
#' having low or high expected co-occurrences in the first place.
#' 

## Plot pairwise profile and observed-expected plot alongside each 
## other

grid.arrange(p10,p11, ncol = 2,
             top = textGrob("Vertebrate cooccurrence at ponds in the UK", 
                            gp = gpar(cex = 2), 
                            vjust = 0.75))



#' ---
#' 
#' ## Environmental predictors
#' 

## Adapt vert_binary dataframe
## Remove columns created for previous analysis
vert_spp <- vert_binary[,1:61]

## Create column containing the total number of species detected in
## each sample i.e. species richness
vert_spp$Sp_richness <- rowSums(vert_spp[,2:61])

## Import habitat metadata for eDNA samples
metadata <- read.csv("../Data/HabitatMetadata.csv", header=TRUE)


#'
#' I want to merge the habitat metadata with the metabarcoding data 
#' using the 'sample' column from the metadata dataframe and the 
#' 'sample' column from the metabarcoding dataframe.
#'

habitat.dat <- merge(metadata, vert_spp, all=TRUE)

## Check merge was successful
head(habitat.dat[,1:6])
head(habitat.dat[,90:97])

## Remove samples which did not have associated pond metadata
## (ADAS UK Ltd samples and 4 samples from FERA)
habitat.dat <- habitat.dat[-c(505:532),]

## Now, check dataframe structure
summary(habitat.dat)
names(habitat.dat)
head(habitat.dat)
str(habitat.dat)


## Exploratory plots for GCN presence/absence

plot(habitat.dat$Triturus_cristatus ~ habitat.dat$County)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Sp_richness)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Max_Depth)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Approx_Circ)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Max_Width)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Max_Length)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Pond_Area)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Permanance)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Quality)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Pond_Substrate)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Inflow)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Outflow)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Pollution)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Overhang)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Shading)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Amphibians)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Waterfowl)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Fish)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Pond_Density)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Woodland)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Rough_Grass)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Scrub_Hedge)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Ruderals)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Terrestrial_Other)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Overall_Habitat_Score)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Macrophytes)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$HSI)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$HSI_Band)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Bufo_bufo)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Lissotriton_vulgaris)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Phasianus_colchicus)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Sciurus_carolinensis)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Pungitius_pungitius)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Gasterosteus_aculeatus)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Fulica_atra)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Sus_scrofa)
plot(habitat.dat$Triturus_cristatus ~ habitat.dat$Gallinula_chloropus)


## corvif function and pairplot with the Pearson correlation 
## coefficients were imported earlier

## Matrix of Spearman's rank correlations for dataset
## Use Spearman's rank as no assumptions are made about linearity in a
## relationship between two variables
cor(habitat.dat[, 10:19], method = "spearman")

## Booth et al. (1994) and Zuur et al. (2009) suggest that correlations
## between pairs of variables with magnitudes ± 0.5 indicate high
## collinearity

## View pairwise plots of variables
pairs(habitat.dat[, 10:19])

## Below produces a pairwise plot but the strength of colinearity 
## strength is indicated by the size of the text
plot(habitat.dat[, 10:19], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)


#' 
#' There appears to be collinearity between pond circumference, pond
#' length, pond width, and pond area.
#' 
#' Pond area encompasses length and width thus accounting for the same variance 
#' in the data as these variables. Therefore, it is logical to remove pond 
#' circumference, pond length and pond width to counter high collinearity between 
#' explanatory variables.
#' 
#' Shading (% total pond margin shaded) and terrestrial overhang (% pond
#' overhung by trees and shrubs) also appear to be collinear.
#' As terrestrial overhang encompasses the entire pond, wheras shading
#' considers only the pond margin, terrestrial overhang will be retained
#' as an explanatory variable. Shading is a know driver of biodiversity
#' in ponds (Sayer et al. 2012).
#' 
#' Similarly, there is no need to include overall habitat score and HSI
#' band as these are encompassed by HSI score.
#' 
#' 
#' Now check the variance inflation factors (VIFs) of variables to 
#' assess the extent of remaining collinearity.
#' 
#' This is done using a binomial GLM and logit link function containing
#' all explanatory variables to the GCN presence/absence data, then 
#' calculating the VIFs for each variable from the resulting model.
#' 

glmGCN <- glm(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes + HSI 
              + Bufo_bufo + Lissotriton_vulgaris 
              + Phasianus_colchicus + Sciurus_carolinensis 
              + Pungitius_pungitius + Gasterosteus_aculeatus 
              + Fulica_atra + Sus_scrofa + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow 
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score + HSI_Band,
              family = binomial,
              data = habitat.dat)

## Use vif() function in car package to calculate VIFs
vif(glmGCN)
  
            
#'
#' Many VIF values are below 10 but HSI score and HSI band are extremely
#' high. This may be because many of the explanatory variables included 
#' in the model are used to calculate HSI score and thus explain the 
#' same variation. However, given HSI score is highly correlated with 
#' great crested newt presence (Oldham *et al.* 2000), it is important to 
#' test this relationship. This will be tested separately after the 
#' final model has been achieved.
#' 
#' Remove HSI/HSI band and recalculate VIFs.
#' 

glmGCN <- glm(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris 
              + Phasianus_colchicus + Sciurus_carolinensis 
              + Pungitius_pungitius + Gasterosteus_aculeatus 
              + Fulica_atra + Sus_scrofa + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow 
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score,
              family = binomial,
              data = habitat.dat)

vif(glmGCN)

## All VIF values are now below 5.


## Check possible relationships between response variable and nominal
## variables (factors) in a design plot.

factors <- data.frame(Triturus_cristatus = habitat.dat$Triturus_cristatus,
                      Permanance = habitat.dat$Permanance,
                      Quality = habitat.dat$Quality,
                      Pond_Substrate = habitat.dat$Pond_Substrate,
                      Inflow = habitat.dat$Inflow,
                      Outflow = habitat.dat$Outflow,
                      Pollution = habitat.dat$Pollution,
                      Amphibians = habitat.dat$Amphibians,
                      Waterfowl = habitat.dat$Waterfowl,
                      Fish = habitat.dat$Fish,
                      Woodland = habitat.dat$Woodland,
                      Rough_Grass = habitat.dat$Rough_Grass,
                      Scrub_Hedge = habitat.dat$Scrub_Hedge,
                      Ruderals = habitat.dat$Ruderals,
                      Terrestrial_Other = habitat.dat$Terrestrial_Other,
                      Overall_Habitat_Score = habitat.dat$Overall_Habitat_Score)

plot.design(Triturus_cristatus ~ Permanance + Quality + Pond_Substrate 
            + Inflow + Outflow + Pollution + Amphibians + Waterfowl 
            + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
            + Terrestrial_Other + Overall_Habitat_Score, 
            data = factors, axes = T, xtick = T)


#' 
#' The highest mean values were observed with pond substrate,
#' amphibian presence and permanance.
#' 

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris 
              + Phasianus_colchicus + Sciurus_carolinensis 
              + Pungitius_pungitius + Gasterosteus_aculeatus 
              + Fulica_atra + Sus_scrofa + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)

GCN_tree <- rpart(f1, data = habitat.dat, method = "class",
                  minsplit=5, cp = 0.001)

par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - L. vulgaris presence
#' - species richness
#' - max depth
#' - fish presence
#' - pond density
#' - toad presence
#' - pond area
#' - amphibian presence
#' - waterfowl presence
#' - other good terrestrial habitat
#' - pond substrate
#' - overhang
#' - overall habitat score
#' - grey squirrel presence
#' - scrub/hedge
#' - three-spined stickleback presence
#' - outflow
#' - ruderals
#' - macrophytes
#' - water quality
#' - pig presence
#' - rough grass
#' - pond permanance
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 3 or 27 is optimal i.e. no. explanatory variables. 
#' The best tree will be subjective.
#' 
#' L. vulgaris and T. cristatus known to have shared ecology including
#' foraging and habitat preference.
#' 
#' Species richness may provide indication of pond health and capability
#' to support more species.
#' 
#' Depth may play a role in GCN occupancy as more fish species are 
#' likely to inhabit deep ponds and less breeding substrate available
#' for GCN.
#' 
#' Large predatory fish species feed on newt larvae.
#' 
#' Pond density can influence migration and population viability.
#' 
#' Larger ponds are likely to have more fish species but may also 
#' support greater GCN populations.
#' 
#' Presence of other amphibians could play positive or negative role
#' in GCN occupancy as prey availability will be increased but 
#' predation and competition for resources will also increase.
#' 
#' Similarly, some fish species predate GCN and GCN are in competition
#' with others for food.
#' 
#' Waterfowl species may have shared habitat/prey preference as GCN.
#' 
#' GCN are known to spend majority of their life on land so terrestrial
#' habitat will greatly influence pond occupancy e.g. overhang of 
#' ponds and grassland.
#' 
#' Pond substrate will influence macroinvertebrate prey availability.
#' 
#' Overhang is an indication of terrestrial plant life encroaching 
#' on pond.
#' 
#' Great crested newt are depending upon rough grass habitat for terrestrial
#' migration and scrub/hedge and ruderal habitat for hibernacula.
#' 
#' Grey squirrel inhabit wooded areas which may not be suitable for 
#' GCN.
#' 
#' Outflow may indicate whether pollutants and excess nutriets are 
#' flushed out.
#' 
#' Macrophyte cover is a known driver of biodiversity in ponds and
#' influence GCN occupancy.
#' 
#' Water quality is an indication of pollution and runoff
#' 
#' Permanant ponds will be more suitable for GCN breeding but may
#' contain more fish species. 
#' 
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - pollution
#' - inflow
#' - woodland
#' 
#' 
#' Therefore, variables to be taken forward in model selection are:
#' 
#' - smooth newt presence
#' - species richness
#' - max depth
#' - fish presence
#' - pond density
#' - toad presence
#' - pond area
#' - amphibian presence
#' - waterfowl presence
#' - other good terrestrial habitat
#' - pond substrate
#' - overhang
#' - overall terrestrial habitat quality
#' - grey squirrel presence
#' - scrub/hedge
#' - three-spined stickleback presence
#' - outflow
#' - ruderals
#' - macrophytes
#' - water quality
#' - pig presence
#' - rough grass
#' - pond permanance
#' 


#'
#' Several variables occur more than once in the tree which indicates 
#' weak non-linear relationships with the response variable thus a 
#' GAM may be more appropriate to model the data.
#' 
#' A GAM with binomial distribution and logistic link function will be
#' used to relate GCN presence/absence with the explanatory variables.
#' It is actually a Bernoulli distribution thus overdispersion cannot 
#' occur which is a major problem in ecological distribution studies.
#' 
#' Forward selection will be applied to find the optimal set of 
#' explanatory variables based on model AIC values.

## Use package mgcv for gam modelling

GCN.gam1 <- gam(Triturus_cristatus ~ Lissotriton_vulgaris, data=habitat.dat, family=binomial)
GCN.gam2 <- gam(Triturus_cristatus ~ s(Sp_richness), data=habitat.dat, family=binomial)
GCN.gam3 <- gam(Triturus_cristatus ~ s(Max_Depth), data=habitat.dat, family=binomial)
GCN.gam4 <- gam(Triturus_cristatus ~ Fish, data=habitat.dat, family=binomial)
GCN.gam5 <- gam(Triturus_cristatus ~ s(Pond_Density), data=habitat.dat, family=binomial)
GCN.gam6 <- gam(Triturus_cristatus ~ Bufo_bufo, data=habitat.dat, family=binomial)
GCN.gam7 <- gam(Triturus_cristatus ~ s(Pond_Area), data=habitat.dat, family=binomial)
GCN.gam8 <- gam(Triturus_cristatus ~ Amphibians, data=habitat.dat, family=binomial)
GCN.gam9 <- gam(Triturus_cristatus ~ Waterfowl, data=habitat.dat, family=binomial)
GCN.gam10 <- gam(Triturus_cristatus ~ Pond_Substrate, data=habitat.dat, family=binomial)
GCN.gam11 <- gam(Triturus_cristatus ~ s(Overhang), data=habitat.dat, family=binomial)
GCN.gam12 <- gam(Triturus_cristatus ~ Sciurus_carolinensis, data=habitat.dat, family=binomial)
GCN.gam13 <- gam(Triturus_cristatus ~ Scrub_Hedge, data=habitat.dat, family=binomial)
GCN.gam14 <- gam(Triturus_cristatus ~ Gasterosteus_aculeatus, data=habitat.dat, family=binomial)
GCN.gam15 <- gam(Triturus_cristatus ~ Outflow, data=habitat.dat, family=binomial)
GCN.gam16 <- gam(Triturus_cristatus ~ Ruderals, data=habitat.dat, family=binomial)
GCN.gam17 <- gam(Triturus_cristatus ~ s(Macrophytes), data=habitat.dat, family=binomial)
GCN.gam18 <- gam(Triturus_cristatus ~ Rough_Grass, data=habitat.dat, family=binomial)
GCN.gam19 <- gam(Triturus_cristatus ~ Quality, data=habitat.dat, family=binomial)
GCN.gam20 <- gam(Triturus_cristatus ~ Permanance, data=habitat.dat, family=binomial)
GCN.gam21 <- gam(Triturus_cristatus ~ Terrestrial_Other, data=habitat.dat, family=binomial)
GCN.gam22 <- gam(Triturus_cristatus ~ Overall_Habitat_Score, data=habitat.dat, family=binomial)
GCN.gam23 <- gam(Triturus_cristatus ~ Sus_scrofa, data=habitat.dat, family=binomial)


AIC(GCN.gam1,GCN.gam2,GCN.gam3,GCN.gam4,GCN.gam5,GCN.gam6, GCN.gam7,
    GCN.gam8,GCN.gam9,GCN.gam10,GCN.gam11,GCN.gam12,GCN.gam13,
    GCN.gam14,GCN.gam15,GCN.gam16,GCN.gam17,GCN.gam18,GCN.gam19,
    GCN.gam20,GCN.gam21,GCN.gam22,GCN.gam23)


## Species richness identified as best explanatory variable and most 
## optimal model

plot(GCN.gam2)

summary(GCN.gam1)
summary(GCN.gam2)
summary(GCN.gam3)
summary(GCN.gam4)
summary(GCN.gam5)
summary(GCN.gam6)
summary(GCN.gam7)
summary(GCN.gam8)
summary(GCN.gam9)
summary(GCN.gam10)
summary(GCN.gam11)
summary(GCN.gam12)
summary(GCN.gam13)
summary(GCN.gam14)
summary(GCN.gam15)
summary(GCN.gam16)
summary(GCN.gam17)
summary(GCN.gam18)
summary(GCN.gam19)
summary(GCN.gam20)
summary(GCN.gam21)
summary(GCN.gam22)
summary(GCN.gam23)


#' 
#' However, the other variables are biologically important for GCN and
#' several have an estimated 1 degree of freedom for the smoothers
#' which is equivalent to a linear relationship.
#' 
#' Therefore, GLM should be applied instead of GAM. GLM is also a 
#' parametric method. 
#' 
#' However, spatial autocorrrelation is common in ecological studies 
#' of species presence/absence as sites that are located within an 
#' animal's ranging capability are likely to be inhabited.
#' 
#' In the case of GCN, individuals can migrate distances of 2km when 
#' they have finished breeding. Therefore, occurrence of GCN is likely
#' in ponds that are closely located to one another in a given area.
#' 
#' Furthermore, if pond distances are substantially smaller than the 
#' spatial extent of the study area, this can also lead to spatial 
#' auto-correlation.
#' 
#' The underlying spatial patterns of habitat surrounding ponds could
#' also be a source of spatial autocorrelation but the explanatory 
#' variables which describe the habitat should account for this.
#' 
#' One way to assess spatial autocorrelation is through correlograms 
#' of the data. These are graphical representations of spatial 
#' correlation between locations at a range of lag distances. Positive 
#' spatial correlation indicates spatial autocorrelation between data 
#' points. Negative spatial correlation may also indicate problems but 
#' is very uncommon.
#' 
#' Spline correlograms are smoothed correlograms using a spline 
#' function (Zuur *et al.* 2009).
#' 

## The ncf package will be used to produce spline correlograms
## Location data is required 

Correlog <- spline.correlog(x=habitat.dat$Easting,
                            y=habitat.dat$Northing,
                            z=habitat.dat$Triturus_cristatus, 
                            xmax=10000)

plot.spline.correlog(Correlog)


#' 
#' A spline correlogram with 95% pointwise bootstrap CIs and maximum 
#' lag distance of 10km has been produced.
#' 
#' Some positive spatial auto-correlation is present but only
#' at short lag distances of less than or around 1km. This suggests
#' that spatial auto-correlation may be a problem for sites that are
#' located near to one another.
#' 
#' Check for spatial auto-correlation in Pearson residuals of logistic
#' regression model containing all explanatory variables fitted to the
#' presence/absence data.
#' 

glmGCN <- glm(Triturus_cristatus ~ Lissotriton_vulgaris + Sp_richness 
              + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
              + Amphibians + Waterfowl + Pond_Substrate + Overhang
              + Sciurus_carolinensis + Scrub_Hedge + Gasterosteus_aculeatus
              + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
              + Permanance,
              family = binomial,
              data = habitat.dat)


Correlog_glmGCN <- spline.correlog(x=habitat.dat$Easting,
                                   y=habitat.dat$Northing,
                                   z=residuals(glmGCN, type="pearson"),
                                   xmax=10000)

plot.spline.correlog(Correlog_glmGCN, ylim=c(-5, 5))


#' 
#' A GLMM can account for dependencies within sites and is appropriate
#' where nested model selection will be used and the data are structured 
#' hierarchically. Ponds are nested within counties thus a mixed model
#' is necessary to account for spatial dependencies within sites. Each
#' sample represents a different pond, or a different site, and thus
#' sample must be treated as a random effect.
#' 
#' 
#' Apply GLMM to variables most likely to influence GCN presence/absence
#' based on exploratory analysis
#' 
#' 
#' Sample (individual pond) will be random factor as interested in GCN 
#' presence/absence across all ponds as a whole not in specific ponds.
#' 
#' Check that a GLMM will adequately account for spatial auto-correlation
#' that is present by fitting a model including all explanatory 
#' variables and examining a spline correlogram of the Pearson 
#' residuals.
#' 
#' glmmML estimates model parameters by maximum likelihood and allows 
#' AIC values to be calculated. However, lmer in lme4 can also be 
#' used.
#' 
#' The cluster argument indicates the grouping level for the random 
#' effect.
#' 

modelML <- glmmML(formula=Triturus_cristatus ~ Lissotriton_vulgaris 
                  + Sp_richness + Max_Depth + Fish + Pond_Density 
                  + Bufo_bufo + Pond_Area + Amphibians + Waterfowl 
                  + Pond_Substrate + Overhang + Sciurus_carolinensis 
                  + Scrub_Hedge + Gasterosteus_aculeatus + Outflow 
                  + Ruderals + Macrophytes + Rough_Grass + Quality
                  + Permanance + Sus_Scrofa + Overall_Habitat_Score,
                  cluster = sample,
                  family = binomial,
                  data = habitat.dat)


## Check for spatial autocorrelation of the Pearson residuals

Correlog.modelML <- spline.correlog(x=habitat.dat$Easting,
                                    y=habitat.dat$Northing,
                                    z=pres.glmmML(model = modelML,
                                                  data = habitat.dat), 
                                    na.rm=TRUE, 
                                    xmax=10000)

plot.spline.correlog(Correlog.modelML)


#' 
#' Spatial autocorrelation has been removed at short lag distances
#' thus both models are successfully accomodating some of the 
#' spatial autocorrelation within sites. However, there is no real
#' difference between a GLM and GLMM therefore mixed effect modeling
#' is not necessarily required to remove spatial autocorrelation. 
#' To examine GCN occupancy at the pondscape instead of individual 
#' ponds though, I have to model sample as a random effect.
#' 
#' Now a suitable set of explanatory variables and modelling framework
#' has been identified, I now need to identify which variables are 
#' important determinants of GCN distribution and choose a suitable,
#' parsimonious approximating model to make predictions. This is
#' done using an information-theoretic approach using AIC criteria
#' to evaluate model fit. Low AIC models are more parsimonious than
#' high AIC models.
#' 
#' 
#' First, group explanatory variables into different functional groups:
#' 
#' - Pond properties (pond substrate, max depth, pond area, 
#'   water quality, outflow, pond density, permanence, macrophytes)
#'   
#' - Terrestrial habitat (terrestrial overhang, rough grass, scrub/hedge,
#'   ruderals, overall terrestrial habitat quality)
#'   
#' - Species (smooth newt, toad, amphibian, fish, waterfowl, 
#'   species richness, grey squirrel, three-spined stickleback, pig)
#'   


## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist <- data.frame(habitat.dat[,c(1,10,14:16,18,20:22,24,26:28,30:34,
                                     44,59,65,88,94:95,97)])


## Model with all terms: 

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo 
                + Pond_Area + Amphibians + Waterfowl + Pond_Substrate 
                + Overhang + Sciurus_carolinensis + Scrub_Hedge 
                + Gasterosteus_aculeatus + Outflow + Ruderals 
                + Macrophytes + Rough_Grass + Quality
                + Permanance + Terrestrial_Other,
                family = binomial,
                data = GCNdist)

summary(model0)   # AIC = 530.2
drop1(model0, test = "Chi")


## Pond properties:

pondprop <- glmer(Triturus_cristatus ~ (1|sample) 
                  + Max_Depth + Pond_Density 
                  + Pond_Area + Pond_Substrate + Outflow
                  + Macrophytes + Quality + Permanance,
                  family = binomial,
                  data = GCNdist)

summary(pondprop)  # AIC = 604.3, not as good as first model
drop1(pondprop, test="Chi")


## Terrestrial habitat:

terr <- glmer(Triturus_cristatus ~ (1|sample) 
              + Overhang + Rough_Grass + Ruderals + Scrub_Hedge
              + Terrestrial_Other + Overall_Habitat_Score,
              family = binomial,
              data = GCNdist)

summary(terr)   # AIC = 599.4, again not as good as first model
drop1(terr, test = "Chi")


## Species

species <- glmer(Triturus_cristatus ~ (1|sample) + Lissotriton_vulgaris 
                 + Sp_richness + Amphibians + Fish + Waterfowl
                 + Sciurus_carolinensis + Gasterosteus_aculeatus
                 + Bufo_bufo + Sus_scrofa,
                 family = binomial,
                 data = GCNdist)

summary(species)   # AIC = 498.1 but this does not provide information on habitat
drop1(species, test = "Chi")



#' 
#' The AIC comparison of alternate models implies that the explanatory
#' variables best explain the GCN presence/absence data when they are
#' combined in a single model. However, no estimation of standard 
#' error implies overparameterisation.
#' 
#' Apply step-wise down nested model selection using glmer instead.
#' 

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Amphibians + Waterfowl + Pond_Substrate + Overhang
                + Sciurus_carolinensis + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Permanance + Terrestrial_Other,
                family = binomial,
                data = GCNdist)

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Amphibians + Waterfowl + Pond_Substrate + Overhang
                + Sciurus_carolinensis + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Permanance + Terrestrial_Other + Overall_Habitat_Score
                + Sus_scrofa,
                family = binomial,
                data = GCNdist)


#'
#' Use drop1() function to select order that variables should be 
#' dropped from the model. Function computes all the single terms 
#' in the scope argument that can be added to or dropped from the 
#' model, fits those models and computes a table of the changes in 
#' fit.
#' 

drop1(model0)


## Remove pond substrate as this variable has the lowest AIC value

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Amphibians + Waterfowl + Overhang
                + Sciurus_carolinensis + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Permanance + Terrestrial_Other,
                family = binomial,
                data = GCNdist)

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Amphibians + Waterfowl + Overhang
                + Sciurus_carolinensis + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Permanance + Terrestrial_Other + Overall_Habitat_Score
                + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model1, model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit

drop1(model1)

## Remove amphibian presence from model next as this has the lowest AIC
## value

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Permanance + Terrestrial_Other,
                family = binomial,
                data = GCNdist)

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Permanance + Terrestrial_Other + Overall_Habitat_Score
                + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model2, model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibian presence has resulted in better model fit.

drop1(model2)

## Remove permanance from model next as this has the lowest AIC 
## value

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Terrestrial_Other,
                family = binomial,
                data = GCNdist)

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass + Quality
                + Terrestrial_Other + Overall_Habitat_Score
                + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model3, model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond permanance has resulted in better model fit

drop1(model3)

## Remove water quality from model next as this has the lowest AIC 
## value

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes + Rough_Grass
                + Terrestrial_Other + Overall_Habitat_Score
                + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model4, model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of water quality has resulted in better model fit

drop1(model4)

## Remove rough grass from model next as this has the lowest AIC 
## value

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Fish + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes
                + Terrestrial_Other + Overall_Habitat_Score
                + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model5, model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of rough grass has resulted in better model fit

drop1(model5)

## Remove fish from model next as this has the lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes
                + Terrestrial_Other + Overall_Habitat_Score
                + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model6, model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish has resulted in better model fit

drop1(model6)

## Remove overall terrestrial habitat quality from model next 
## as this has the lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Scrub_Hedge + Gasterosteus_aculeatus
                + Outflow + Ruderals + Macrophytes
                + Terrestrial_Other + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model7, model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of overall terrestial habitat has resulted in better 
## model fit

drop1(model7)

## Remove scrub/hedge from model next as this has the lowest AIC 
## value

model8 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Overhang + Sciurus_carolinensis 
                + Gasterosteus_aculeatus + Outflow + Ruderals 
                + Macrophytes + Terrestrial_Other + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model8, model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit

drop1(model8)

## Remove terrestrial overhang from model next as this has the lowest 
## AIC value

model9 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Pond_Density + Bufo_bufo + Pond_Area
                + Waterfowl + Sciurus_carolinensis 
                + Gasterosteus_aculeatus + Outflow + Ruderals 
                + Macrophytes + Terrestrial_Other + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model9, model8)

## p-value from LRT not significant and reduction in AIC value
## so removal of terrestrial overhang has resulted in better model fit

drop1(model9)

## Remove pond density from model next as this has the lowest 
## AIC value

model10 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Max_Depth + Bufo_bufo + Pond_Area
                + Waterfowl + Sciurus_carolinensis 
                + Gasterosteus_aculeatus + Outflow + Ruderals 
                + Macrophytes + Terrestrial_Other + Sus_scrofa,
                family = binomial,
                data = GCNdist)

anova(model10, model9)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond density has resulted in better model fit

drop1(model10)

## Remove macrophytes from model next as this has the lowest 
## AIC value

model11 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Max_Depth + Bufo_bufo + Pond_Area
                 + Waterfowl + Sciurus_carolinensis 
                 + Gasterosteus_aculeatus + Outflow + Ruderals 
                 + Terrestrial_Other + Sus_scrofa,
                 family = binomial,
                 data = GCNdist)

anova(model11, model10)

## p-value from LRT not significant and reduction in AIC value.
## Reduction is <2 but better to remove variable to avoid 
## overparameterisation. Removal of macrophytes has resulted in better 
## model fit.

drop1(model11)

## Remove waterfowl presence from model next as this has the lowest 
## AIC value

model12 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Max_Depth + Bufo_bufo + Pond_Area
                 + Sciurus_carolinensis + Gasterosteus_aculeatus 
                 + Outflow + Ruderals + Terrestrial_Other 
                 + Sus_scrofa,
                 family = binomial,
                 data = GCNdist)

anova(model12, model11)

## p-value from LRT not significant but increase in AIC value.
## However, reduction is <2 so better to remove variable to avoid 
## overparameterisation. Removal of waterfowl presence has resulted 
## in better model fit.

drop1(model12)

## Remove pig presence from model next as this has the lowest 
## AIC value

model13 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Max_Depth + Bufo_bufo + Pond_Area
                 + Sciurus_carolinensis + Gasterosteus_aculeatus 
                 + Outflow + Ruderals + Terrestrial_Other,
                 family = binomial,
                 data = GCNdist)

anova(model13, model12)


## p-value from LRT not significant minimal change in AIC value.
## Reduction in AIC is <2 so better to remove variable to avoid 
## overparameterisation. Removal of pig presence has resulted 
## in better model fit.

drop1(model13)

## Remove pond depth from model next as this has the lowest 
## AIC value

model14 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Bufo_bufo + Pond_Area
                 + Sciurus_carolinensis + Gasterosteus_aculeatus 
                 + Outflow + Ruderals + Terrestrial_Other,
                 family = binomial,
                 data = GCNdist)

anova(model14, model13)

## p-value from LRT significant and increase in AIC value (>2)
## Removal of pond depth has reduced model fit.

drop1(model13)

## Remove outflow to see if final model reached as this has 
## lowest AIC value after pond depth

model15 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Max_Depth + Bufo_bufo + Pond_Area
                 + Sciurus_carolinensis + Gasterosteus_aculeatus 
                 + Ruderals + Terrestrial_Other,
                 family = binomial,
                 data = GCNdist)

anova(model15, model13)

## p-value from LRT significant and increase in AIC value (>2)
## Removal of outflow has reduced model fit.


#'
#' Variation in GCN occupancy is best explained by species richness,
#' pond area, pond depth, ruderals, outflow, good terrestrial habitat,
#' smooth newt presence, common toad presence, grey squirrel presence, 
#' and three-spined stickleback presence.
#' 


# Final model: 

GCNmodel <- glmer(Triturus_cristatus ~ (1|sample)
                  + Lissotriton_vulgaris + Bufo_bufo
                  + Sciurus_carolinensis + Gasterosteus_aculeatus
                  + Sp_richness + Max_Depth + Pond_Area + Outflow 
                  + Ruderals + Terrestrial_Other,
                  family = binomial,
                  data = GCNdist)

summary(GCNmodel)
anova(GCNmodel)
drop1(GCNmodel, test = "Chi")   # for binomial/integer y models, 
                                # statistics for significance of each term

display(GCNmodel)
se.ranef(GCNmodel)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(GCNmodel)

## Overdispersion test (chi-square)

1-pchisq(424.838, df=490)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(GCNmodel)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(GCNmodel, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(GCNmodel)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(GCNmodel)


## Plot the residuals against each x vairable
plot(sresid ~ GCNdist$Lissotriton_vulgaris)
plot(sresid ~ GCNdist$Sp_richness)
plot(sresid ~ GCNdist$Max_Depth)
plot(sresid ~ GCNdist$Bufo_bufo)
plot(sresid ~ GCNdist$Pond_Area)
plot(sresid ~ GCNdist$Sciurus_carolinensis)
plot(sresid ~ GCNdist$Gasterosteus_aculeatus)
plot(sresid ~ GCNdist$Outflow)
plot(sresid ~ GCNdist$Ruderals)
plot(sresid ~ GCNdist$Terrestrial_Other)


## Check model residuals
chkres(GCNmodel)

## Plot the fitted data against the observed data
plot(GCNdist$Triturus_cristatus ~ fitted(GCNmodel))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(GCNmodel)

## 38.59% explained


#'
#' An additional test of model fit is the Hosmer and Lemeshow Goodness 
#' of Fit Test, which is similar to a chi-square test and indicates the 
#' extent to which the model provides better fit than a null model with 
#' no predictors.
#' 
#' If chi-square goodness of fit is not significant, then the model has 
#' adequate fit. Similarly, if the test is significant, the model does not 
#' adequately fit the data.
#' 

hoslem.test(GCNdist$Triturus_cristatus, fitted(GCNmodel))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data


## Examine model coefficients
library(coefplot)
coefplot(GCNmodel)



## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(GCNmodel, newdata=GCNdist, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(GCNdist, fit) 



## PLOT 12: Smooth newt presence and GCN presence/absence 

p12 <- ggplot(box.dat) + ggtitle('(a)')
p12 <- p12 + geom_jitter(aes(x=factor(Lissotriton_vulgaris), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.3, height=0.3, cex=1)
p12 <- p12 + geom_boxplot(aes(x=factor(Lissotriton_vulgaris), y=fit), alpha=0.7, outlier.colour = "blue")
p12 <- p12 + scale_colour_gradient(name="Great crested newt occupancy",
                                   low="grey45", high="orange")
p12 <- p12 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p12 <- p12 + labs(x = "Smooth newt",
                  y = "Probability of GCN") 
p12 <- p12 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold", hjust=0),
                 text = element_text(size=24),
                 legend.position = "none")
p12



## PLOT 13: Common toad presence and GCN presence/absence 

p13 <- ggplot(box.dat) + ggtitle('(b)')
p13 <- p13 + geom_jitter(aes(x=factor(Bufo_bufo), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.3, height=0.3, cex=1)
p13 <- p13 + geom_boxplot(aes(x=factor(Bufo_bufo), y=fit), alpha=0.7, outlier.colour = "blue")
p13 <- p13 + scale_colour_gradient(low="grey45", high="orange")
p13 <- p13 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p13 <- p13 + labs(x = "Common toad",y="") 
p13 <- p13 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p13



## PLOT 14: Three-spined sticklebank presence and GCN presence/absence 

p14 <- ggplot(box.dat) + ggtitle('(c)')
p14 <- p14 + geom_jitter(aes(x=factor(Gasterosteus_aculeatus), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.3, height=0.3, cex=1)
p14 <- p14 + geom_boxplot(aes(x=factor(Gasterosteus_aculeatus), y=fit), alpha=0.7, outlier.colour = "blue")
p14 <- p14 + scale_colour_gradient(low="grey45", high="orange")
p14 <- p14 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p14 <- p14 + labs(x = "Three-spined stickleback", 
                  y = "") 
p14 <- p14 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p14



## PLOT 15: Grey squirrel presence and GCN presence/absence 

p15 <- ggplot(box.dat) + ggtitle('(d)')
p15 <- p15 + geom_jitter(aes(x=factor(Sciurus_carolinensis), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.3, height=0.3, cex=1)
p15 <- p15 + geom_boxplot(aes(x=factor(Sciurus_carolinensis), y=fit), alpha=0.7, outlier.colour = "blue")
p15 <- p15 + scale_colour_gradient(low="grey45", high="orange")
p15 <- p15 + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p15 <- p15 + labs(x = "Grey squirrel", 
                  y = "Probability of GCN")
p15 <- p15 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p15



## PLOT 16: Outflow and GCN presence/absence 

p16 <- ggplot(box.dat) + ggtitle('(e)')
p16 <- p16 + geom_jitter(aes(x=Outflow, y=Triturus_cristatus, colour=GCNdist$Triturus_cristatus), width=0.3, height=0.3, cex=1)
p16 <- p16 + geom_boxplot(aes(x=Outflow, y=fit), alpha=0.7, outlier.colour = "blue")
p16 <- p16 + scale_colour_gradient(low="grey45", high="orange")
p16 <- p16 + labs(x = "Outflow", 
                  y="") 
p16 <- p16 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p16



## PLOT 17: Ruderals and GCN presence/absence 

f.Ruderals <- ordered(GCNdist$Ruderals, levels = c("None", "Some","Important"))
f.Ruderals

p17 <- ggplot(box.dat) + ggtitle('(f)')
p17 <- p17 + geom_jitter(aes(x=f.Ruderals, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.3, height=0.3, cex=1)
p17 <- p17 + geom_boxplot(aes(x=f.Ruderals, y=fit), alpha=0.7, outlier.colour = "blue")
p17 <- p17 + scale_colour_gradient(low="grey45", high="orange")
p17 <- p17 + labs(x = "Ruderal habitat", y = "")
p17 <- p17 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p17



## PLOT 18: Terrestrial other and GCN presence/absence 

f.Terrestrial_Other <- ordered(GCNdist$Terrestrial_Other, levels = c("None", "Some","Important"))
f.Terrestrial_Other

p18 <- ggplot(box.dat) + ggtitle('(g)')
p18 <- p18 + geom_jitter(aes(x=f.Terrestrial_Other, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.3, height=0.3, cex=1)
p18 <- p18 + geom_boxplot(aes(x=f.Terrestrial_Other, y=fit), alpha=0.7, outlier.colour = "blue")
p18 <- p18 + scale_colour_gradient(low="grey45", high="orange")
p18 <- p18 + labs(x = "Other terrestrial habitat",
                  y = "Probability of GCN") 
p18 <- p18 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p18



## PLOT 19: Species richness and GCN presence/absence

## Create a range of species richness values which increase by 0.01985
## to predict GCN occupancy as species richness increases
range <- seq(from=min(GCNdist$Sp_richness), to=max(GCNdist$Sp_richness), by=0.01985)


## Create new data frame where only species richness changes
## Factors are set to first level in dataset
## Continuous variables are set as the mean values in this study

d1 <- data.frame(sample=rep('F1', length(range)),
                 Lissotriton_vulgaris=rep(0, length(range)),
                 Sp_richness=range,
                 Max_Depth=rep(1.344147, length(range)),
                 Bufo_bufo=rep(0, length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Outflow=rep('Present', length(range)),
                 Ruderals=rep('Important', length(range)),
                 Terrestrial_Other=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNmodel, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.spp <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p19 <- ggplot() + ggtitle('(h)')
p19 <- p19 + geom_jitter(aes(x=GCNdist$Sp_richness, y=GCNdist$Triturus_cristatus, colour=GCNdist$Triturus_cristatus), cex=1, height=0.3)
p19 <- p19 + geom_line(aes(x=dat.spp$Sp_richness, y=dat.spp$fit), size = 1)
p19 <- p19 + geom_ribbon(aes(x=dat.spp$Sp_richness, ymin = dat.spp$cil, ymax = dat.spp$ciu), alpha = 0.25)
p19 <- p19 + scale_colour_gradient(low="grey45", high="orange")
p19 <- p19 + labs(x = "Species richness", y = "")
p19 <- p19 + coord_cartesian(xlim=c(0,10))
p19 <- p19 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p19



## PLOT 20: Pond area and GCN presence/absence

## Create a range of pond area values which increase by 18.6012
## to predict GCN occupancy as pond area increases
range <- seq(from=min(GCNdist$Pond_Area), to=max(GCNdist$Pond_Area), by=18.6012)


## Create new data frame where only pond area changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d2 <- data.frame(sample=rep('F1', length(range)),
                 Lissotriton_vulgaris=rep(0, length(range)),
                 Sp_richness=rep(4, length(range)),
                 Max_Depth=rep(1.344147, length(range)),
                 Bufo_bufo=rep(0, length(range)),
                 Pond_Area=range,
                 Sciurus_carolinensis=rep(0, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Outflow=rep('Present', length(range)),
                 Ruderals=rep('Important', length(range)),
                 Terrestrial_Other=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNmodel, newdata=d2, 
                            se.fit = TRUE, print.matrix=T)) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.PA <- cbind(d2, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between pond area and GCN presence/absence 
## with predicted values from model

p20 <- ggplot() + ggtitle('(i)')
p20 <- p20 + geom_jitter(aes(x=GCNdist$Pond_Area, y=GCNdist$Triturus_cristatus, colour=GCNdist$Triturus_cristatus), height=0.3, cex=1)
p20 <- p20 + geom_line(aes(x=dat.PA$Pond_Area, y=dat.PA$fit), size = 1)
p20 <- p20 + geom_ribbon(aes(x=dat.PA$Pond_Area, ymin = dat.PA$cil, ymax = dat.PA$ciu), alpha = 0.25)
p20 <- p20 + scale_colour_gradient(low="grey45", high="orange")
p20 <- p20 + coord_cartesian(xlim=c(0,10000))
p20 <- p20 + labs(x = "Pond area ("~m^2~")", 
                  y = "")
p20 <- p20 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p20



## PLOT 21: Max. depth and GCN presence/absence

## Create a range of max depth values which increase by 0.01587302
## to predict GCN occupancy as max depth increases
range <- seq(from=min(GCNdist$Max_Depth), to=max(GCNdist$Max_Depth), by=0.01587302)


## Create new data frame where only max depth changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Lissotriton_vulgaris=rep(0, length(range)),
                 Sp_richness=rep(4, length(range)),
                 Max_Depth=range,
                 Bufo_bufo=rep(0, length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Outflow=rep('Present', length(range)),
                 Ruderals=rep('Important', length(range)),
                 Terrestrial_Other=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNmodel, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.PD <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between max depth and GCN presence/absence 
## with predicted values from model

p21 <- ggplot() + ggtitle('(j)')
p21 <- p21 + geom_jitter(aes(x=GCNdist$Max_Depth, y=GCNdist$Triturus_cristatus, colour=GCNdist$Triturus_cristatus), height=0.3, cex=1)
p21 <- p21 + geom_line(aes(x=dat.PD$Max_Depth, y=dat.PD$fit), size = 1)
p21 <- p21 + geom_ribbon(aes(x=dat.PD$Max_Depth, ymin = dat.PD$cil, ymax = dat.PD$ciu), alpha = 0.25)
p21 <- p21 + scale_colour_gradient(low="grey45", high="orange")
p21 <- p21 + coord_cartesian(xlim=c(0,8))
p21 <- p21 + labs(x = "Pond depth (m)",
                  y = "Probability of GCN")
p21 <- p21 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p21




## SHOW ALL PLOTS
## grid.arrange() is a function within the gridExtra package that allows to
## arrange how ggplots are viewed, similar to par(mfrow=c())
## Within this function, I can use arrangeGrob() to specify how many plots I
## want on each row

p12_21 <- grid.arrange(arrangeGrob(p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, 
                                   nrow=5, ncol=2),
                     mylegend, nrow=2,heights=c(10, 1))





#'
#' Now, test HSI score for correlation with GCN presence.
#' 

## Run simple linear regression model
## Use GLM as it is meant to be linear relationship

HSIm <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
              family = binomial,
              data = habitat.dat)


## Compare to null model
HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat.dat)

anova(HSIm, HSInull)


## Model output:
summary(HSIm)
anova(HSIm)
drop1(HSIm, test="Chi")  # for binomial/integer y models, 
                         # statistics for significance of each term

display(HSIm)
se.ranef(HSIm)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSIm)

## Overdispersion test (chi-square)

1-pchisq(581.761, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSIm)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSIm, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSIm)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSIm)


## Plot the residuals against each x vairable
plot(sresid ~ habitat.dat$HSI)


## Check model residuals
chkres(HSIm)

## Plot the fitted data against the observed data
plot(habitat.dat$Triturus_cristatus ~ fitted(HSIm))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSIm)

## 4.99% explained


#'
#' An additional test of model fit is the Hosmer and Lemeshow Goodness 
#' of Fit Test, which is similar to a chi-square test and indicates the 
#' extent to which the model provides better fit than a null model with 
#' no predictors.
#' 
#' If chi-square goodness of fit is not significant, then the model has 
#' adequate fit. Similarly, if the test is significant, the model does not 
#' adequately fit the data.
#' 

hoslem.test(habitat.dat$Triturus_cristatus, fitted(HSIm))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 22: HSI score and GCN presence/absence

## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as HSI score increases
range <- seq(from=min(habitat.dat$HSI), to=max(habitat.dat$HSI), by=0.00142)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d4 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


# Get predictions for this new dataset
fit <- data.frame(predictSE(HSIm, newdata=d4, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d4, fit, ciu, cil) # This is now the data frame


# Plot showing relationship between HSI score and GCN presence/absence 
# with predicted values from model

p22 <- ggplot() + ggtitle('(k)')
p22 <- p22 + geom_jitter(aes(x=habitat.dat$HSI, y=habitat.dat$Triturus_cristatus, colour=habitat.dat$Triturus_cristatus), height=0.3, cex=1)
p22 <- p22 + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p22 <- p22 + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p22 <- p22 + scale_colour_gradient(low="grey45", high="orange")
p22 <- p22 + coord_cartesian(xlim=c(0,1))
p22 <- p22 + labs(x = "Habitat Suitability Index (HSI) score",
                  y = "")
p22 <- p22 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p22


p22mod <- grid.arrange(arrangeGrob(p22, nrow=1, ncol=1),
                       mylegend, heights=c(10, 1))


## Combine HSI score plot with plots from GLMM of biotic and abiotic
## determinants of GCN
## Load 'ggpubr' package
library(ggpubr)
all_plot <- ggarrange(ggarrange(p12,p13,p14,p15,p16,p17,p18,p19,p20, 
                                ncol = 3, nrow = 3),
                      ggarrange(p21,p22, ncol = 2, nrow = 1),
                      heights = c(2, 0.7), nrow = 2) 
p <- grid.arrange(arrangeGrob(all_plot), mylegend, heights=c(15, 1))





#' ---
#' 
#' ## Predictors of Species Richness
#' 
#' Follow a similar approach as predictors of GCN occupancy to identiy
#' predictors of pond species richness.
#' 
#' As the same environmental metadata is being used, there is no need 
#' to examine collinearity between continuous variables again. We know
#' that only pond area should be used to represent pond size and 
#' terrestrial overhang used to represent all terrestrial habitat 
#' variables and instead of shading. Similar, only HSI score should
#' be included as overall habitat score and HSI band are encompassed
#' by this variable.
#' 

## Exploratory plots for species richness

plot(habitat.dat$Sp_richness ~ habitat.dat$Triturus_cristatus)
plot(habitat.dat$Sp_richness ~ habitat.dat$Max_Depth)
plot(habitat.dat$Sp_richness ~ habitat.dat$Pond_Area)
plot(habitat.dat$Sp_richness ~ habitat.dat$Permanance)
plot(habitat.dat$Sp_richness ~ habitat.dat$Quality)
plot(habitat.dat$Sp_richness ~ habitat.dat$Pond_Substrate)
plot(habitat.dat$Sp_richness ~ habitat.dat$Inflow)
plot(habitat.dat$Sp_richness ~ habitat.dat$Outflow)
plot(habitat.dat$Sp_richness ~ habitat.dat$Pollution)
plot(habitat.dat$Sp_richness ~ habitat.dat$Overhang)
plot(habitat.dat$Sp_richness ~ habitat.dat$Amphibians)
plot(habitat.dat$Sp_richness ~ habitat.dat$Waterfowl)
plot(habitat.dat$Sp_richness ~ habitat.dat$Fish)
plot(habitat.dat$Sp_richness ~ habitat.dat$Pond_Density)
plot(habitat.dat$Sp_richness ~ habitat.dat$Woodland)
plot(habitat.dat$Sp_richness ~ habitat.dat$Rough_Grass)
plot(habitat.dat$Sp_richness ~ habitat.dat$Scrub_Hedge)
plot(habitat.dat$Sp_richness ~ habitat.dat$Ruderals)
plot(habitat.dat$Sp_richness ~ habitat.dat$Terrestrial_Other)
plot(habitat.dat$Sp_richness ~ habitat.dat$Overall_Habitat_Score)
plot(habitat.dat$Sp_richness ~ habitat.dat$Macrophytes)
plot(habitat.dat$Sp_richness ~ habitat.dat$HSI)
plot(habitat.dat$Sp_richness ~ habitat.dat$HSI_Band)


## corvif function and pairplot with the Pearson correlation 
## coefficients were imported earlier

## Matrix of Spearman's rank correlations for dataset
## Use Spearman's rank as no assumptions are made about linearity in a
## relationship between two variables
cor(habitat.dat[, 10:19], method = "spearman")


## Booth et al. (1994) and Zuur et al. (2009) suggest that correlations
## between pairs of variables with magnitudes ± 0.5 indicate high
## collinearity

## View pairwise plots of variables
pairs(habitat.dat[, 10:19])


## Below produces a pairwise plot but the strength of colinearity 
## strength is indicated by the size of the text

plot(habitat.dat[, 10:19], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)


#' 
#' There appears to be collinearity between pond circumference, pond
#' length, pond width, and pond area.
#' 
#' Pond area encompasses length and width thus taking the same 
#' measurements and accounting for the same variance in the data as 
#' these variables. Therefore, it is logical to remove pond circumference, 
#' pond length and pond width to counter high collinearity between 
#' explanatory variables.
#' 
#' Shading (% total pond margin shaded) and terrestrial overhang (% pond
#' overhung by trees and shrubs) also appear to be collinear.
#' As terrestrial overhang encompasses the entire pond, wheras shading
#' considers only the pond margin, terrestrial overhang will be retained
#' as an explanatory variable. Shading is a known driver of biodiversity
#' in ponds (Sayer et al. 2012).
#' 
#' Similarly, there is no need to include overall habitat score and HSI
#' band as these are encompassed by HSI score.
#' 
#' 
#' Now check the variance inflation factors (VIFs) of variables to 
#' assess the extent of remaining collinearity.
#' 
#' This is done using a binomial GLM and logit link function containing
#' all explanatory variables to the GCN presence/absence data, then 
#' calculating the VIFs for each variable from the resulting model.
#' 

glmSR <- glm(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes + HSI + HSI_Band
              + Permanance + Quality + Pond_Substrate + Inflow 
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score,
              family = poisson,
              data = habitat.dat)


## Use vif() function in car package to calculate VIFs
vif(glmSR)


#'
#' Many VIF values are below 10 but HSI score, HSI band, pond substrate 
#' are extremely high. This may be because many of the explanatory 
#' variables included in the model are used to calculate HSI score and 
#' thus explain the same variation. 
#' 
#' However, given HSI score is highly correlated with great crested 
#' newt presence (Oldham et al. 2000), it is important to test this
#' variable against species richness. This will be tested separately 
#' after the final habitat model has been achieved.
#' 
#' Pond substrate will be removed due to high VIF.
#' 


## Check possible relationships between response variable and nominal
## variables (factors) in a design plot.

factors <- data.frame(Sp_richness = habitat.dat$Sp_richness,
                      Permanance = habitat.dat$Permanance,
                      Quality = habitat.dat$Quality,
                      Inflow = habitat.dat$Inflow,
                      Outflow = habitat.dat$Outflow,
                      Pollution = habitat.dat$Pollution,
                      Amphibians = habitat.dat$Amphibians,
                      Waterfowl = habitat.dat$Waterfowl,
                      Fish = habitat.dat$Fish,
                      Woodland = habitat.dat$Woodland,
                      Rough_Grass = habitat.dat$Rough_Grass,
                      Scrub_Hedge = habitat.dat$Scrub_Hedge,
                      Ruderals = habitat.dat$Ruderals,
                      Terrestrial_Other = habitat.dat$Terrestrial_Other,
                      Overall_Habitat_Score = habitat.dat$Overall_Habitat_Score)


plot.design(Sp_richness ~ Permanance + Quality + Inflow + Outflow 
            + Pollution + Amphibians + Waterfowl + Fish + Woodland 
            + Rough_Grass + Scrub_Hedge + Ruderals + Terrestrial_Other
            + Overall_Habitat_Score, 
            data = factors, axes = T, xtick = T)


#' 
#' The highest mean values were observed with amphibian and fish 
#' presence.
#' 


## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat.dat, method = "class",
                  minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - overhang
#' - amphibian presence
#' - rough grass habitat
#' - pond density
#' - max depth
#' - pond area
#' - woodland
#' - overall terrestrial habitat quality
#' - ruderals
#' - pollution
#' - fish presence
#' - waterfowl presence
#' - outflow
#' - pond permanence
#' - terrestrial other
#' - water quality
#' - macrophytes
#' - scrub/hedge
#' - inflow
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 5 is optimal i.e. 5 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Overhang, woodland, scrub/hedge, ruderals, rough grass, terrestrial
#' other and overall terrestrial habitat describe terrestrial habitat 
#' in between ponds i.e. the pondscape. 
#' 
#' Presence of amphibians, fish or waterfowl may indicate which groups 
#' provide contain greatest number of species. If habitat is good for 
#' other amphibians (bioindicators), then it is likely to be good for a 
#' number of species.
#' 
#' Terrestrial habitat will greatly influence species richness by 
#' provision of more ecological niches e.g. overhang of ponds and 
#' grassland.
#' 
#' Pond density can influence migration and population viability.
#' 
#' Deeper ponds may inhabit more fish species.
#' 
#' Pond area is complex; larger ponds are likely to be inhabited by
#' more species but may contain more predatory species reducing 
#' number of prey species.
#' 
#' Pollution will affect number of prey species that can utilise a 
#' pond and subsequently number of vertebrate species.
#' 
#' Fish presence is indicative of number of fish species and consequently 
#' macrophyte availability which influences macroinvertebrate
#' species richness, prey source for many aquatic species.
#' 
#' Macrophyte cover is a known driver of biodiversity in ponds and
#' influence species richness of macroinvertebrates.
#' 
#' Outflow may indicate whether pollutants and excess nutrients are flushed 
#' out.
#' 
#' Water quality is an indication of pollution and runoff
#' 
#' Waterfowl presence may be indicative of number of waterfowl species
#' and consequently macrophyte availability which influences macroinvertebrate
#' species richness, prey source for many aquatic species.
#' 
#' Pond inflow may contain pollutants affecting quality, trophic status,
#' and volume.
#' 
#' If ponds are permanant, they are more likely to be used by migratory
#' species.
#' 
#' 
#' Therefore, take all variables forward in model selection.
#' 


#'
#' Several variables occur more than once in the tree which indicates 
#' weak non-linear relationships with the response variable thus a 
#' GAM may be more appropriate to model the data.
#' 
#' However, already know from GCN model that some variables possessed
#' linear relationships and spatial autocorrelation was present in the
#' data and removed by mixed effect modelling.
#' 
#' A GLMM can account for dependencies within sites and is appropriate
#' where nested model selection will be used and the data are structured 
#' hierarchically. Ponds are nested within counties thus a mixed model
#' is necessary to account for spatial dependencies within sites. Each
#' sample represents a different pond, or a different site, and thus
#' sample must be treated as a random effect.
#' 
#' 
#' Apply GLMM to variables most likely to influence species richness
#' based on exploratory analysis
#' 
#' Sample (individual pond) will be random factor as interested in
#' species richness across all ponds as a whole not in specific ponds.
#' 
#' 
#' Identify which variables are important determinants of species richness 
#' and choose a suitable parsimonious approximating model to make predictions. 
#' using an information-theoretic approach (AIC criteria) to evaluate model 
#' fit. Low AIC models are more parsimonious than high AIC models.
#' 
#' 
#' First, group explanatory variables into different functional groups:
#' 
#' - Pond properties (max depth, pond area, pond density, pollution,
#'   macrophytes, outflow, water quality, inflow, permanance)
#'   
#' - Terrestrial habitat (terrestrial overhang, rough grass habitat,
#'   woodland, ruderals, terrestrial other, scrub/hedge, overall
#'   terrestrial habitat)
#'   
#' - Species (amphibian, fish, waterfowl)
#'   


## Create new dataframe containing only explanatory variables to be
## modelled

SRdist <- data.frame(habitat.dat[,c(1,10,14:16,18,20:21,23:34,97)])


## Model with all terms: 

model0 <- glmer(Sp_richness ~ (1|sample) + Overhang + Amphibians 
                + Rough_Grass + Pond_Density + Max_Depth + Pond_Area 
                + Woodland + Ruderals + Pollution + Fish 
                + Terrestrial_Other + Macrophytes + Outflow 
                + Permanance + Inflow + Quality + Scrub_Hedge 
                + Waterfowl + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

summary(model0) # AIC = 2107
drop1(model0, test = "Chi")


## Pond properties:

pondprop <- glmer(Sp_richness ~ (1|sample) + Pond_Density + Max_Depth 
                  + Pond_Area + Pollution + Macrophytes + Outflow 
                  + Permanance + Inflow + Quality,
                  family = poisson,
                  data = SRdist)

summary(pondprop) # AIC = 2086.5, better than first model
drop1(pondprop, test = "Chi")


## Terrestrial habitat:

terr <- glmer(Sp_richness ~ (1|sample) + Overhang + Rough_Grass
              + Woodland + Ruderals + Terrestrial_Other + Scrub_Hedge
              + Overall_Habitat_Score,
              family = poisson,
              data = SRdist)

summary(terr)   # AIC = 2091.4, again better than first model but worse than
                # model containing pond properties
drop1(terr, test = "Chi")


## Species

species <- glmer(Sp_richness ~  (1|sample) 
                 + Amphibians + Waterfowl + Fish,
                 family = poisson,
                 data = SRdist)

summary(species)   # AIC = 2114.4, worse than first model containing all
                   # variables
drop1(species, test = "Chi")


## Try combining pond properties and terrestrial habitat to see if this is
## better model than model with all variables

habitat <- glmer(Sp_richness ~ (1|sample) + Pond_Density + Max_Depth 
                 + Pond_Area + Pollution + Macrophytes + Outflow 
                 + Quality + Inflow + Permanance + Overhang
                 + Rough_Grass + Woodland + Ruderals 
                 + Terrestrial_Other + Scrub_Hedge
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist)

summary(habitat)   # AIC = 2082.8, better than previous models


#' 
#' The AIC comparison of alternate models implies that the explanatory
#' variables best explain the species richness data when they are split
#' into functional groups, and that the model containing pond properties
#' and terrestrial habitat variables best explains species richness.
#' 
#' Apply step-wise down nested model selection on the overall habitat 
#' model.
#' 

model0 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Max_Depth + Pond_Area + Pollution 
                + Macrophytes + Outflow + Quality + Inflow 
                + Permanance + Overhang + Rough_Grass + Woodland 
                + Ruderals + Terrestrial_Other + Scrub_Hedge
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)


#'
#' Use drop1() function to select order that variables should be 
#' dropped from the model. Function computes all the single terms 
#' in the scope argument that can be added to or dropped from the 
#' model, fits those models and computes a table of the changes in 
#' fit.
#' 

drop1(model0)

## Remove woodland as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Max_Depth + Pond_Area + Pollution 
                + Macrophytes + Outflow + Quality + Inflow 
                + Permanance + Overhang + Rough_Grass
                + Ruderals + Terrestrial_Other + Scrub_Hedge
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model1, model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit

drop1(model1)

## Remove permanance as this has next lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Max_Depth + Pond_Area + Pollution 
                + Macrophytes + Outflow + Quality + Inflow 
                + Overhang + Rough_Grass
                + Ruderals + Terrestrial_Other + Scrub_Hedge
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model2, model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of permanance has resulted in better model fit

drop1(model2)

## Remove scrub/hedge as this has next lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Max_Depth + Pond_Area + Pollution 
                + Macrophytes + Outflow + Quality + Inflow 
                + Overhang + Rough_Grass
                + Ruderals + Terrestrial_Other
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model3, model2)


## p-value from LRT not significant and minimal reduction in AIC value 
## (<2). In this instance, it is better to choose model with less 
## parameters to avoid overparameterisation. Removal of scrub/hedge
## has resulted in better model fit

drop1(model3)

## Remove ruderal habitat as this has next lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Max_Depth + Pond_Area + Pollution 
                + Macrophytes + Outflow + Quality + Inflow 
                + Overhang + Rough_Grass
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model4, model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderal habitat has resulted in better model fit

drop1(model4)

## Remove pond depth as this now has lowest AIC value 

model5 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Pond_Area + Pollution 
                + Macrophytes + Outflow + Quality + Inflow 
                + Overhang + Rough_Grass
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model5, model4)

## p-value from LRT not significant and minimal reduction in AIC value 
## (<2). In this instance, it is better to choose model with less 
## parameters to avoid overparameterisation. Removal of pond depth
## has resulted in better model fit

drop1(model5)

## Remove inflow as this now has lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Pond_Area + Pollution 
                + Macrophytes + Outflow + Quality
                + Overhang + Rough_Grass
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model6, model5)

## p-value from LRT not significant and minimal reduction in AIC value 
## (<2). In this instance, it is better to choose model with less 
## parameters to avoid overparameterisation. Removal of inflow
## has improved model fit

drop1(model6)

## Remove pond area next as this has lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Pollution + Macrophytes 
                + Outflow + Quality + Overhang + Rough_Grass
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model7, model6)

## p-value from LRT not significant and minimal reduction in AIC value 
## (<2). In this instance, it is better to choose model with less 
## parameters to avoid overparameterisation. Removal of pond area
## has improved model fit

drop1(model7)

## Remove other terrestrial habitat next as this now has lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Pollution + Macrophytes 
                + Outflow + Quality + Overhang + Rough_Grass
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model8, model7)

## p-value from LRT not significant and minimal reduction in AIC value 
## (<2). In this instance, it is better to choose model with less 
## parameters to avoid overparameterisation. Removal of other terrestrial
## habitat has improved model fit

drop1(model8)

## Remove water quality next as this now has lowest AIC value

model9 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Pollution + Macrophytes 
                + Outflow + Overhang + Rough_Grass
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model9, model8)

## p-value from LRT not significant and minimal increase in AIC value 
## (<2). In this instance, it is better to choose model with less 
## parameters to avoid overparameterisation. Removal of water quality
## has improved model fit

drop1(model9)

## Remove pollution from model as this now has lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)  
                + Pond_Density + Macrophytes + Outflow 
                + Overhang + Rough_Grass + Overall_Habitat_Score,
                family = poisson,
                data = SRdist)

anova(model10, model9)

## p-value from LRT not significant and minimal increase in AIC value 
## (<2). In this instance, it is better to choose model with less 
## parameters to avoid overparameterisation. Removal of pollution
## has improved model fit

drop1(model10)

## Remove macrophyte cover from model next as this has lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)  
                 + Pond_Density + Outflow 
                 + Overhang + Rough_Grass + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist)

anova(model11, model10)

## p-value from LRT significant and increase in AIC value (>2).
## Removal of macrophyte cover has reduced model fit.

drop1(model10)

## Remove pond density next as this has lowest AIC value after 
## macrophyte cover

model12 <- glmer(Sp_richness ~ (1|sample)  
                 + Macrophytes + Outflow 
                 + Overhang + Rough_Grass + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist)

anova(model12, model10)

## p-value from LRT significant and increase in AIC value (>2).
## Removal of pond density has reduced model fit.



## Final model: 

SRmodel <- glmer(Sp_richness ~ (1|sample)  
                 + Pond_Density + Macrophytes + Outflow 
                 + Overhang + Rough_Grass + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist)

summary(SRmodel)
anova(SRmodel)
drop1(SRmodel, test = "Chi")   # for binomial/integer y models, 
                               # statistics for significance of each term

display(SRmodel)
se.ranef(SRmodel)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRmodel)

## Overdispersion test (chi-square)

1-pchisq(509.55, df=494)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRmodel)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRmodel, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRmodel)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRmodel)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist$Overhang)
plot(sresid ~ SRdist$Outflow)
plot(sresid ~ SRdist$Rough_Grass)


## Check model residuals
chkres(SRmodel)

## Plot the fitted data against the observed data
plot(SRdist$Sp_richness ~ fitted(SRmodel))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRmodel)

## 8.94% explained


#'
#' An additional test of model fit is the Hosmer and Lemeshow Goodness 
#' of Fit Test, which is similar to a chi-square test and indicates the 
#' extent to which the model provides better fit than a null model with 
#' no predictors.
#' 
#' If chi-square goodness of fit is not significant, then the model has 
#' adequate fit. Similarly, if the test is significant, the model does not 
#' adequately fit the data.
#' 

hoslem.test(SRdist$Sp_richness, fitted(SRmodel))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: habitat.dat
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRmodel, newdata=SRdist, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist, fit) 


## PLOT 23: Species richness and pond outflow

p23 <- ggplot(box.dat) + ggtitle('(a)')
p23 <- p23 + geom_jitter(aes(x=factor(Outflow), y=Sp_richness), colour="grey45", width=0.3, height=0.3, cex=1)
p23 <- p23 + geom_boxplot(aes(x=factor(Outflow), y=fit), alpha=0.7, outlier.colour = "red")
p23 <- p23 + labs(y="Species richness", x = "Outflow") 
p23 <- p23 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p23



## PLOT 24: Rough grass habitat and species richness

f.Rough_Grass <- ordered(SRdist$Rough_Grass, levels = c("None", "Some","Important"))
f.Rough_Grass

p24 <- ggplot(box.dat) + ggtitle('(b)')
p24 <- p24 + geom_jitter(aes(x=f.Rough_Grass, y=Sp_richness), colour="grey45", width=0.3, height=0.3, cex=1)
p24 <- p24 + geom_boxplot(aes(x=f.Rough_Grass, y=fit), alpha=0.7, outlier.colour = "red")
p24 <- p24 + labs(x = "Rough grass habitat", y="") 
p24 <- p24 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p24



## PLOT 25: Overall terrestrial habitat and species richness

f.Overall_Habitat_Score <- ordered(SRdist$Overall_Habitat_Score, 
                                   levels = c("Poor", "Moderate","Good"))
f.Overall_Habitat_Score

p25 <- ggplot(box.dat) + ggtitle('(c)')
p25 <- p25 + geom_jitter(aes(x=f.Overall_Habitat_Score, y=Sp_richness), colour="grey45", width=0.3, height=0.3, cex=1)
p25 <- p25 + geom_boxplot(aes(x=f.Overall_Habitat_Score, y=fit), alpha=0.7, outlier.colour = "red")
p25 <- p25 + labs(x = "Overall terrestrial habitat quality",
                  y="Species richness") 
p25 <- p25 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p25



## PLOT 26: Overhang and species richness

## Create a range of terrestrial overhang values which increase by 0.1985
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist$Overhang), to=max(SRdist$Overhang), by=0.1985)


## Create new data frame where only terrestrial overhang changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d5 <- data.frame(sample=rep('F1', length(range)),
                 Pond_Density=rep(11.75198, length(range)),
                 Outflow=rep('Present', length(range)),
                 Overhang=range,
                 Macrophytes=rep(25.87302, length(range)),
                 Rough_Grass=rep('Some', length(range)),
                 Overall_Habitat_Score=rep('Moderate', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRmodel, newdata=d5, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d5, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species 
## richness with predicted values from model

p26 <- ggplot() + ggtitle('(d)')
p26 <- p26 + geom_jitter(aes(x=SRdist$Overhang, y=SRdist$Sp_richness), colour="grey45", height=0.7, cex=1)
p26 <- p26 + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p26 <- p26 + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p26 <- p26 + coord_cartesian(xlim=c(0,100))
p26 <- p26 + labs(x = "Terrestrial overhang (%)", y = "")
p26 <- p26 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p26



## PLOT 27: macrophyte cover and species richness

## Create a range of macrophyte cover values which increase by 0.1985
## to predict species richness as macrophyte cover increases
range <- seq(from=min(SRdist$Macrophytes), to=max(SRdist$Macrophytes), by=0.1985)


## Create new data frame where only macrophyte cover changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d6 <- data.frame(sample=rep('F1', length(range)),
                 Pond_Density=rep(11.75198, length(range)),
                 Outflow=rep('Present', length(range)),
                 Overhang=rep(39.36508, length(range)),
                 Macrophytes=range,
                 Rough_Grass=rep('Some', length(range)),
                 Overall_Habitat_Score=rep('Moderate', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRmodel, newdata=d6, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.MC <- cbind(d6, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between macrophyte cover and species richness
## with predicted values from model

p27 <- ggplot() + ggtitle('(e)')
p27 <- p27 + geom_jitter(aes(x=SRdist$Macrophytes, y=SRdist$Sp_richness), colour="grey45", height=0.7, cex=1)
p27 <- p27 + geom_line(aes(x=dat.MC$Macrophytes, y=dat.MC$fit), size = 1)
p27 <- p27 + geom_ribbon(aes(x=dat.MC$Overhang, ymin = dat.MC$cil, ymax = dat.MC$ciu), alpha = 0.25)
p27 <- p27 + coord_cartesian(xlim=c(0,100))
p27 <- p27 + labs(x = "Macrophyte cover (%)", 
                  y = "Species richness")
p27 <- p27 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p27



## PLOT 28: pond density and species richness

## Create a range of pond density values which increase by 0.0879
## to predict species richness as pond density increases
range <- seq(from=min(SRdist$Pond_Density), to=max(SRdist$Pond_Density), by=0.0879)


## Create new data frame where only terrestrial overhang changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d7 <- data.frame(sample=rep('F1', length(range)),
                 Pond_Density=range,
                 Outflow=rep('Present', length(range)),
                 Overhang=rep(39.36508, length(range)),
                 Macrophytes=rep(25.87302, length(range)),
                 Rough_Grass=rep('Some', length(range)),
                 Overall_Habitat_Score=rep('Moderate', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRmodel, newdata=d7, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.PD <- cbind(d7, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p28 <- ggplot() + ggtitle('(f)')
p28 <- p28 + geom_jitter(aes(x=SRdist$Pond_Density, y=SRdist$Sp_richness), colour="grey45", height=0.7, cex=1)
p28 <- p28 + geom_line(aes(x=dat.PD$Pond_Density, y=dat.PD$fit), size = 1)
p28 <- p28 + geom_ribbon(aes(x=dat.PD$Pond_Density, ymin = dat.PD$cil, ymax = dat.PD$ciu), alpha = 0.25)
p28 <- p28 + coord_cartesian(xlim=c(0,50))
p28 <- p28 + labs(x = "Pond density", y = "")
p28 <- p28 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p28




#'
#' Now, test HSI score for correlation with vertebrate species richness.
#' 

## Run simple linear regression model
## Use GLM as it is meant to be linear relationship

srHSIm <- glmer(Sp_richness ~ (1|sample) + HSI,
                family = poisson,
                data = habitat.dat)


## Compare to null model
srHSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                   family = poisson,
                   data = habitat.dat)

anova(srHSIm, srHSInull)

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(srHSIm)
anova(srHSIm)
drop1(srHSIm, test="Chi")  # for binomial/integer y models, 
                           # statistics for significance of each term

display(srHSIm)
se.ranef(srHSIm)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(srHSIm)

## Overdispersion test (chi-square)

1-pchisq(473.844, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(srHSIm)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(srHSIm, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(srHSIm)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(srHSIm)


## Plot the residuals against each x vairable
plot(sresid ~ habitat.dat$HSI)


## Check model residuals
chkres(srHSIm)

## Plot the fitted data against the observed data
plot(habitat.dat$Sp_richness ~ fitted(srHSIm))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(srHSIm)

## only 1.1% explained


#'
#' An additional test of model fit is the Hosmer and Lemeshow Goodness 
#' of Fit Test, which is similar to a chi-square test and indicates the 
#' extent to which the model provides better fit than a null model with 
#' no predictors.
#' 
#' If chi-square goodness of fit is not significant, then the model has 
#' adequate fit. Similarly, if the test is significant, the model does not 
#' adequately fit the data.
#' 

hoslem.test(habitat.dat$Sp_richness, fitted(srHSIm))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 29: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat.dat$HSI), to=max(habitat.dat$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d8 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(srHSIm, newdata=d8, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d8, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p29 <- ggplot() + ggtitle('(g)')
p29 <- p29 + geom_jitter(aes(x=habitat.dat$HSI, y=habitat.dat$Sp_richness), colour="gray45", height=0.7, cex=1)
p29 <- p29 + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p29 <- p29 + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p29 <- p29 + scale_colour_gradient(name="Species\nrichness\n",
                                   low="grey20", high="grey80")
p29 <- p29 + coord_cartesian(xlim=c(0,1))
p29 <- p29 + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p29 <- p29 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title = element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position = "none")
p29



## SHOW ALL PLOTS
## ggarrange is a function within the ggpubr package that allows you to
## arrange how ggplots are viewed, similar to par(mfrow=c())
## Within this function, I can specify how many plots I want on each row
library(ggpubr)
ggarrange(ggarrange(p23,p24,p25,p26,p27,p28, 
                    ncol = 2, nrow = 3, align = "v"),
          p29, heights = c(2, 0.7), nrow = 2)
                      
