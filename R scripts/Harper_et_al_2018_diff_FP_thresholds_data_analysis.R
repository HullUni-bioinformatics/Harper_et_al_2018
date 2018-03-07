#' ---
#' Title: "Threshold analysis for 'Understanding biodiversity at the pondscape using environmental DNA: a focus great crested newts'"
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
#' Species-specific false positive sequence thresholds were applied to the
#' data analysed to address ecological questions performed. Specifically,
#' species associations between great crested newt (GCN) and other 
#' vertebrates, biotic and abiotic determinants of GCN, and determinants
#' of vertebrate species richness. 
#' 
#' The effect of different false positive sequence thresholds on these 
#' analyses will now be investigated. These thresholds were applied in a 
#' separate R script and subsequent datasets saved as csv files.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. Then,
#' load functions for model selection, plotting model residuals and testing 
#' model fit.
#' 

## Clear memory
rm(list=ls())

## Direct R to folder containing packages on workstation
.libPaths(c("C:\\R\\rlib", .libPaths("rlib")))

## set working directory to the location of the script
# install.packages("rstudioapi") # first time for each computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory
getwd()


## Load required packages
p <- c("ggplot2","grid","gridExtra","lme4","glmmADMB","MASS","car","scales",
       "AICcmodavg","gtools","reshape2","plyr","dplyr","arm",
       "sp","RVAideMemoire","ResourceSelection","bbmle","RColorBrewer", 
       "MuMIn","rpart","mgcv","ncf","glmmML","digest","LMERConvenienceFunctions")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies = TRUE)

lapply(p, require, character.only = TRUE)


## Load custom functions
f <- c("CheckResidsFunction.R","OverdispersalFunction.R",
       "CheckConvergenceFunction.R", "HighstatLibV6.R",
       "MyLibrary.R", "glmmML_pres_function.R")

lapply(f, source)


#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()


#'
#' ## Metabarcoding datasets
#' 
#' Next, import the different threshold datasets produced.
#' 

spp_threshold <- read.csv("../Data/MetabarcodingData.csv",
                          header=TRUE)

## Convert species specific thresholds dataset to presence-absence data
spp_threshold[,17:79][spp_threshold[,17:79] > 0] <- 1

## Import other datasets
no_threshold <- read.csv("../Data/data_no_threshold.csv",
                         header=TRUE)

threshold_0.0005 <- read.csv("../Data/data_threshold_0.0005.csv",
                             header=TRUE)

threshold_0.001 <- read.csv("../Data/data_threshold_0.001.csv",
                             header=TRUE)

threshold_0.005 <- read.csv("../Data/data_threshold_0.005.csv",
                             header=TRUE)

threshold_0.01 <- read.csv("../Data/data_threshold_0.01.csv",
                            header=TRUE)

threshold_0.05 <- read.csv("../Data/data_threshold_0.05.csv",
                           header=TRUE)

threshold_0.1 <- read.csv("../Data/data_threshold_0.1.csv",
                          header=TRUE)

threshold_0.3 <- read.csv("../Data/data_threshold_0.3.csv",
                          header=TRUE)


## From each data frame:
## Remove metadata so we only have site by species matrix
## Remove Human as this cannot be ruled out as laboratory contamination
## Remove Clupea harengus, only present in PCR controls 
## Also remove cichlid from species-specific data frame

spp_threshold <- spp_threshold[,-c(2:16,30,44,69,80)]
no_threshold <- no_threshold[,-c(2:14,28,42)]
threshold_0.0005 <- threshold_0.0005[,-c(2:14,28,42)]
threshold_0.001 <- threshold_0.001[,-c(2:14,28,42)]
threshold_0.005 <- threshold_0.005[,-c(2:14,28,42)]
threshold_0.01 <- threshold_0.01[,-c(2:14,28,42)]
threshold_0.05 <- threshold_0.05[,-c(2:14,28,42)]
threshold_0.1 <- threshold_0.1[,-c(2:14,28,42)]
threshold_0.3 <- threshold_0.3[,-c(2:14,28,42)]


#'
#' Examine the effect of other vertebrate species on GCN occupancy.
#' First, I need to create columns containing the number of species
#' in each vertebrate group found in each eDNA sample. Do this for 
#' each data frame.
#' 

## total up number of fish species detected in each sample
spp_threshold$Fish <- rowSums(spp_threshold[,c(3,6:7,12,15:16,19,
                                               24:25,40,45,49,52,
                                               56)])

no_threshold$Fish <- rowSums(no_threshold[,c(3,6:7,12,15:16,19,
                                              24:25,40,45,49,52,
                                              56)])

threshold_0.0005$Fish <- rowSums(threshold_0.0005[,c(3,6:7,12,15:16,19,
                                                      24:25,40,45,49,52,
                                                      56)])

threshold_0.001$Fish <- rowSums(threshold_0.001[,c(3,6:7,12,15:16,19,
                                                    24:25,40,45,49,52,
                                                    56)])


threshold_0.005$Fish <- rowSums(threshold_0.005[,c(3,6:7,12,15:16,19,
                                                   24:25,40,45,49,52,
                                                   56)])

threshold_0.01$Fish <- rowSums(threshold_0.01[,c(3,6:7,12,15:16,19,
                                                 24:25,40,45,49,52,
                                                 56)])

threshold_0.05$Fish <- rowSums(threshold_0.05[,c(3,6:7,12,15:16,19,
                                                 24:25,40,45,49,52,
                                                 56)])

threshold_0.1$Fish <- rowSums(threshold_0.1[,c(3,6:7,12,15:16,19,
                                               24:25,40,45,49,52,
                                               56)])

threshold_0.3$Fish <- rowSums(threshold_0.3[,c(3,6:7,12,15:16,19,
                                               24:25,40,45,49,52,
                                               56)])


## Total up number of amphibian species other than GCN detected in 
## each sample
spp_threshold$Amphibian <- rowSums(spp_threshold[,c(9,29:30,43,50)])

no_threshold$Amphibian <- rowSums(no_threshold[,c(9,29:30,43,50)])

threshold_0.0005$Amphibian <- rowSums(threshold_0.0005[,c(9,29:30,43,50)])

threshold_0.001$Amphibian <- rowSums(threshold_0.001[,c(9,29:30,43,50)])

threshold_0.005$Amphibian <- rowSums(threshold_0.005[,c(9,29:30,43,50)])

threshold_0.01$Amphibian <- rowSums(threshold_0.01[,c(9,29:30,43,50)])

threshold_0.05$Amphibian <- rowSums(threshold_0.05[,c(9,29:30,43,50)])

threshold_0.1$Amphibian <- rowSums(threshold_0.1[,c(9,29:30,43,50)])

threshold_0.3$Amphibian <- rowSums(threshold_0.3[,c(9,29:30,43,50)])


## Total up number of waterfowl species detected in each sample
spp_threshold$Waterfowl <- rowSums(spp_threshold[,c(2,4,21:22,26)])

no_threshold$Waterfowl <- rowSums(no_threshold[,c(2,4,21:22,26)])

threshold_0.0005$Waterfowl <- rowSums(threshold_0.0005[,c(2,4,21:22,26)])

threshold_0.001$Waterfowl <- rowSums(threshold_0.001[,c(2,4,21:22,26)])

threshold_0.005$Waterfowl <- rowSums(threshold_0.005[,c(2,4,21:22,26)])

threshold_0.01$Waterfowl <- rowSums(threshold_0.01[,c(2,4,21:22,26)])

threshold_0.05$Waterfowl <- rowSums(threshold_0.05[,c(2,4,21:22,26)])

threshold_0.1$Waterfowl <- rowSums(threshold_0.1[,c(2,4,21:22,26)])

threshold_0.3$Waterfowl <- rowSums(threshold_0.3[,c(2,4,21:22,26)])


## Total up number of terrestrial bird species detected in each sample
spp_threshold$Terrestrial.Bird <- rowSums(spp_threshold[,c(10,13,17,23,27,
                                                           32,39,44,46,48,
                                                           54,57:58)])
  
no_threshold$Terrestrial.Bird <- rowSums(no_threshold[,c(10,13,17,23,27,
                                                         32,39,44,46,48,
                                                         54,57:58)])

threshold_0.0005$Terrestrial.Bird <- rowSums(threshold_0.0005[,c(10,13,17,23,27,
                                                                 32,39,44,46,48,
                                                                 54,57:58)])

threshold_0.001$Terrestrial.Bird <- rowSums(threshold_0.001[,c(10,13,17,23,27,
                                                               32,39,44,46,48,
                                                               54,57:58)])

threshold_0.005$Terrestrial.Bird <- rowSums(threshold_0.005[,c(10,13,17,23,27,
                                                               32,39,44,46,48,
                                                               54,57:58)])

threshold_0.01$Terrestrial.Bird <- rowSums(threshold_0.01[,c(10,13,17,23,27,
                                                             32,39,44,46,48,
                                                             54,57:58)])

threshold_0.05$Terrestrial.Bird <- rowSums(threshold_0.05[,c(10,13,17,23,27,
                                                             32,39,44,46,48,
                                                             54,57:58)])

threshold_0.1$Terrestrial.Bird <- rowSums(threshold_0.1[,c(10,13,17,23,27,
                                                           32,39,44,46,48,
                                                           54,57:58)])

threshold_0.3$Terrestrial.Bird <- rowSums(threshold_0.3[,c(10,13,17,23,27,
                                                           32,39,44,46,48,
                                                           54,57:58)])


## Total up number of mammal species detected in each sample
spp_threshold$Mammal <- rowSums(spp_threshold[,c(5,8,11,14,18,20,28,
                                                 31,33:38,41:42,47,
                                                 51,53,55,59,61)])

no_threshold$Mammal <- rowSums(no_threshold[,c(5,8,11,14,18,20,28,
                                               31,33:38,41:42,47,
                                               51,53,55,59,61)])

threshold_0.0005$Mammal <- rowSums(threshold_0.0005[,c(5,8,11,14,18,20,28,
                                                       31,33:38,41:42,47,
                                                       51,53,55,59,61)])

threshold_0.001$Mammal <- rowSums(threshold_0.001[,c(5,8,11,14,18,20,28,
                                                     31,33:38,41:42,47,
                                                     51,53,55,59,61)])

threshold_0.005$Mammal <- rowSums(threshold_0.005[,c(5,8,11,14,18,20,28,
                                                     31,33:38,41:42,47,
                                                     51,53,55,59,61)])

threshold_0.01$Mammal <- rowSums(threshold_0.01[,c(5,8,11,14,18,20,28,
                                                   31,33:38,41:42,47,
                                                   51,53,55,59,61)])

threshold_0.05$Mammal <- rowSums(threshold_0.05[,c(5,8,11,14,18,20,28,
                                                   31,33:38,41:42,47,
                                                   51,53,55,59,61)])

threshold_0.1$Mammal <- rowSums(threshold_0.1[,c(5,8,11,14,18,20,28,
                                                 31,33:38,41:42,47,
                                                 51,53,55,59,61)])

threshold_0.3$Mammal <- rowSums(threshold_0.3[,c(5,8,11,14,18,20,28,
                                                 31,33:38,41:42,47,
                                                 51,53,55,59,61)])


## Plot number of fish species present in ponds with and without GCN
## No threshold:
p1a <- ggplot(no_threshold, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1a <- p1a + geom_bar(stat="count", alpha=.5, position='fill')
p1a <- p1a + scale_y_continuous(labels = percent_format())
p1a <- p1a + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1a <- p1a + labs(title="(a) No threshold")
p1a <- p1a + labs(x="Number of fish species", y="Ponds (N = 532)")
p1a <- p1a + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1a <- p1a + guides(fill=guide_legend(reverse=T))
p1a


## 0.05%
p1b <- ggplot(threshold_0.0005, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1b <- p1b + geom_bar(stat="count", alpha=.5, position='fill')
p1b <- p1b + scale_y_continuous(labels = percent_format())
p1b <- p1b + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1b <- p1b + labs(title="(b) 0.05%")
p1b <- p1b + labs(x="Number of fish species", y="Ponds (N = 532)")
p1b <- p1b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1b <- p1b + guides(fill=guide_legend(reverse=T))
p1b


## 0.1%
p1c <- ggplot(threshold_0.001, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1c <- p1c + geom_bar(stat="count", alpha=.5, position='fill')
p1c <- p1c + scale_y_continuous(labels = percent_format())
p1c <- p1c + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1c <- p1c + labs(title="(c) 0.1%")
p1c <- p1c + labs(x="Number of fish species", y="Ponds (N = 532)")
p1c <- p1c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1c <- p1c + guides(fill=guide_legend(reverse=T))
p1c


## 0.5%
p1d <- ggplot(threshold_0.005, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1d <- p1d + geom_bar(stat="count", alpha=.5, position='fill')
p1d <- p1d + scale_y_continuous(labels = percent_format())
p1d <- p1d + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1d <- p1d + labs(title="(d) 0.5%")
p1d <- p1d + labs(x="Number of fish species", y="Ponds (N = 532)")
p1d <- p1d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1d <- p1d + guides(fill=guide_legend(reverse=T))
p1d


## 1%
p1e <- ggplot(threshold_0.01, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1e <- p1e + geom_bar(stat="count", alpha=.5, position='fill')
p1e <- p1e + scale_y_continuous(labels = percent_format())
p1e <- p1e + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1e <- p1e + labs(title="(e) 1%")
p1e <- p1e + labs(x="Number of fish species", y="Ponds (N = 532)")
p1e <- p1e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1e <- p1e + guides(fill=guide_legend(reverse=T))
p1e


## 5%
p1f <- ggplot(threshold_0.05, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1f <- p1f + geom_bar(stat="count", alpha=.5, position='fill')
p1f <- p1f + scale_y_continuous(labels = percent_format())
p1f <- p1f + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1f <- p1f + labs(title="(f) 5%")
p1f <- p1f + labs(x="Number of fish species", y="Ponds (N = 532)")
p1f <- p1f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1f <- p1f + guides(fill=guide_legend(reverse=T))
p1f


## 10%
p1g <- ggplot(threshold_0.1, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1g <- p1g + geom_bar(stat="count", alpha=.5, position='fill')
p1g <- p1g + scale_y_continuous(labels = percent_format())
p1g <- p1g + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1g <- p1g + labs(title="(g) 10%")
p1g <- p1g + labs(x="Number of fish species", y="Ponds (N = 532)")
p1g <- p1g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1g <- p1g + guides(fill=guide_legend(reverse=T))
p1g


## 30%
p1h <- ggplot(threshold_0.3, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1h <- p1h + geom_bar(stat="count", alpha=.5, position='fill')
p1h <- p1h + scale_y_continuous(labels = percent_format())
p1h <- p1h + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p1h <- p1h + labs(title="(h) 30%")
p1h <- p1h + labs(x="Number of fish species", y="Ponds (N = 532)")
p1h <- p1h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p1h <- p1h + guides(fill=guide_legend(reverse=T))
p1h


## species specific thresholds:
p1 <- ggplot(spp_threshold, aes(x=factor(Fish), fill=factor(Triturus_cristatus))) 
p1 <- p1 + geom_bar(stat="count", alpha=.5, position='fill')
p1 <- p1 + scale_y_continuous(labels = percent_format())
p1 <- p1 + scale_fill_manual(values=c("grey48","orange"),
                             breaks=c("0","1"),
                             labels=c("GCN negative",
                                      "GCN positive")) 
p1 <- p1 + labs(title="(i) Species-specific")
p1 <- p1 + labs(x="Number of fish species", y="Ponds (N = 532)")
p1 <- p1 + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.title = element_blank(),
                 text = element_text(size=18))
p1 <- p1 + guides(fill=guide_legend(reverse=T))
p1


# SHOW ALL PLOTS
# grid.arrange() is a function within the gridEXTRA package that allows to
# arrange how ggplots are viewed, similar to par(mfrow=c())
# Within this function, I can use arrangeGrob() to specify how many plots I
# want on each row
# I want to display my factors on the top row and my continuous predictors
# on the bottom row
grid.arrange(arrangeGrob(p1a,p1b,p1c,p1d,p1e, nrow=1, ncol=5),
             arrangeGrob(p1f,p1g,p1h,p1, nrow=1, ncol=4), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of fish species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Plot number of other amphibian species present in ponds with and without 
## GCN
## No threshold:
p2a <- ggplot(no_threshold, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2a <- p2a + geom_bar(stat="count", alpha=.5, position='fill')
p2a <- p2a + scale_y_continuous(labels = percent_format())
p2a <- p2a + scale_fill_manual(values=c("grey48","orange"),
                             breaks=c("0","1"),
                             labels=c("GCN negative",
                                      "GCN positive")) 
p2a <- p2a + labs(title="(a) No threshold")
p2a <- p2a + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2a <- p2a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=18))
p2a <- p2a + guides(fill=guide_legend(reverse=T))
p2a


## 0.05%
p2b <- ggplot(threshold_0.0005, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2b <- p2b + geom_bar(stat="count", alpha=.5, position='fill')
p2b <- p2b + scale_y_continuous(labels = percent_format())
p2b <- p2b + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p2b <- p2b + labs(title="(b) 0.05%")
p2b <- p2b + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2b <- p2b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p2b <- p2b + guides(fill=guide_legend(reverse=T))
p2b


## 0.1%
p2c <- ggplot(threshold_0.001, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2c <- p2c + geom_bar(stat="count", alpha=.5, position='fill')
p2c <- p2c + scale_y_continuous(labels = percent_format())
p2c <- p2c + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p2c <- p2c + labs(title="(c) 0.1%")
p2c <- p2c + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2c <- p2c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p2c <- p2c + guides(fill=guide_legend(reverse=T))
p2c


## 0.5%
p2d <- ggplot(threshold_0.005, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2d <- p2d + geom_bar(stat="count", alpha=.5, position='fill')
p2d <- p2d + scale_y_continuous(labels = percent_format())
p2d <- p2d + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p2d <- p2d + labs(title="(d) 0.5%")
p2d <- p2d + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2d <- p2d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p2d <- p2d + guides(fill=guide_legend(reverse=T))
p2d


## 1%
p2e <- ggplot(threshold_0.01, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2e <- p2e + geom_bar(stat="count", alpha=.5, position='fill')
p2e <- p2e + scale_y_continuous(labels = percent_format())
p2e <- p2e + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p2e <- p2e + labs(title="(e) 1%")
p2e <- p2e + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2e <- p2e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p2e <- p2e + guides(fill=guide_legend(reverse=T))
p2e


## 5%
p2f <- ggplot(threshold_0.05, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2f <- p2f + geom_bar(stat="count", alpha=.5, position='fill')
p2f <- p2f + scale_y_continuous(labels = percent_format())
p2f <- p2f + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p2f <- p2f + labs(title="(f) 5%")
p2f <- p2f + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2f <- p2f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p2f <- p2f + guides(fill=guide_legend(reverse=T))
p2f


## 10%
p2g <- ggplot(threshold_0.1, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2g <- p2g + geom_bar(stat="count", alpha=.5, position='fill')
p2g <- p2g + scale_y_continuous(labels = percent_format())
p2g <- p2g + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p2g <- p2g + labs(title="(g) 10%")
p2g <- p2g + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2g <- p2g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p2g <- p2g + guides(fill=guide_legend(reverse=T))
p2g


## 30%
p2h <- ggplot(threshold_0.3, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus)))
p2h <- p2h + geom_bar(stat="count", alpha=.5, position='fill')
p2h <- p2h + scale_y_continuous(labels = percent_format())
p2h <- p2h + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p2h <- p2h + labs(title="(h) 30%")
p2h <- p2h + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2h <- p2h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p2h <- p2h + guides(fill=guide_legend(reverse=T))
p2h


## Species-specific
p2 <- ggplot(spp_threshold, aes(x=factor(Amphibian), fill=factor(Triturus_cristatus))) 
p2 <- p2 + geom_bar(stat="count", alpha=.5, position='fill')
p2 <- p2 + scale_y_continuous(labels = percent_format())
p2 <- p2 + scale_fill_manual(values=c("grey48","orange"),
                             breaks=c("0","1"),
                             labels=c("GCN negative",
                                      "GCN positive")) 
p2 <- p2 + labs(title="(i) Species-specific")
p2 <- p2 + labs(x="Number of other amphibian species", y="Ponds (N = 532)")
p2 <- p2 + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.title = element_blank(),
                 text = element_text(size=18))
p2 <- p2 + guides(fill=guide_legend(reverse=T))
p2


## SHOW ALL PLOTS
grid.arrange(arrangeGrob(p2a,p2b,p2c,p2d,p2e, nrow=1, ncol=5),
             arrangeGrob(p2f,p2g,p2h,p2, nrow=1, ncol=4), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of other amphibian species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Plot number of waterfowl species in relation to GCN occupancy
## No threshold
p3a <- ggplot(no_threshold, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3a <- p3a + geom_bar(stat="count", alpha=.5, position='fill')
p3a <- p3a + scale_y_continuous(labels = percent_format())
p3a <- p3a + scale_fill_manual(values=c("grey48","orange"),
                             breaks=c("0","1"),
                             labels=c("GCN negative",
                                      "GCN positive")) 
p3a <- p3a + labs(title="(a) No threshold")
p3a <- p3a + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3a <- p3a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=18))
p3a <- p3a + guides(fill=guide_legend(reverse=T))
p3a


## 0.05%
p3b <- ggplot(threshold_0.0005, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3b <- p3b + geom_bar(stat="count", alpha=.5, position='fill')
p3b <- p3b + scale_y_continuous(labels = percent_format())
p3b <- p3b + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3b <- p3b + labs(title="(b) 0.05%")
p3b <- p3b + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3b <- p3b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p3b <- p3b + guides(fill=guide_legend(reverse=T))
p3b


## 0.1%
p3c <- ggplot(threshold_0.001, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3c <- p3c + geom_bar(stat="count", alpha=.5, position='fill')
p3c <- p3c + scale_y_continuous(labels = percent_format())
p3c <- p3c + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3c <- p3c + labs(title="(c) 0.1%")
p3c <- p3c + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3c <- p3c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p3c <- p3c + guides(fill=guide_legend(reverse=T))
p3c


## 0.5%
p3d <- ggplot(threshold_0.005, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3d <- p3d + geom_bar(stat="count", alpha=.5, position='fill')
p3d <- p3d + scale_y_continuous(labels = percent_format())
p3d <- p3d + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3d <- p3d + labs(title="(d) 0.5%")
p3d <- p3d + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3d <- p3d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p3d <- p3d + guides(fill=guide_legend(reverse=T))
p3d


## 1%
p3e <- ggplot(threshold_0.01, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3e <- p3e + geom_bar(stat="count", alpha=.5, position='fill')
p3e <- p3e + scale_y_continuous(labels = percent_format())
p3e <- p3e + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3e <- p3e + labs(title="(e) 1%")
p3e <- p3e + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3e <- p3e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p3e <- p3e + guides(fill=guide_legend(reverse=T))
p3e


## 5%
p3f <- ggplot(threshold_0.05, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3f <- p3f + geom_bar(stat="count", alpha=.5, position='fill')
p3f <- p3f + scale_y_continuous(labels = percent_format())
p3f <- p3f + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3f <- p3f + labs(title="(f) 5%")
p3f <- p3f + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3f <- p3f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p3f <- p3f + guides(fill=guide_legend(reverse=T))
p3f


## 10%
p3g <- ggplot(threshold_0.1, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3g <- p3g + geom_bar(stat="count", alpha=.5, position='fill')
p3g <- p3g + scale_y_continuous(labels = percent_format())
p3g <- p3g + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3g <- p3g + labs(title="(g) 10%")
p3g <- p3g + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3g <- p3g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p3g <- p3g + guides(fill=guide_legend(reverse=T))
p3g


## 30%
p3h <- ggplot(threshold_0.3, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3h <- p3h + geom_bar(stat="count", alpha=.5, position='fill')
p3h <- p3h + scale_y_continuous(labels = percent_format())
p3h <- p3h + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3h <- p3h + labs(title="(h) 30%")
p3h <- p3h + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3h <- p3h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p3h <- p3h + guides(fill=guide_legend(reverse=T))
p3h
p3h


## Species-specific
p3 <- ggplot(spp_threshold, aes(x=factor(Waterfowl), fill=factor(Triturus_cristatus)))
p3 <- p3 + geom_bar(stat="count", alpha=.5, position='fill')
p3 <- p3 + scale_y_continuous(labels = percent_format())
p3 <- p3 + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p3 <- p3 + labs(title="(i) Species-specific")
p3 <- p3 + labs(x="Number of waterfowl species", y="Ponds (N = 532)")
p3 <- p3 + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.title = element_blank(),
                   text = element_text(size=18))
p3 <- p3 + guides(fill=guide_legend(reverse=T))
p3


## SHOW ALL PLOTS
grid.arrange(arrangeGrob(p3a,p3b,p3c,p3d,p3e, nrow=1, ncol=5),
             arrangeGrob(p3f,p3g,p3h,p3, nrow=1, ncol=4), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of waterfowl species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Plot number of terrestrial bird species in relation to GCN
## No threshold
p4a <- ggplot(no_threshold, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4a <- p4a + geom_bar(stat="count", alpha=.5, position='fill')
p4a <- p4a + scale_y_continuous(labels = percent_format())
p4a <- p4a + scale_fill_manual(values=c("grey48","orange"),
                             breaks=c("0","1"),
                             labels=c("GCN negative",
                                      "GCN positive")) 
p4a <- p4a + labs(title="(a) No threshold")
p4a <- p4a + labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4a <- p4a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=18))
p4a <- p4a + guides(fill=guide_legend(reverse=T))
p4a


## 0.05%
p4b <- ggplot(threshold_0.0005, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4b <- p4b + geom_bar(stat="count", alpha=.5, position='fill')
p4b <- p4b + scale_y_continuous(labels = percent_format())
p4b <- p4b + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p4b <- p4b + labs(title="(b) 0.05%")
p4b <- p4b + labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4b <- p4b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p4b <- p4b + guides(fill=guide_legend(reverse=T))
p4b


## 0.1%
p4c <- ggplot(threshold_0.001, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4c <- p4c + geom_bar(stat="count", alpha=.5, position='fill')
p4c <- p4c + scale_y_continuous(labels = percent_format())
p4c <- p4c + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p4c <- p4c + labs(title="(c) 0.1%")
p4c <- p4c + labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4c <- p4c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p4c <- p4c + guides(fill=guide_legend(reverse=T))
p4c


## 0.5%
p4d <- ggplot(threshold_0.005, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4d <- p4d + geom_bar(stat="count", alpha=.5, position='fill')
p4d <- p4d + scale_y_continuous(labels = percent_format())
p4d <- p4d + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p4d <- p4d + labs(title="(d) 0.5%")
p4d <- p4d + labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4d <- p4d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p4d <- p4d + guides(fill=guide_legend(reverse=T))
p4d


## 1%
p4e <- ggplot(threshold_0.01, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4e <- p4e+ geom_bar(stat="count", alpha=.5, position='fill')
p4e <- p4e+ scale_y_continuous(labels = percent_format())
p4e <- p4e+ scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive")) 
p4e <- p4e+ labs(title="(e) 1%")
p4e <- p4e+ labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4e <- p4e+ theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p4e <- p4e+ guides(fill=guide_legend(reverse=T))
p4e


## 5%
p4f <- ggplot(threshold_0.05, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4f <- p4f+ geom_bar(stat="count", alpha=.5, position='fill')
p4f <- p4f+ scale_y_continuous(labels = percent_format())
p4f <- p4f+ scale_fill_manual(values=c("grey48","orange"),
                              breaks=c("0","1"),
                              labels=c("GCN negative",
                                       "GCN positive")) 
p4f <- p4f+ labs(title="(f) 5%")
p4f <- p4f+ labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4f <- p4f+ theme(panel.background = element_rect(fill = 'grey98'),
                  plot.title = element_text(face="bold", hjust=0),
                  legend.position = "none",
                  text = element_text(size=18))
p4f <- p4f+ guides(fill=guide_legend(reverse=T))
p4f


## 10%
p4g <- ggplot(threshold_0.1, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4g <- p4g+ geom_bar(stat="count", alpha=.5, position='fill')
p4g <- p4g+ scale_y_continuous(labels = percent_format())
p4g <- p4g+ scale_fill_manual(values=c("grey48","orange"),
                              breaks=c("0","1"),
                              labels=c("GCN negative",
                                       "GCN positive")) 
p4g <- p4g+ labs(title="(g) 10%")
p4g <- p4g+ labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4g <- p4g+ theme(panel.background = element_rect(fill = 'grey98'),
                  plot.title = element_text(face="bold", hjust=0),
                  legend.position = "none",
                  text = element_text(size=18))
p4g <- p4g+ guides(fill=guide_legend(reverse=T))
p4g


## 30%
p4h <- ggplot(threshold_0.3, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4h <- p4h + geom_bar(stat="count", alpha=.5, position='fill')
p4h <- p4h + scale_y_continuous(labels = percent_format())
p4h <- p4h + scale_fill_manual(values=c("grey48","orange"),
                              breaks=c("0","1"),
                              labels=c("GCN negative",
                                       "GCN positive")) 
p4h <- p4h + labs(title="(h) 30%")
p4h <- p4h + labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4h <- p4h + theme(panel.background = element_rect(fill = 'grey98'),
                  plot.title = element_text(face="bold", hjust=0),
                  legend.position = "none",
                  text = element_text(size=18))
p4h <- p4h + guides(fill=guide_legend(reverse=T))
p4h


## Species specific
p4 <- ggplot(spp_threshold, aes(x=factor(Terrestrial.Bird), fill=factor(Triturus_cristatus)))
p4 <- p4 + geom_bar(stat="count", alpha=.5, position='fill')
p4 <- p4 + scale_y_continuous(labels = percent_format())
p4 <- p4 + scale_fill_manual(values=c("grey48","orange"),
                              breaks=c("0","1"),
                              labels=c("GCN negative",
                                       "GCN positive")) 
p4 <- p4 + labs(title="(i) Species-specific")
p4 <- p4 + labs(x="Number of terrestrial bird species", y="Ponds (N = 532)")
p4 <- p4 + theme(panel.background = element_rect(fill = 'grey98'),
                  plot.title = element_text(face="bold", hjust=0),
                  legend.title = element_blank(),
                  text = element_text(size=18))
p4 <- p4 + guides(fill=guide_legend(reverse=T))
p4


## SHOW ALL PLOTS

grid.arrange(arrangeGrob(p4a,p4b,p4c,p4d,p4e, nrow=1, ncol=5),
             arrangeGrob(p4f,p4g,p4h,p4, nrow=1, ncol=4), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of woodland bird species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Plot number of mammal species in relation to GCN occupancy
## No threshold
p5a <- ggplot(no_threshold, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5a <- p5a + geom_bar(stat="count", alpha=.5, position='fill')
p5a <- p5a + scale_y_continuous(labels = percent_format())
p5a <- p5a + scale_fill_manual(values=c("grey48","orange"),
                             breaks=c("0","1"),
                             labels=c("GCN negative",
                                      "GCN positive"))
p5a <- p5a + labs(title="(a) No threshold")
p5a <- p5a + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5a <- p5a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 legend.position = "none",
                 text = element_text(size=18))
p5a <- p5a + guides(fill=guide_legend(reverse=T))
p5a


## 0.05%
p5b <- ggplot(threshold_0.0005, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5b <- p5b + geom_bar(stat="count", alpha=.5, position='fill')
p5b <- p5b + scale_y_continuous(labels = percent_format())
p5b <- p5b + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5b <- p5b + labs(title="(b) 0.05%")
p5b <- p5b + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5b <- p5b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p5b <- p5b + guides(fill=guide_legend(reverse=T))
p5b


## 0.1%
p5c <- ggplot(threshold_0.001, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5c <- p5c + geom_bar(stat="count", alpha=.5, position='fill')
p5c <- p5c + scale_y_continuous(labels = percent_format())
p5c <- p5c + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5c <- p5c + labs(title="(c) 0.1%")
p5c <- p5c + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5c <- p5c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p5c <- p5c + guides(fill=guide_legend(reverse=T))
p5c


## 0.5%
p5d <- ggplot(threshold_0.005, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5d <- p5d + geom_bar(stat="count", alpha=.5, position='fill')
p5d <- p5d + scale_y_continuous(labels = percent_format())
p5d <- p5d + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5d <- p5d + labs(title="(d) 0.5%")
p5d <- p5d + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5d <- p5d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p5d <- p5d + guides(fill=guide_legend(reverse=T))
p5d


## 1%
p5e <- ggplot(threshold_0.01, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5e <- p5e + geom_bar(stat="count", alpha=.5, position='fill')
p5e <- p5e + scale_y_continuous(labels = percent_format())
p5e <- p5e + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5e <- p5e + labs(title="(e) 1%")
p5e <- p5e + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5e <- p5e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p5e <- p5e + guides(fill=guide_legend(reverse=T))
p5e


## 5%
p5f <- ggplot(threshold_0.05, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5f <- p5f + geom_bar(stat="count", alpha=.5, position='fill')
p5f <- p5f + scale_y_continuous(labels = percent_format())
p5f <- p5f + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5f <- p5f + labs(title="(f) 5%")
p5f <- p5f + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5f <- p5f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p5f <- p5f + guides(fill=guide_legend(reverse=T))
p5f


## 10%
p5g <- ggplot(threshold_0.1, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5g <- p5g + geom_bar(stat="count", alpha=.5, position='fill')
p5g <- p5g + scale_y_continuous(labels = percent_format())
p5g <- p5g + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5g <- p5g + labs(title="(g) 10%")
p5g <- p5g + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5g <- p5g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p5g <- p5g + guides(fill=guide_legend(reverse=T))
p5g


## 30%
p5h <- ggplot(threshold_0.3, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5h <- p5h + geom_bar(stat="count", alpha=.5, position='fill')
p5h <- p5h + scale_y_continuous(labels = percent_format())
p5h <- p5h + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5h <- p5h + labs(title="(h) 30%")
p5h <- p5h + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5h <- p5h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.position = "none",
                   text = element_text(size=18))
p5h <- p5h + guides(fill=guide_legend(reverse=T))
p5h


## Species-specific
p5 <- ggplot(spp_threshold, aes(x=factor(Mammal), fill=factor(Triturus_cristatus)))
p5 <- p5 + geom_bar(stat="count", alpha=.5, position='fill')
p5 <- p5 + scale_y_continuous(labels = percent_format())
p5 <- p5 + scale_fill_manual(values=c("grey48","orange"),
                               breaks=c("0","1"),
                               labels=c("GCN negative",
                                        "GCN positive"))
p5 <- p5 + labs(title="(i) Species-specific")
p5 <- p5 + labs(x="Number of mammal species", y="Ponds (N = 532)")
p5 <- p5 + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   legend.title = element_blank(),
                   text = element_text(size=18))
p5 <- p5 + guides(fill=guide_legend(reverse=T))
p5


## SHOW ALL PLOTS
grid.arrange(arrangeGrob(p5a,p5b,p5c,p5d,p5e, nrow=1, ncol=5),
             arrangeGrob(p5f,p5g,p5h,p5, nrow=1, ncol=4), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of mammal species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Statistically test relationships between number of species and 
## GCN occupancy


####################
# Species-specific #
####################

spp_specific <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                      + Waterfowl + Terrestrial.Bird + Mammal,
                      family = binomial,
                      data = spp_threshold)

summary(spp_specific)  
drop1(spp_specific, test = "Chi")

# amphibian***, waterfowl*** and fish* significant


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(spp_specific)

## Overdispersion test (chi-square)

1-pchisq(588.343, df=525)   # final model overdispersed

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(spp_specific)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(spp_specific)


## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(spp_threshold$Triturus_cristatus, fitted(spp_specific))

## Model does not fit well as p-value is significant 
## i.e. significant difference between the model and the observed 
## data



################
# No threshold #
################

no_thresh <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                   + Waterfowl + Terrestrial.Bird + Mammal,
                   family = binomial,
                   data = no_threshold)

summary(no_thresh)
drop1(no_thresh, test = "Chi")

# amphibian***, waterfowl*** and fish* significant


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(no_thresh)

## Overdispersion test (chi-square)

1-pchisq(636.552, df=525)   # final model overdispersed

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(no_thresh)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(no_thresh)


## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(no_threshold$Triturus_cristatus, fitted(no_thresh))

## Model does not fit well as p-value is significant 
## i.e. significant difference between the model and the observed 
## data



#########
# 0.05% #
#########

thresh0005 <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                   + Waterfowl + Terrestrial.Bird + Mammal,
                   family = binomial,
                   data = threshold_0.0005)

summary(thresh0005)
drop1(thresh0005, test = "Chi") 

# amphibian**, waterfowl** and fish* still significant but less so


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(thresh0005)

## Overdispersion test (chi-square)

1-pchisq(600.96, df=525)   # final model overdispersed

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(thresh0005)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(thresh0005)


## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(threshold_0.0005$Triturus_cristatus, fitted(thresh0005))

## Model does not fit well as p-value is significant 
## i.e. significant difference between the model and the observed 
## data



########
# 0.1% #
########

thresh001 <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                    + Waterfowl + Terrestrial.Bird + Mammal,
                    family = binomial,
                    data = threshold_0.001)

summary(thresh001)
drop1(thresh001, test = "Chi") 

# amphibian**, waterfowl** significant, fish not anymore


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(thresh001)

## Overdispersion test (chi-square)

1-pchisq(591.369, df=525)   # final model overdispersed

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(thresh001)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(thresh001)


## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(threshold_0.001$Triturus_cristatus, fitted(thresh001))

## Model does not fit well as p-value is significant 
## i.e. significant difference between the model and the observed 
## data



########
# 0.5% #
########

thresh005 <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                   + Waterfowl + Terrestrial.Bird + Mammal,
                   family = binomial,
                   data = threshold_0.005)

summary(thresh005)
drop1(thresh005, test = "Chi")

## Fish*, amphibian*, waterfowl*** and woodland bird* significant


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(thresh005)

## Overdispersion test (chi-square)

1-pchisq(565.088, df=525)   # final model not overdispersed

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(thresh005)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(thresh005)


## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(threshold_0.005$Triturus_cristatus, fitted(thresh005))

## Model does fit well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



######
# 1% #
######

thresh01 <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                   + Waterfowl + Terrestrial.Bird + Mammal,
                   family = binomial,
                   data = threshold_0.01)

summary(thresh01)
drop1(thresh01, test = "Chi")

## Fish**, amphibian**, waterfowl*** and woodland bird** significant


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(thresh01)

## Overdispersion test (chi-square)

1-pchisq(527.058, df=525)   # final model not overdispersed

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(thresh01)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(thresh01)


## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(threshold_0.01$Triturus_cristatus, fitted(thresh01))

## Model does not fit well as p-value is significant 
## i.e. significant difference between the model and the observed 
## data



######
# 5% #
######

thresh05 <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                  + Waterfowl + Terrestrial.Bird + Mammal,
                  family = binomial,
                  data = threshold_0.05)

## GLMM could not fit



#######
# 10% #
#######

thresh1 <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian 
                  + Waterfowl + Terrestrial.Bird + Mammal,
                  family = binomial,
                  data = threshold_0.1)

summary(thresh1)
drop1(thresh1, test = "Chi")

## only amphibian*** and mammal** significant


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(thresh1)

## Overdispersion test (chi-square)

1-pchisq(0.807, df=525)   # final model not overdispersed

## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(thresh1)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(thresh1)


## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(threshold_0.1$Triturus_cristatus, fitted(thresh1))

## Model does fit well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



#######
# 30% #
#######

thresh3 <- glmer(Triturus_cristatus ~ (1|sample) + Fish + Amphibian
                 + Waterfowl + Terrestrial.Bird + Mammal,
                 family = binomial,
                 data = threshold_0.3)

## GLMM couldn't fit


#'
#' Significance of explanatory variables was broadly consistent between
#' species-specific thresholds, no threshold, 0.05% and 0.1%. Above this,
#' different variables are found to be significant and some that were 
#' significant at lower thresholds are no longer significant.
#' 

## Plot model output
## Obtain predicted values for the full data set: vert_binary
## se.fit = TRUE will obtain the standard error for each of these 
## predictions

fit_spp <- data.frame(predictSE(spp_specific, newdata=spp_threshold,
                                se.fit = TRUE, print.matrix=T))

fit_NT <- data.frame(predictSE(no_thresh, newdata=no_threshold,
                               se.fit = TRUE, print.matrix=T))

fit_0.0005 <- data.frame(predictSE(thresh0005, newdata=threshold_0.0005,
                                   se.fit = TRUE, print.matrix=T))

fit_0.001 <- data.frame(predictSE(thresh001, newdata=threshold_0.001,
                                  se.fit = TRUE, print.matrix=T))

fit_0.005 <- data.frame(predictSE(thresh005, newdata=threshold_0.005,
                                  se.fit = TRUE, print.matrix=T))

fit_0.01 <- data.frame(predictSE(thresh01, newdata=threshold_0.01,
                                 se.fit = TRUE, print.matrix=T))

fit_0.1 <- data.frame(predictSE(thresh1, newdata=threshold_0.1, 
                                se.fit = TRUE, print.matrix=T))


## Create new data sets with fitted values and original data
spp.dat <- cbind(spp_threshold, fit_spp) 
NT.dat <- cbind(no_threshold, fit_NT) 
T0005.dat <- cbind(threshold_0.0005, fit_0.0005) 
T001.dat <- cbind(threshold_0.001, fit_0.001) 
T005.dat <- cbind(threshold_0.005, fit_0.005)
T01.dat <- cbind(threshold_0.01, fit_0.01) 
T1.dat <- cbind(threshold_0.1, fit_0.1)


## Plots showing relationship between great crested newt presence and 
## species number of each vertebrate group using predicted values from
## models based on different threshold datasets

## Fish
p1ma <- ggplot(NT.dat) 
p1ma <- p1ma + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p1ma <- p1ma + scale_colour_gradient(low="gray65", high="orange")
p1ma <- p1ma + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p1ma <- p1ma + labs(title="(a) No threshold")
p1ma <- p1ma + labs(y="Probability of great crested newt", x = "Number of fish species") 
p1ma <- p1ma + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=18))
p1ma


p1mb <- ggplot(T0005.dat) 
p1mb <- p1mb + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p1mb <- p1mb + scale_colour_gradient(low="gray65", high="orange")
p1mb <- p1mb + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p1mb <- p1mb + labs(title="(b) 0.05%")
p1mb <- p1mb + labs(y="Probability of great crested newt", x = "Number of fish species") 
p1mb <- p1mb + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p1mb


p1mc <- ggplot(T001.dat) 
p1mc <- p1mc + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p1mc <- p1mc + scale_colour_gradient(low="gray65", high="orange")
p1mc <- p1mc + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p1mc <- p1mc + labs(title="(c) 0.1%")
p1mc <- p1mc + labs(y="Probability of great crested newt", x = "Number of fish species") 
p1mc <- p1mc + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p1mc


p1md <- ggplot(T005.dat) 
p1md <- p1md + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p1md <- p1md + scale_colour_gradient(low="gray65", high="orange")
p1md <- p1md + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p1md <- p1md + labs(title="(d) 0.5%")
p1md <- p1md + labs(y="Probability of great crested newt", x = "Number of fish species") 
p1md <- p1md + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p1md


p1me <- ggplot(T01.dat) 
p1me <- p1me + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p1me <- p1me + scale_colour_gradient(low="gray65", high="orange")
p1me <- p1me + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p1me <- p1me + labs(title="(e) 1%")
p1me <- p1me + labs(y="Probability of great crested newt", x = "Number of fish species") 
p1me <- p1me + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p1me


p1mg <- ggplot(T1.dat) 
p1mg <- p1mg + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p1mg <- p1mg + scale_colour_gradient(low="gray65", high="orange")
p1mg <- p1mg + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p1mg <- p1mg + labs(title="(f) 10%")
p1mg <- p1mg + labs(y="Probability of great crested newt", x = "Number of fish species") 
p1mg <- p1mg + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p1mg


p1m <- ggplot(spp.dat) 
p1m <- p1m + geom_jitter(aes(x=factor(Fish), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p1m <- p1m + scale_colour_gradient(low="gray65", high="orange")
p1m <- p1m + geom_boxplot(aes(x=factor(Fish), y=fit), outlier.colour = "blue", alpha=0.7)
p1m <- p1m + labs(title="(g) Species-specific")
p1m <- p1m + labs(y="Probability of great crested newt", x = "Number of fish species") 
p1m <- p1m + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p1m


grid.arrange(arrangeGrob(p1ma,p1mb,p1mc,p1md, nrow=1, ncol=4),
             arrangeGrob(p1me,p1mg,p1m, nrow=1, ncol=3), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of fish species and great crested newt occupancy (GLMM)",
                          gp=gpar(cex=2)))



## Amphibians

p2ma <- ggplot(NT.dat) 
p2ma <- p2ma + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p2ma <- p2ma + scale_colour_gradient(low="gray65", high="orange")
p2ma <- p2ma + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p2ma <- p2ma + labs(title="(a) No threshold")
p2ma <- p2ma + labs(y="Probability of great crested newt", x = "Number of other amphibian species") 
p2ma <- p2ma + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p2ma


p2mb <- ggplot(T0005.dat) 
p2mb <- p2mb + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p2mb <- p2mb + scale_colour_gradient(low="gray65", high="orange")
p2mb <- p2mb + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p2mb <- p2mb + labs(title="(b) 0.05%")
p2mb <- p2mb + labs(y="Probability of great crested newt", x = "Number of other amphibian species") 
p2mb <- p2mb + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p2mb


p2mc <- ggplot(T001.dat) 
p2mc <- p2mc + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p2mc <- p2mc + scale_colour_gradient(low="gray65", high="orange")
p2mc <- p2mc + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p2mc <- p2mc + labs(title="(c) 0.1%")
p2mc <- p2mc + labs(y="Probability of great crested newt", x = "Number of other amphibian species") 
p2mc <- p2mc + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p2mc


p2md <- ggplot(T005.dat) 
p2md <- p2md + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p2md <- p2md + scale_colour_gradient(low="gray65", high="orange")
p2md <- p2md + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p2md <- p2md + labs(title="(d) 0.5%")
p2md <- p2md + labs(y="Probability of great crested newt", x = "Number of other amphibian species") 
p2md <- p2md + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p2md


p2me <- ggplot(T01.dat) 
p2me <- p2me + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p2me <- p2me + scale_colour_gradient(low="gray65", high="orange")
p2me <- p2me + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p2me <- p2me + labs(title="(e) 1%")
p2me <- p2me + labs(y="Probability of great crested newt", x = "Number of other amphibian species") 
p2me <- p2me + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p2me


p2mg <- ggplot(T1.dat) 
p2mg <- p2mg + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p2mg <- p2mg + scale_colour_gradient(low="gray65", high="orange")
p2mg <- p2mg + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p2mg <- p2mg + labs(title="(f) 10%")
p2mg <- p2mg + labs(y="Probability of great crested newt", x = "Number of other amphibian species") 
p2mg <- p2mg + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p2mg


p2m <- ggplot(spp.dat) 
p2m <- p2m + geom_jitter(aes(x=factor(Amphibian), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p2m <- p2m + scale_colour_gradient(low="gray65", high="orange")
p2m <- p2m + geom_boxplot(aes(x=factor(Amphibian), y=fit), outlier.colour = "blue", alpha=0.7)
p2m <- p2m + labs(title="(g) Species-specific")
p2m <- p2m + labs(y="Probability of great crested newt", x = "Number of other amphibian species") 
p2m <- p2m + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=18))
p2m


grid.arrange(arrangeGrob(p2ma,p2mb,p2mc,p2md, nrow=1, ncol=4),
             arrangeGrob(p2me,p2mg,p2m, nrow=1, ncol=3), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of other amphibian species and great crested newt occupancy (GLMM)",
                          gp=gpar(cex=2)))


## Waterfowl
p3ma <- ggplot(NT.dat) 
p3ma <- p3ma + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p3ma <- p3ma + scale_colour_gradient(low="gray65", high="orange")
p3ma <- p3ma + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3ma <- p3ma + labs(title="(a) No threshold")
p3ma <- p3ma + labs(y="Probability of great crested newt", x = "Number of waterfowl species") 
p3ma <- p3ma + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p3ma


p3mb <- ggplot(T0005.dat) 
p3mb <- p3mb + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p3mb <- p3mb + scale_colour_gradient(low="gray65", high="orange")
p3mb <- p3mb + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3mb <- p3mb + labs(title="(b) 0.05%")
p3mb <- p3mb + labs(y="Probability of great crested newt", x = "Number of waterfowl species") 
p3mb <- p3mb + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p3mb


p3mc <- ggplot(T001.dat) 
p3mc <- p3mc + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p3mc <- p3mc + scale_colour_gradient(low="gray65", high="orange")
p3mc <- p3mc + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3mc <- p3mc + labs(title="(c) 0.1%")
p3mc <- p3mc + labs(y="Probability of great crested newt", x = "Number of waterfowl species") 
p3mc <- p3mc + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p3mc


p3md <- ggplot(T005.dat) 
p3md <- p3md + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p3md <- p3md + scale_colour_gradient(low="gray65", high="orange")
p3md <- p3md + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3md <- p3md + labs(title="(d) 0.5%")
p3md <- p3md + labs(y="Probability of great crested newt", x = "Number of waterfowl species") 
p3md <- p3md + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p3md


p3me <- ggplot(T01.dat) 
p3me <- p3me + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p3me <- p3me + scale_colour_gradient(low="gray65", high="orange")
p3me <- p3me + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3me <- p3me + labs(title="(e) 1%")
p3me <- p3me + labs(y="Probability of great crested newt", x = "Number of waterfowl species") 
p3me <- p3me + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p3me


p3mg <- ggplot(T1.dat) 
p3mg <- p3mg + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p3mg <- p3mg + scale_colour_gradient(low="gray65", high="orange")
p3mg <- p3mg + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3mg <- p3mg + labs(title="(f) 10%")
p3mg <- p3mg + labs(y="Probability of great crested newt", x = "Number of waterfowl species") 
p3mg <- p3mg + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p3mg


p3m <- ggplot(spp.dat) 
p3m <- p3m + geom_jitter(aes(x=factor(Waterfowl), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p3m <- p3m + scale_colour_gradient(low="gray65", high="orange")
p3m <- p3m + geom_boxplot(aes(x=factor(Waterfowl), y=fit), outlier.colour = "blue", alpha=0.7)
p3m <- p3m + labs(title="(g) Species-specific")
p3m <- p3m + labs(y="Probability of great crested newt", x = "Number of waterfowl species") 
p3m <- p3m + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=18))
p3m


grid.arrange(arrangeGrob(p3ma,p3mb,p3mc,p3md, nrow=1, ncol=4),
             arrangeGrob(p3me,p3mg,p3m, nrow=1, ncol=3), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of waterfowl species and great crested newt occupancy (GLMM)",
                          gp=gpar(cex=2)))


## Terrestrial birds
p4ma <- ggplot(NT.dat) 
p4ma <- p4ma + geom_jitter(aes(x=factor(Terrestrial.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p4ma <- p4ma + scale_colour_gradient(low="gray65", high="orange")
p4ma <- p4ma + geom_boxplot(aes(x=factor(Terrestrial.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4ma <- p4ma + labs(title="(a) No threshold")
p4ma <- p4ma + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p4ma <- p4ma + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p4ma


p4mb <- ggplot(T0005.dat) 
p4mb <- p4mb + geom_jitter(aes(x=factor(Terrestrial.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p4mb <- p4mb + scale_colour_gradient(low="gray65", high="orange")
p4mb <- p4mb + geom_boxplot(aes(x=factor(Terrestrial.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4mb <- p4mb + labs(title="(b) 0.05%")
p4mb <- p4mb + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p4mb <- p4mb + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p4mb


p4mc <- ggplot(T001.dat) 
p4mc <- p4mc + geom_jitter(aes(x=factor(Terrestrial.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p4mc <- p4mc + scale_colour_gradient(low="gray65", high="orange")
p4mc <- p4mc + geom_boxplot(aes(x=factor(Terrestrial.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4mc <- p4mc + labs(title="(c) 0.1%")
p4mc <- p4mc + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p4mc <- p4mc + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p4mc


p4md <- ggplot(T005.dat) 
p4md <- p4md + geom_jitter(aes(x=factor(Terrestrial.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p4md <- p4md + scale_colour_gradient(low="gray65", high="orange")
p4md <- p4md + geom_boxplot(aes(x=factor(Terrestrial.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4md <- p4md + labs(title="(d) 0.5%")
p4md <- p4md + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p4md <- p4md + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p4md


p4me <- ggplot(T01.dat) 
p4me <- p4me + geom_jitter(aes(x=factor(Terrestrial.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p4me <- p4me + scale_colour_gradient(low="gray65", high="orange")
p4me <- p4me + geom_boxplot(aes(x=factor(Terrestrial.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4me <- p4me + labs(title="(e) 1%")
p4me <- p4me + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p4me <- p4me + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p4me


p4mg <- ggplot(T1.dat) 
p4mg <- p4mg + geom_jitter(aes(x=factor(Terrestrial.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p4mg <- p4mg + scale_colour_gradient(low="gray65", high="orange")
p4mg <- p4mg + geom_boxplot(aes(x=factor(Terrestrial.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4mg <- p4mg + labs(title="(f) 10%")
p4mg <- p4mg + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p4mg <- p4mg + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p4mg


p4m <- ggplot(spp.dat) 
p4m <- p4m + geom_jitter(aes(x=factor(Terrestrial.Bird), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p4m <- p4m + scale_colour_gradient(low="gray65", high="orange")
p4m <- p4m + geom_boxplot(aes(x=factor(Terrestrial.Bird), y=fit), outlier.colour = "blue", alpha=0.7)
p4m <- p4m + labs(title="(g) Species-specific")
p4m <- p4m + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p4m <- p4m + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=18))
p4m


grid.arrange(arrangeGrob(p4ma,p4mb,p4mc,p4md, nrow=1, ncol=4),
             arrangeGrob(p4me,p4mg,p4m, nrow=1, ncol=3), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of woodland bird species and great crested newt occupancy (GLMM)",
                          gp=gpar(cex=2)))


## Mammals
p5ma <- ggplot(NT.dat) 
p5ma <- p5ma + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p5ma <- p5ma + scale_colour_gradient(low="gray65", high="orange")
p5ma <- p5ma + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5ma <- p5ma + labs(title="(a) No threshold")
p5ma <- p5ma + labs(y="Probability of great crested newt", x = "Number of mammal species") 
p5ma <- p5ma + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p5ma


p5mb <- ggplot(T0005.dat) 
p5mb <- p5mb + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p5mb <- p5mb + scale_colour_gradient(low="gray65", high="orange")
p5mb <- p5mb + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5mb <- p5mb + labs(title="(b) 0.05%")
p5mb <- p5mb + labs(y="Probability of great crested newt", x = "Number of mammal species") 
p5mb <- p5mb + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p5mb


p5mc <- ggplot(T001.dat) 
p5mc <- p5mc + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p5mc <- p5mc + scale_colour_gradient(low="gray65", high="orange")
p5mc <- p5mc + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5mc <- p5mc + labs(title="(c) 0.1%")
p5mc <- p5mc + labs(y="Probability of great crested newt", x = "Number of mammal species") 
p5mc <- p5mc + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p5mc


p5md <- ggplot(T005.dat) 
p5md <- p5md + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p5md <- p5md + scale_colour_gradient(low="gray65", high="orange")
p5md <- p5md + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5md <- p5md + labs(title="(d) 0.5%")
p5md <- p5md + labs(y="Probability of great crested newt", x = "Number of mammal species") 
p5md <- p5md + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p5md


p5me <- ggplot(T01.dat) 
p5me <- p5me + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p5me <- p5me + scale_colour_gradient(low="gray65", high="orange")
p5me <- p5me + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5me <- p5me + labs(title="(e) 1%")
p5me <- p5me + labs(y="Probability of great crested newt", x = "Number of mammal species") 
p5me <- p5me + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p5me


p5mg <- ggplot(T1.dat) 
p5mg <- p5mg + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p5mg <- p5mg + scale_colour_gradient(low="gray65", high="orange")
p5mg <- p5mg + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5mg <- p5mg + labs(title="(f) 10%")
p5mg <- p5mg + labs(y="Probability of great crested newt", x = "Number of mammal species") 
p5mg <- p5mg + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     legend.position="none",
                     text = element_text(size=18))
p5mg


p5m <- ggplot(spp.dat) 
p5m <- p5m + geom_jitter(aes(x=factor(Mammal), y=Triturus_cristatus, colour=Triturus_cristatus), height=0.7, width=0.5, cex=1, alpha=0.7)
p5m <- p5m + scale_colour_gradient(low="gray65", high="orange")
p5m <- p5m + geom_boxplot(aes(x=factor(Mammal), y=fit), outlier.colour = "blue", alpha=0.7)
p5m <- p5m + labs(title="(g) Species-specific")
p5m <- p5m + labs(y="Probability of great crested newt", x = "Number of woodland bird species") 
p5m <- p5m + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title =element_text(face="bold", hjust=0),
                   legend.position="none",
                   text = element_text(size=18))
p5m


grid.arrange(arrangeGrob(p5ma,p5mb,p5mc,p5md, nrow=1, ncol=4),
             arrangeGrob(p5me,p5mg,p5m, nrow=1, ncol=3), 
             nrow=2,
             top=textGrob("Effect of sequence thresholds on number of mammal species and great crested newt occupancy (GLMM)",
                          gp=gpar(cex=2)))



#'
#' Now, show individual species in each threshold dataset. 
#' 


## total up number times an assignment occurred across GCN positive
## ponds
GCNpos_SS <- subset(spp_threshold, spp_threshold$Triturus_cristatus == 1)
GCNpos_NT <- subset(no_threshold, no_threshold$Triturus_cristatus == 1)
GCNpos_0.0005 <- subset(threshold_0.0005, threshold_0.0005$Triturus_cristatus == 1)
GCNpos_0.001 <- subset(threshold_0.001, threshold_0.001$Triturus_cristatus == 1)
GCNpos_0.005 <- subset(threshold_0.005, threshold_0.005$Triturus_cristatus == 1)
GCNpos_0.01 <- subset(threshold_0.01, threshold_0.01$Triturus_cristatus == 1)
GCNpos_0.05 <- subset(threshold_0.05, threshold_0.05$Triturus_cristatus == 1)
GCNpos_0.1 <- subset(threshold_0.1, threshold_0.1$Triturus_cristatus == 1)
GCNpos_0.3 <- subset(threshold_0.1, threshold_0.1$Triturus_cristatus == 1)

GCNneg_SS <- subset(spp_threshold, spp_threshold$Triturus_cristatus == 0)
GCNneg_NT <- subset(no_threshold, no_threshold$Triturus_cristatus == 0)
GCNneg_0.0005 <- subset(threshold_0.0005, threshold_0.0005$Triturus_cristatus == 0)
GCNneg_0.001 <- subset(threshold_0.001, threshold_0.001$Triturus_cristatus == 0)
GCNneg_0.005 <- subset(threshold_0.005, threshold_0.005$Triturus_cristatus == 0)
GCNneg_0.01 <- subset(threshold_0.01, threshold_0.01$Triturus_cristatus == 0)
GCNneg_0.05 <- subset(threshold_0.05, threshold_0.05$Triturus_cristatus == 0)
GCNneg_0.1 <- subset(threshold_0.1, threshold_0.1$Triturus_cristatus == 0)
GCNneg_0.3 <- subset(threshold_0.1, threshold_0.1$Triturus_cristatus == 0)


## Create dataframes containing only species in GCN positive ponds
GCNpos_SSspp <- data.frame(GCNpos_SS[,c(2:59,61)])
GCNpos_NTspp <- data.frame(GCNpos_NT[,c(2:59,61)])
GCNpos_0.0005spp <- data.frame(GCNpos_0.0005[,c(2:59,61)])
GCNpos_0.001spp <- data.frame(GCNpos_0.001[,c(2:59,61)])
GCNpos_0.005spp <- data.frame(GCNpos_0.005[,c(2:59,61)])
GCNpos_0.01spp <- data.frame(GCNpos_0.01[,c(2:59,61)])
GCNpos_0.05spp <- data.frame(GCNpos_0.05[,c(2:59,61)])
GCNpos_0.1spp <- data.frame(GCNpos_0.1[,c(2:59,61)])
GCNpos_0.3spp <- data.frame(GCNpos_0.3[,c(2:59,61)])


## Calculate number of ponds each species found in
GCNpos_SSspp <- rbind(GCNpos_SSspp, colSums(GCNpos_SSspp))
GCNpos_NTspp <- rbind(GCNpos_NTspp, colSums(GCNpos_NTspp))
GCNpos_0.0005spp <- rbind(GCNpos_0.0005spp, colSums(GCNpos_0.0005spp))
GCNpos_0.001spp <- rbind(GCNpos_0.001spp, colSums(GCNpos_0.001spp))
GCNpos_0.005spp <- rbind(GCNpos_0.005spp, colSums(GCNpos_0.005spp))
GCNpos_0.01spp <- rbind(GCNpos_0.01spp, colSums(GCNpos_0.01spp))
GCNpos_0.05spp <- rbind(GCNpos_0.05spp, colSums(GCNpos_0.05spp))
GCNpos_0.1spp <- rbind(GCNpos_0.1spp, colSums(GCNpos_0.1spp))
GCNpos_0.3spp <- rbind(GCNpos_0.3spp, colSums(GCNpos_0.3spp))


## Create dataframes containing only species in GCN negative ponds
GCNneg_SSspp <- data.frame(GCNneg_SS[,c(2:59,61)])
GCNneg_NTspp <- data.frame(GCNneg_NT[,c(2:59,61)])
GCNneg_0.0005spp <- data.frame(GCNneg_0.0005[,c(2:59,61)])
GCNneg_0.001spp <- data.frame(GCNneg_0.001[,c(2:59,61)])
GCNneg_0.005spp <- data.frame(GCNneg_0.005[,c(2:59,61)])
GCNneg_0.01spp <- data.frame(GCNneg_0.01[,c(2:59,61)])
GCNneg_0.05spp <- data.frame(GCNneg_0.05[,c(2:59,61)])
GCNneg_0.1spp <- data.frame(GCNneg_0.1[,c(2:59,61)])
GCNneg_0.3spp <- data.frame(GCNneg_0.3[,c(2:59,61)])


## Calculate number of ponds each species found in
GCNneg_SSspp <- rbind(GCNneg_SSspp, colSums(GCNneg_SSspp))
GCNneg_NTspp <- rbind(GCNneg_NTspp, colSums(GCNneg_NTspp))
GCNneg_0.0005spp <- rbind(GCNneg_0.0005spp, colSums(GCNneg_0.0005spp))
GCNneg_0.001spp <- rbind(GCNneg_0.001spp, colSums(GCNneg_0.001spp))
GCNneg_0.005spp <- rbind(GCNneg_0.005spp, colSums(GCNneg_0.005spp))
GCNneg_0.01spp <- rbind(GCNneg_0.01spp, colSums(GCNneg_0.01spp))
GCNneg_0.05spp <- rbind(GCNneg_0.05spp, colSums(GCNneg_0.05spp))
GCNneg_0.1spp <- rbind(GCNneg_0.1spp, colSums(GCNneg_0.1spp))
GCNneg_0.3spp <- rbind(GCNneg_0.3spp, colSums(GCNneg_0.3spp))


## Bind rows containing pond totals from each dataframe
totals_SS <- data.frame(rbind(GCNpos_SSspp[150,], GCNneg_SSspp[384,]))
rownames(totals_SS) <- c("P", "N")

totals_NT <- data.frame(rbind(GCNpos_NTspp[183,], GCNneg_NTspp[351,]))
rownames(totals_NT) <- c("P", "N")

totals_0.0005 <- data.frame(rbind(GCNpos_0.0005spp[149,], GCNneg_0.0005spp[385,]))
rownames(totals_0.0005) <- c("P", "N")

totals_0.001 <- data.frame(rbind(GCNpos_0.001spp[144,], GCNneg_0.001spp[390,]))
rownames(totals_0.001) <- c("P", "N")

totals_0.005 <- data.frame(rbind(GCNpos_0.005spp[133,], GCNneg_0.005spp[401,]))
rownames(totals_0.005) <- c("P", "N")

totals_0.01 <- data.frame(rbind(GCNpos_0.01spp[122,], GCNneg_0.01spp[412,]))
rownames(totals_0.01) <- c("P", "N")

totals_0.05 <- data.frame(rbind(GCNpos_0.05spp[82,], GCNneg_0.05spp[452,]))
rownames(totals_0.05) <- c("P", "N")

totals_0.1 <- data.frame(rbind(GCNpos_0.1spp[57,], GCNneg_0.1spp[477,]))
rownames(totals_0.1) <- c("P", "N")

totals_0.3 <- data.frame(rbind(GCNpos_0.3spp[57,], GCNneg_0.3spp[477,]))
rownames(totals_0.3) <- c("P", "N")


## Transpose data frameS
totals_SS <- t(totals_SS)
totals_NT <- t(totals_NT)
totals_0.0005 <- t(totals_0.0005)
totals_0.001 <- t(totals_0.001)
totals_0.005 <- t(totals_0.005)
totals_0.01 <- t(totals_0.01)
totals_0.05 <- t(totals_0.05)
totals_0.1 <- t(totals_0.1)
totals_0.3 <- t(totals_0.3)


## Melt dataframes to create new factor column containing P and N
ponds_spp_SS <- melt(totals_SS, id=c("P","N"))
colnames(ponds_spp_SS) <- c("Species", "GCN", "Ponds")

ponds_spp_NT <- melt(totals_NT, id=c("P","N"))
colnames(ponds_spp_NT) <- c("Species", "GCN", "Ponds")

ponds_spp_0.0005 <- melt(totals_0.0005, id=c("P","N"))
colnames(ponds_spp_0.0005) <- c("Species", "GCN", "Ponds")

ponds_spp_0.001 <- melt(totals_0.001, id=c("P","N"))
colnames(ponds_spp_0.001) <- c("Species", "GCN", "Ponds")

ponds_spp_0.005 <- melt(totals_0.005, id=c("P","N"))
colnames(ponds_spp_0.005) <- c("Species", "GCN", "Ponds")

ponds_spp_0.01 <- melt(totals_0.01, id=c("P","N"))
colnames(ponds_spp_0.01) <- c("Species", "GCN", "Ponds")

ponds_spp_0.05 <- melt(totals_0.05, id=c("P","N"))
colnames(ponds_spp_0.05) <- c("Species", "GCN", "Ponds")

ponds_spp_0.1 <- melt(totals_0.1, id=c("P","N"))
colnames(ponds_spp_0.1) <- c("Species", "GCN", "Ponds")

ponds_spp_0.3 <- melt(totals_0.3, id=c("P","N"))
colnames(ponds_spp_0.3) <- c("Species", "GCN", "Ponds")


## To complete the dataframe, need to create another factor column
## indicating what vertebrate groups species belong to

ponds_spp_SS$Group <- factor(ifelse(ponds_spp_SS$Species %in% ponds_spp_SS[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                 ifelse(ponds_spp_SS$Species %in% ponds_spp_SS[c(8,28:29,42,49),1], "Amphibian",
                                        ifelse(ponds_spp_SS$Species %in% ponds_spp_SS[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                               ifelse(ponds_spp_SS$Species %in% ponds_spp_SS[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                      "NA")))))

ponds_spp_NT$Group <- factor(ifelse(ponds_spp_NT$Species %in% ponds_spp_NT[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                    ifelse(ponds_spp_NT$Species %in% ponds_spp_NT[c(8,28:29,42,49),1], "Amphibian",
                                           ifelse(ponds_spp_NT$Species %in% ponds_spp_NT[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                  ifelse(ponds_spp_NT$Species %in% ponds_spp_NT[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                         "NA")))))

ponds_spp_0.0005$Group <- factor(ifelse(ponds_spp_0.0005$Species %in% ponds_spp_0.0005[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                    ifelse(ponds_spp_0.0005$Species %in% ponds_spp_0.0005[c(8,28:29,42,49),1], "Amphibian",
                                           ifelse(ponds_spp_0.0005$Species %in% ponds_spp_0.0005[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                  ifelse(ponds_spp_0.0005$Species %in% ponds_spp_0.0005[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                         "NA")))))

ponds_spp_0.001$Group <- factor(ifelse(ponds_spp_0.001$Species %in% ponds_spp_0.001[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                        ifelse(ponds_spp_0.001$Species %in% ponds_spp_0.001[c(8,28:29,42,49),1], "Amphibian",
                                               ifelse(ponds_spp_0.001$Species %in% ponds_spp_0.001[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                      ifelse(ponds_spp_0.001$Species %in% ponds_spp_0.001[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                             "NA")))))

ponds_spp_0.005$Group <- factor(ifelse(ponds_spp_0.005$Species %in% ponds_spp_0.005[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                       ifelse(ponds_spp_0.005$Species %in% ponds_spp_0.005[c(8,28:29,42,49),1], "Amphibian",
                                              ifelse(ponds_spp_0.005$Species %in% ponds_spp_0.005[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                     ifelse(ponds_spp_0.005$Species %in% ponds_spp_0.005[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                            "NA")))))

ponds_spp_0.01$Group <- factor(ifelse(ponds_spp_0.01$Species %in% ponds_spp_0.01[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                       ifelse(ponds_spp_0.01$Species %in% ponds_spp_0.01[c(8,28:29,42,49),1], "Amphibian",
                                              ifelse(ponds_spp_0.01$Species %in% ponds_spp_0.01[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                     ifelse(ponds_spp_0.01$Species %in% ponds_spp_0.01[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                            "NA")))))

ponds_spp_0.05$Group <- factor(ifelse(ponds_spp_0.05$Species %in% ponds_spp_0.05[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                      ifelse(ponds_spp_0.05$Species %in% ponds_spp_0.05[c(8,28:29,42,49),1], "Amphibian",
                                             ifelse(ponds_spp_0.05$Species %in% ponds_spp_0.05[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                    ifelse(ponds_spp_0.05$Species %in% ponds_spp_0.05[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                           "NA")))))

ponds_spp_0.1$Group <- factor(ifelse(ponds_spp_0.1$Species %in% ponds_spp_0.1[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                      ifelse(ponds_spp_0.1$Species %in% ponds_spp_0.1[c(8,28:29,42,49),1], "Amphibian",
                                             ifelse(ponds_spp_0.1$Species %in% ponds_spp_0.1[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                    ifelse(ponds_spp_0.1$Species %in% ponds_spp_0.1[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                           "NA")))))

ponds_spp_0.3$Group <- factor(ifelse(ponds_spp_0.3$Species %in% ponds_spp_0.3[c(2,5:6,11,14:15,18,23:24,39,44,48,51,55),1], "Fish",
                                     ifelse(ponds_spp_0.3$Species %in% ponds_spp_0.3[c(8,28:29,42,49),1], "Amphibian",
                                            ifelse(ponds_spp_0.3$Species %in% ponds_spp_0.3[c(1,3,9,12,16,20:22,25:26,31,38,43,45,47,53,56:57),1], "Bird",
                                                   ifelse(ponds_spp_0.3$Species %in% ponds_spp_0.3[c(4,7,10,13,17,19,27,30,32:37,40:41,46,50,52,54,58:59),1], "Mammal",
                                                          "NA")))))


## Remove underscores from species names
ponds_spp_SS$Species <- gsub('_', ' ', ponds_spp_SS$Species)
ponds_spp_NT$Species <- gsub('_', ' ', ponds_spp_NT$Species)
ponds_spp_0.0005$Species <- gsub('_', ' ', ponds_spp_0.0005$Species)
ponds_spp_0.001$Species <- gsub('_', ' ', ponds_spp_0.001$Species)
ponds_spp_0.005$Species <- gsub('_', ' ', ponds_spp_0.005$Species)
ponds_spp_0.01$Species <- gsub('_', ' ', ponds_spp_0.01$Species)
ponds_spp_0.05$Species <- gsub('_', ' ', ponds_spp_0.05$Species)
ponds_spp_0.1$Species <- gsub('_', ' ', ponds_spp_0.1$Species)
ponds_spp_0.3$Species <- gsub('_', ' ', ponds_spp_0.3$Species)


## Plot proportion of fish species in GCN and non-GCN ponds across
## different threshold datasets
p6a <- ggplot(subset(ponds_spp_NT, Group == "Fish"))
p6a <- p6a + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6a <- p6a + scale_y_continuous(labels = percent_format())
p6a <- p6a + scale_fill_manual(values=c("orange","grey48"),
                             breaks=c("P","N"),
                             labels=c("GCN positive",
                                      "GCN negative"))
p6a <- p6a + labs(title="(a) No threshold")
p6a <- p6a + labs(x="Species", y="Number of Ponds (N = 532)")
p6a <- p6a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p6a <- p6a + guides(fill=guide_legend(reverse=T))
p6a


p6b <- ggplot(subset(ponds_spp_0.0005, Group == "Fish"))
p6b <- p6b + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6b <- p6b + scale_y_continuous(labels = percent_format())
p6b <- p6b + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6b <- p6b + labs(title="(b) 0.05%")
p6b <- p6b + labs(x="Species", y="Number of Ponds (N = 532)")
p6b <- p6b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6b <- p6b + guides(fill=guide_legend(reverse=T))
p6b


p6c <- ggplot(subset(ponds_spp_0.001, Group == "Fish"))
p6c <- p6c + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6c <- p6c + scale_y_continuous(labels = percent_format())
p6c <- p6c + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6c <- p6c + labs(title="(c) 0.1%")
p6c <- p6c + labs(x="Species", y="Number of Ponds (N = 532)")
p6c <- p6c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6c <- p6c + guides(fill=guide_legend(reverse=T))
p6c


p6d <- ggplot(subset(ponds_spp_0.005, Group == "Fish"))
p6d <- p6d + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6d <- p6d + scale_y_continuous(labels = percent_format())
p6d <- p6d + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6d <- p6d + labs(title="(d) 0.5%")
p6d <- p6d + labs(x="Species", y="Number of Ponds (N = 532)")
p6d <- p6d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6d <- p6d + guides(fill=guide_legend(reverse=T))
p6d


p6e <- ggplot(subset(ponds_spp_0.01, Group == "Fish"))
p6e <- p6e + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6e <- p6e + scale_y_continuous(labels = percent_format())
p6e <- p6e + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6e <- p6e + labs(title="(e) 1%")
p6e <- p6e + labs(x="Species", y="Number of Ponds (N = 532)")
p6e <- p6e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6e <- p6e + guides(fill=guide_legend(reverse=T))
p6e


p6f <- ggplot(subset(ponds_spp_0.05, Group == "Fish"))
p6f <- p6f + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6f <- p6f + scale_y_continuous(labels = percent_format())
p6f <- p6f + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6f <- p6f + labs(title="(f) 5%")
p6f <- p6f + labs(x="Species", y="Number of Ponds (N = 532)")
p6f <- p6f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6f <- p6f + guides(fill=guide_legend(reverse=T))
p6f


p6g <- ggplot(subset(ponds_spp_0.1, Group == "Fish"))
p6g <- p6g + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6g <- p6g + scale_y_continuous(labels = percent_format())
p6g <- p6g + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6g <- p6g + labs(title="(g) 10%")
p6g <- p6g + labs(x="Species", y="Number of Ponds (N = 532)")
p6g <- p6g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6g <- p6g + guides(fill=guide_legend(reverse=T))
p6g


p6h <- ggplot(subset(ponds_spp_0.3, Group == "Fish"))
p6h <- p6h + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6h <- p6h + scale_y_continuous(labels = percent_format())
p6h <- p6h + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6h <- p6h + labs(title="(h) 30%")
p6h <- p6h + labs(x="Species", y="Number of Ponds (N = 532)")
p6h <- p6h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6h <- p6h + guides(fill=guide_legend(reverse=T))
p6h


p6 <- ggplot(subset(ponds_spp_SS, Group == "Fish"))
p6 <- p6 + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p6 <- p6 + scale_y_continuous(labels = percent_format())
p6 <- p6 + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p6 <- p6 + labs(title="(i) Species-specific")
p6 <- p6 + labs(x="Species", y="Number of Ponds (N = 532)")
p6 <- p6 + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p6 <- p6 + guides(fill=guide_legend(reverse=T))
p6


grid.arrange(arrangeGrob(p6a,p6b,p6c, nrow=1, ncol=3),
             arrangeGrob(p6d,p6e,p6f, nrow=1, ncol=3),
             arrangeGrob(p6g,p6h,p6, nrow=1, ncol=3), 
             nrow=3,
             top=textGrob("Effect of sequence thresholds on fish species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Plot proportion of amphibian species in GCN and non-GCN ponds
p7a <- ggplot(subset(ponds_spp_NT, Group == "Amphibian"))
p7a <- p7a + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7a <- p7a + scale_y_continuous(labels = percent_format())
p7a <- p7a + scale_fill_manual(values=c("orange","grey48"),
                             breaks=c("P","N"),
                             labels=c("GCN positive",
                                      "GCN negative"))
p7a <- p7a + labs(title="(a) No threshold")
p7a <- p7a + labs(x="Species", y="Number of Ponds (N = 532)")
p7a <- p7a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p7a <- p7a + guides(fill=guide_legend(reverse=T))
p7a


p7b <- ggplot(subset(ponds_spp_0.0005, Group == "Amphibian"))
p7b <- p7b + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7b <- p7b + scale_y_continuous(labels = percent_format())
p7b <- p7b + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p7b <- p7b + labs(title="(b) 0.05%")
p7b <- p7b + labs(x="Species", y="Number of Ponds (N = 532)")
p7b <- p7b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p7b <- p7b + guides(fill=guide_legend(reverse=T))
p7b


p7c <- ggplot(subset(ponds_spp_0.001, Group == "Amphibian"))
p7c <- p7c + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7c <- p7c + scale_y_continuous(labels = percent_format())
p7c <- p7c + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p7c <- p7c + labs(title="(c) 0.1%")
p7c <- p7c + labs(x="Species", y="Number of Ponds (N = 532)")
p7c <- p7c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p7c <- p7c + guides(fill=guide_legend(reverse=T))
p7c


p7d <- ggplot(subset(ponds_spp_0.005, Group == "Amphibian"))
p7d <- p7d + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7d <- p7d + scale_y_continuous(labels = percent_format())
p7d <- p7d + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p7d <- p7d + labs(title="(d) 0.5%")
p7d <- p7d + labs(x="Species", y="Number of Ponds (N = 532)")
p7d <- p7d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p7d <- p7d + guides(fill=guide_legend(reverse=T))
p7d


p7e <- ggplot(subset(ponds_spp_0.01, Group == "Amphibian"))
p7e <- p7e + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7e <- p7e + scale_y_continuous(labels = percent_format())
p7e <- p7e + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p7e <- p7e + labs(title="(e) 1%")
p7e <- p7e + labs(x="Species", y="Number of Ponds (N = 532)")
p7e <- p7e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p7e <- p7e + guides(fill=guide_legend(reverse=T))
p7e


p7f <- ggplot(subset(ponds_spp_0.05, Group == "Amphibian"))
p7f <- p7f + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7f <- p7f + scale_y_continuous(labels = percent_format())
p7f <- p7f + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p7f <- p7f + labs(title="(f) 5%")
p7f <- p7f + labs(x="Species", y="Number of Ponds (N = 532)")
p7f <- p7f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p7f <- p7f + guides(fill=guide_legend(reverse=T))
p7f


p7g <- ggplot(subset(ponds_spp_0.1, Group == "Amphibian"))
p7g <- p7g + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7g <- p7g + scale_y_continuous(labels = percent_format())
p7g <- p7g + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p7g <- p7g + labs(title="(g) 10%")
p7g <- p7g + labs(x="Species", y="Number of Ponds (N = 532)")
p7g <- p7g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p7g <- p7g + guides(fill=guide_legend(reverse=T))
p7g


p7h <- ggplot(subset(ponds_spp_0.3, Group == "Amphibian"))
p7h <- p7h + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7h <- p7h + scale_y_continuous(labels = percent_format())
p7h <- p7h + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p7h <- p7h + labs(title="(h) 30%")
p7h <- p7h + labs(x="Species", y="Number of Ponds (N = 532)")
p7h <- p7h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p7h <- p7h + guides(fill=guide_legend(reverse=T))
p7h


p7 <- ggplot(subset(ponds_spp_SS, Group == "Amphibian"))
p7 <- p7 + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p7 <- p7 + scale_y_continuous(labels = percent_format())
p7 <- p7 + scale_fill_manual(values=c("orange","grey48"),
                             breaks=c("P","N"),
                             labels=c("GCN positive",
                                      "GCN negative"))
p7 <- p7 + labs(title="(i) Species-specific")
p7 <- p7 + labs(x="Species", y="Number of Ponds (N = 532)")
p7 <- p7 + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p7 <- p7 + guides(fill=guide_legend(reverse=T))
p7


grid.arrange(arrangeGrob(p7a,p7b,p7c, nrow=1, ncol=3),
             arrangeGrob(p7d,p7e,p7f, nrow=1, ncol=3),
             arrangeGrob(p7g,p7h,p7, nrow=1, ncol=3), 
             nrow=3,
             top=textGrob("Effect of sequence thresholds on other amphibian species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Plot bird species in GCN and non-GCN ponds
p8a <- ggplot(subset(ponds_spp_NT, Group == "Bird"))
p8a <- p8a + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8a <- p8a + scale_y_continuous(labels = percent_format())
p8a <- p8a + scale_fill_manual(values=c("orange","grey48"),
                             breaks=c("P","N"),
                             labels=c("GCN positive",
                                      "GCN negative"))
p8a <- p8a + labs(title="(a) No threshold")
p8a <- p8a + labs(x="Species", y="Number of Ponds (N = 532)")
p8a <- p8a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p8a <- p8a + guides(fill=guide_legend(reverse=T))
p8a


p8b <- ggplot(subset(ponds_spp_0.0005, Group == "Bird"))
p8b <- p8b + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8b <- p8b + scale_y_continuous(labels = percent_format())
p8b <- p8b + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p8b <- p8b + labs(title="(b) 0.05%")
p8b <- p8b + labs(x="Species", y="Number of Ponds (N = 532)")
p8b <- p8b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p8b <- p8b + guides(fill=guide_legend(reverse=T))
p8b


p8d <- ggplot(subset(ponds_spp_0.001, Group == "Bird"))
p8c <- p8c + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8c <- p8c + scale_y_continuous(labels = percent_format())
p8c <- p8c + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p8c <- p8c + labs(title="(c) 0.1%")
p8c <- p8c + labs(x="Species", y="Number of Ponds (N = 532)")
p8c <- p8c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p8c <- p8c + guides(fill=guide_legend(reverse=T))
p8c


p8d <- ggplot(subset(ponds_spp_0.005, Group == "Bird"))
p8d <- p8d + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8d <- p8d + scale_y_continuous(labels = percent_format())
p8d <- p8d + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p8d <- p8d + labs(title="(d) 0.5%")
p8d <- p8d + labs(x="Species", y="Number of Ponds (N = 532)")
p8d <- p8d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p8d <- p8d + guides(fill=guide_legend(reverse=T))
p8d


p8e <- ggplot(subset(ponds_spp_0.01, Group == "Bird"))
p8e <- p8e + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8e <- p8e + scale_y_continuous(labels = percent_format())
p8e <- p8e + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p8e <- p8e + labs(title="(e) 1%")
p8e <- p8e + labs(x="Species", y="Number of Ponds (N = 532)")
p8e <- p8e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p8e <- p8e + guides(fill=guide_legend(reverse=T))
p8e


p8f <- ggplot(subset(ponds_spp_0.05, Group == "Bird"))
p8f <- p8f + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8f <- p8f + scale_y_continuous(labels = percent_format())
p8f <- p8f + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p8f <- p8f + labs(title="(f) 5%")
p8f <- p8f + labs(x="Species", y="Number of Ponds (N = 532)")
p8f <- p8f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p8f <- p8f + guides(fill=guide_legend(reverse=T))
p8f


p8g <- ggplot(subset(ponds_spp_0.1, Group == "Bird"))
p8g <- p8g + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8g <- p8g + scale_y_continuous(labels = percent_format())
p8g <- p8g + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p8g <- p8g + labs(title="(g) 10%")
p8g <- p8g + labs(x="Species", y="Number of Ponds (N = 532)")
p8g <- p8g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p8g <- p8g + guides(fill=guide_legend(reverse=T))
p8g


p8h <- ggplot(subset(ponds_spp_0.3, Group == "Bird"))
p8h <- p8h + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8h <- p8h + scale_y_continuous(labels = percent_format())
p8h <- p8h + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p8h <- p8h + labs(title="(h) 30%")
p8h <- p8h + labs(x="Species", y="Number of Ponds (N = 532)")
p8h <- p8h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p8h <- p8h + guides(fill=guide_legend(reverse=T))
p8h


p8 <- ggplot(subset(ponds_spp_SS, Group == "Bird"))
p8 <- p8 + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p8 <- p8 + scale_y_continuous(labels = percent_format())
p8 <- p8 + scale_fill_manual(values=c("orange","grey48"),
                             breaks=c("P","N"),
                             labels=c("GCN positive",
                                      "GCN negative"))
p8 <- p8 + labs(title="(i) Species-specific")
p8 <- p8 + labs(x="Species", y="Number of Ponds (N = 532)")
p8 <- p8 + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p8 <- p8 + guides(fill=guide_legend(reverse=T))
p8


grid.arrange(arrangeGrob(p8a,p8b,p8c, nrow=1, ncol=3),
             arrangeGrob(p8d,p8e,p8f, nrow=1, ncol=3),
             arrangeGrob(p8g,p8h,p8, nrow=1, ncol=3), 
             nrow=3,
             top=textGrob("Effect of sequence thresholds on bird species and great crested newt occupancy",
                          gp=gpar(cex=2)))



## Plot mammal species in GCN and non-GCN ponds
p9a <- ggplot(subset(ponds_spp_NT, Group == "Mammal"))
p9a <- p9a + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9a <- p9a + scale_y_continuous(labels = percent_format())
p9a <- p9a + scale_fill_manual(values=c("orange","grey48"),
                             breaks=c("P","N"),
                             labels=c("GCN positive",
                                      "GCN negative"))
p9a <- p9a + labs(title="(a) No threshold")
p9a <- p9a + labs(x="Species", y="Number of Ponds (N = 532)")
p9a <- p9a + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p9a <- p9a + guides(fill=guide_legend(reverse=T))
p9a


p9b <- ggplot(subset(ponds_spp_0.0005, Group == "Mammal"))
p9b <- p9b + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9b <- p9b + scale_y_continuous(labels = percent_format())
p9b <- p9b + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p9b <- p9b + labs(title="(b) 0.05%")
p9b <- p9b + labs(x="Species", y="Number of Ponds (N = 532)")
p9b <- p9b + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p9b <- p9b + guides(fill=guide_legend(reverse=T))
p9b


p9c <- ggplot(subset(ponds_spp_0.001, Group == "Mammal"))
p9c <- p9c + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9c <- p9c + scale_y_continuous(labels = percent_format())
p9c <- p9c + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p9c <- p9c + labs(title="(c) 0.1%")
p9c <- p9c + labs(x="Species", y="Number of Ponds (N = 532)")
p9c <- p9c + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p9c <- p9c + guides(fill=guide_legend(reverse=T))
p9c


p9d <- ggplot(subset(ponds_spp_0.005, Group == "Mammal"))
p9d <- p9d + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9d <- p9d + scale_y_continuous(labels = percent_format())
p9d <- p9d + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p9d <- p9d + labs(title="(d) 0.5%")
p9d <- p9d + labs(x="Species", y="Number of Ponds (N = 532)")
p9d <- p9d + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p9d <- p9d + guides(fill=guide_legend(reverse=T))
p9d


p9e <- ggplot(subset(ponds_spp_0.01, Group == "Mammal"))
p9e <- p9e + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9e <- p9e + scale_y_continuous(labels = percent_format())
p9e <- p9e + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p9e <- p9e + labs(title="(e) 1%")
p9e <- p9e + labs(x="Species", y="Number of Ponds (N = 532)")
p9e <- p9e + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p9e <- p9e + guides(fill=guide_legend(reverse=T))
p9e


p9f <- ggplot(subset(ponds_spp_0.05, Group == "Mammal"))
p9f <- p9f + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9f <- p9f + scale_y_continuous(labels = percent_format())
p9f <- p9f + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p9f <- p9f + labs(title="(f) 5%")
p9f <- p9f + labs(x="Species", y="Number of Ponds (N = 532)")
p9f <- p9f + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p9f <- p9f + guides(fill=guide_legend(reverse=T))
p9f


p9g <- ggplot(subset(ponds_spp_0.1, Group == "Mammal"))
p9g <- p9g + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9g <- p9g + scale_y_continuous(labels = percent_format())
p9g <- p9g + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p9g <- p9g + labs(title="(g) 10%")
p9g <- p9g + labs(x="Species", y="Number of Ponds (N = 532)")
p9g <- p9g + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p9g <- p9g + guides(fill=guide_legend(reverse=T))
p9g


p9h <- ggplot(subset(ponds_spp_0.3, Group == "Mammal"))
p9h <- p9h + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9h <- p9h + scale_y_continuous(labels = percent_format())
p9h <- p9h + scale_fill_manual(values=c("orange","grey48"),
                               breaks=c("P","N"),
                               labels=c("GCN positive",
                                        "GCN negative"))
p9h <- p9h + labs(title="(h) 30%")
p9h <- p9h + labs(x="Species", y="Number of Ponds (N = 532)")
p9h <- p9h + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                   legend.title = element_blank(),
                   text = element_text(size=24))
p9h <- p9h + guides(fill=guide_legend(reverse=T))
p9h


p9 <- ggplot(subset(ponds_spp_SS, Group == "Mammal"))
p9 <- p9 + geom_bar(stat="identity", alpha=.5, aes(x=Species, y=Ponds, fill=GCN), position="fill")
p9 <- p9 + scale_y_continuous(labels = percent_format())
p9 <- p9 + scale_fill_manual(values=c("orange","grey48"),
                             breaks=c("P","N"),
                             labels=c("GCN positive",
                                      "GCN negative"))
p9 <- p9 + labs(title="(i) Species-specific")
p9 <- p9 + labs(x="Species", y="Number of Ponds (N = 532)")
p9 <- p9 + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p9 <- p9 + guides(fill=guide_legend(reverse=T))
p9


grid.arrange(arrangeGrob(p9a,p9b,p9c, nrow=1, ncol=3),
             arrangeGrob(p9d,p9e,p9f, nrow=1, ncol=3),
             arrangeGrob(p9g,p9h,p9, nrow=1, ncol=3), 
             nrow=3,
             top=textGrob("Effect of sequence thresholds on mammal species and great crested newt occupancy",
                          gp=gpar(cex=2)))



#'
#' ## Significance of species associations at different thresholds
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

## Create new data frames for each threshold and remove data used for
## previous analysis
ponds_SS <- spp_threshold[,-c(62:66)]
ponds_NT <- no_threshold[,-c(62:66)]
ponds_0.0005 <- threshold_0.0005[,-c(62:66)]
ponds_0.001 <- threshold_0.001[,-c(62:66)]
ponds_0.005 <- threshold_0.005[,-c(62:66)]
ponds_0.01 <- threshold_0.01[,-c(62:66)]
ponds_0.05 <- threshold_0.05[,-c(62:66)]
ponds_0.1 <- threshold_0.1[,-c(62:66)]
ponds_0.3 <- threshold_0.3[,-c(62:66)]

## Transpose data frames so that species are first column
ponds_SS <- t(ponds_SS)
colnames(ponds_SS) <- ponds_SS[1,]
ponds_SS <- ponds_SS[-1,]
ponds_SS <- data.frame(ponds_SS)
ponds_SS <- tibble::rownames_to_column(ponds_SS, "species")

ponds_NT <- t(ponds_NT)
colnames(ponds_NT) <- ponds_NT[1,]
ponds_NT <- ponds_NT[-1,]
ponds_NT <- data.frame(ponds_NT)
ponds_NT <- tibble::rownames_to_column(ponds_NT, "species")

ponds_0.0005 <- t(ponds_0.0005)
colnames(ponds_0.0005) <- ponds_0.0005[1,]
ponds_0.0005 <- ponds_0.0005[-1,]
ponds_0.0005 <- data.frame(ponds_0.0005)
ponds_0.0005 <- tibble::rownames_to_column(ponds_0.0005, "species")

ponds_0.001 <- t(ponds_0.001)
colnames(ponds_0.001) <- ponds_0.001[1,]
ponds_0.001 <- ponds_0.001[-1,]
ponds_0.001 <- data.frame(ponds_0.001)
ponds_0.001 <- tibble::rownames_to_column(ponds_0.001, "species")

ponds_0.005 <- t(ponds_0.005)
colnames(ponds_0.005) <- ponds_0.005[1,]
ponds_0.005 <- ponds_0.005[-1,]
ponds_0.005 <- data.frame(ponds_0.005)
ponds_0.005 <- tibble::rownames_to_column(ponds_0.005, "species")

ponds_0.01 <- t(ponds_0.01)
colnames(ponds_0.01) <- ponds_0.01[1,]
ponds_0.01 <- ponds_0.01[-1,]
ponds_0.01 <- data.frame(ponds_0.01)
ponds_0.01 <- tibble::rownames_to_column(ponds_0.01, "species")

ponds_0.05 <- t(ponds_0.05)
colnames(ponds_0.05) <- ponds_0.05[1,]
ponds_0.05 <- ponds_0.05[-1,]
ponds_0.05 <- data.frame(ponds_0.05)
ponds_0.05 <- tibble::rownames_to_column(ponds_0.05, "species")

ponds_0.1 <- t(ponds_0.1)
colnames(ponds_0.1) <- ponds_0.1[1,]
ponds_0.1 <- ponds_0.1[-1,]
ponds_0.1 <- data.frame(ponds_0.1)
ponds_0.1 <- tibble::rownames_to_column(ponds_0.1, "species")

ponds_0.3 <- t(ponds_0.3)
colnames(ponds_0.3) <- ponds_0.3[1,]
ponds_0.3 <- ponds_0.3[-1,]
ponds_0.3 <- data.frame(ponds_0.3)
ponds_0.3 <- tibble::rownames_to_column(ponds_0.3, "species")


## correct the species names
ponds_SS$species <- gsub('_', ' ', ponds_SS$species)
ponds_NT$species <- gsub('_', ' ', ponds_NT$species)
ponds_0.0005$species <- gsub('_', ' ', ponds_0.0005$species)
ponds_0.001$species <- gsub('_', ' ', ponds_0.001$species)
ponds_0.005$species <- gsub('_', ' ', ponds_0.005$species)
ponds_0.01$species <- gsub('_', ' ', ponds_0.01$species)
ponds_0.05$species <- gsub('_', ' ', ponds_0.05$species)
ponds_0.1$species <- gsub('_', ' ', ponds_0.1$species)
ponds_0.3$species <- gsub('_', ' ', ponds_0.3$species)


## Export new data frames as csv files
write.csv(ponds_SS, "../Data/SS_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_NT, "../Data/NT_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_0.0005, "../Data/0.0005_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_0.001, "../Data/0.001_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_0.005, "../Data/0.005_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_0.01, "../Data/0.01_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_0.05, "../Data/0.05_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_0.1, "../Data/0.1_cooccur_data.csv", 
          row.names=FALSE)

write.csv(ponds_0.3, "../Data/0.3_cooccur_data.csv", 
          row.names=FALSE)


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
#' Do this for the metabarcoding data at all false positive thresholds.
#' 

## Load datasets containing presence-absence data of all species
## detected in each pond by metabarcoding at different thresholds

ponds_SS <- read.csv("../Data/SS_cooccur_data.csv", 
                     row.names=1)


cooccur.ponds_SS <- cooccur(mat = ponds_SS,     # the dataframe
                            type = "spp_site",  # data organized with species as rows and sites as column
                            thresh = TRUE,      # summarize the most important species associations
                            spp_names = TRUE)   # use species names in data

class(cooccur.ponds_SS)
summary(cooccur.ponds_SS)
prob.table(cooccur.ponds_SS)
print(cooccur.ponds_SS)



ponds_NT <- read.csv("../Data/NT_cooccur_data.csv", 
                     row.names=1)


cooccur.ponds_NT <- cooccur(mat = ponds_NT,     # the dataframe
                            type = "spp_site",  # data organized with species as rows and sites as column
                            thresh = TRUE,      # summarize the most important species associations
                            spp_names = TRUE)   # use species names in data

class(cooccur.ponds_NT)
summary(cooccur.ponds_NT)
prob.table(cooccur.ponds_NT)
print(cooccur.ponds_NT)



ponds_0.0005 <- read.csv("../Data/0.0005_cooccur_data.csv", 
                     row.names=1)


cooccur.ponds_0.0005 <- cooccur(mat = ponds_0.0005,     # the dataframe
                            type = "spp_site",  # data organized with species as rows and sites as column
                            thresh = TRUE,      # summarize the most important species associations
                            spp_names = TRUE)   # use species names in data

class(cooccur.ponds_0.0005)
summary(cooccur.ponds_0.0005)
prob.table(cooccur.ponds_0.0005)
print(cooccur.ponds_0.0005)



ponds_0.001 <- read.csv("../Data/0.001_cooccur_data.csv", 
                         row.names=1)


cooccur.ponds_0.001 <- cooccur(mat = ponds_0.001,     # the dataframe
                               type = "spp_site",  # data organized with species as rows and sites as column
                               thresh = TRUE,      # summarize the most important species associations
                               spp_names = TRUE)   # use species names in data

class(cooccur.ponds_0.001)
summary(cooccur.ponds_0.001)
prob.table(cooccur.ponds_0.001)
print(cooccur.ponds_0.001)



ponds_0.005 <- read.csv("../Data/0.005_cooccur_data.csv", 
                        row.names=1)


cooccur.ponds_0.005 <- cooccur(mat = ponds_0.005,     # the dataframe
                               type = "spp_site",  # data organized with species as rows and sites as column
                               thresh = TRUE,      # summarize the most important species associations
                               spp_names = TRUE)   # use species names in data

class(cooccur.ponds_0.005)
summary(cooccur.ponds_0.005)
prob.table(cooccur.ponds_0.005)
print(cooccur.ponds_0.005)



ponds_0.01 <- read.csv("../Data/0.01_cooccur_data.csv", 
                        row.names=1)


cooccur.ponds_0.01 <- cooccur(mat = ponds_0.01,     # the dataframe
                               type = "spp_site",  # data organized with species as rows and sites as column
                               thresh = TRUE,      # summarize the most important species associations
                               spp_names = TRUE)   # use species names in data

class(cooccur.ponds_0.01)
summary(cooccur.ponds_0.01)
prob.table(cooccur.ponds_0.01)
print(cooccur.ponds_0.01)



ponds_0.05 <- read.csv("../Data/0.05_cooccur_data.csv", 
                       row.names=1)


cooccur.ponds_0.05 <- cooccur(mat = ponds_0.05,     # the dataframe
                              type = "spp_site",  # data organized with species as rows and sites as column
                              thresh = TRUE,      # summarize the most important species associations
                              spp_names = TRUE)   # use species names in data

class(cooccur.ponds_0.05)
summary(cooccur.ponds_0.05)
prob.table(cooccur.ponds_0.05)
print(cooccur.ponds_0.05)



ponds_0.1 <- read.csv("../Data/0.1_cooccur_data.csv", 
                       row.names=1)


cooccur.ponds_0.1 <- cooccur(mat = ponds_0.1,     # the dataframe
                              type = "spp_site",  # data organized with species as rows and sites as column
                              thresh = TRUE,      # summarize the most important species associations
                              spp_names = TRUE)   # use species names in data

class(cooccur.ponds_0.1)
summary(cooccur.ponds_0.1)
prob.table(cooccur.ponds_0.1)
print(cooccur.ponds_0.1)



ponds_0.3 <- read.csv("../Data/0.3_cooccur_data.csv", 
                      row.names=1)


cooccur.ponds_0.3 <- cooccur(mat = ponds_0.3,     # the dataframe
                             type = "spp_site",  # data organized with species as rows and sites as column
                             thresh = TRUE,      # summarize the most important species associations
                             spp_names = TRUE)   # use species names in data

class(cooccur.ponds_0.3)
summary(cooccur.ponds_0.3)
prob.table(cooccur.ponds_0.3)
print(cooccur.ponds_0.3)



## Plot species co-occurence at different false positive thresholds

p10a <- plot(cooccur.ponds_NT, main = "a")
p10b <- plot(cooccur.ponds_0.0005)
p10c <- plot(cooccur.ponds_0.001)
p10d <- plot(cooccur.ponds_0.005)
p10e <- plot(cooccur.ponds_0.01)
p10f <- plot(cooccur.ponds_0.05)
p10g <- plot(cooccur.ponds_0.1)
p10h <- plot(cooccur.ponds_0.3)
p10 <- plot(cooccur.ponds_SS)


grid.arrange(arrangeGrob(p10a,p10b, nrow=1, ncol=2),
             arrangeGrob(p10c,p10d, nrow=1, ncol=2),
             arrangeGrob(p10e,p10f, nrow=1, ncol=2),
             arrangeGrob(p10g,p10h,p10, nrow=1, ncol=3), 
             nrow=4,
             top=textGrob("Effect of false positive sequence threshold on species associations",
                          gp=gpar(cex=2)))



## Inspect great crested newt interactions at each threshold
pair(mod = cooccur.ponds_NT, "Triturus cristatus")
pair(mod = cooccur.ponds_0.0005, "Triturus cristatus")
pair(mod = cooccur.ponds_0.001, "Triturus cristatus")
pair(mod = cooccur.ponds_0.005, "Triturus cristatus")
pair(mod = cooccur.ponds_0.01, "Triturus cristatus")
pair(mod = cooccur.ponds_0.05, "Triturus cristatus")
pair(mod = cooccur.ponds_0.1, "Triturus cristatus")
pair(mod = cooccur.ponds_0.3, "Triturus cristatus")
pair(mod = cooccur.ponds_SS, "Triturus cristatus")



## Construct pair profiles at different thresholds

p11a <- pair.profile(cooccur.ponds_NT)
p11b <- pair.profile(cooccur.ponds_0.0005)
p11c <- pair.profile(cooccur.ponds_0.001)
p11d <- pair.profile(cooccur.ponds_0.005)
p11e <- pair.profile(cooccur.ponds_0.01)
p11f <- pair.profile(cooccur.ponds_0.05)
p11g <- pair.profile(cooccur.ponds_0.1)
p11h <- pair.profile(cooccur.ponds_0.3)
p11 <- pair.profile(cooccur.ponds_SS)


grid.arrange(arrangeGrob(p11a,p11b, nrow=1, ncol=2),
             arrangeGrob(p11c,p11d, nrow=1, ncol=2),
             arrangeGrob(p11e,p11f, nrow=1, ncol=2),
             arrangeGrob(p11g,p11h,p11, nrow=1, ncol=3), 
             nrow=4,
             top=textGrob("Effect of false positive sequence thresholds on pairing plots",
                          gp=gpar(cex=2)))



## Plot degree to which species pairs deviate from expected co-occurence


p12a <- obs.v.exp(cooccur.ponds_NT)
p12b <- obs.v.exp(cooccur.ponds_0.0005)
p12c <- obs.v.exp(cooccur.ponds_0.001)
p12d <- obs.v.exp(cooccur.ponds_0.005)
p12e <- obs.v.exp(cooccur.ponds_0.01)
p12f <- obs.v.exp(cooccur.ponds_0.05)
p12g <- obs.v.exp(cooccur.ponds_0.1)
p12h <- obs.v.exp(cooccur.ponds_0.3)
p12 <- obs.v.exp(cooccur.ponds_SS)


grid.arrange(arrangeGrob(p12a,p12b, nrow=1, ncol=2),
             arrangeGrob(p12c,p12d, nrow=1, ncol=2),
             arrangeGrob(p12e,p12f, nrow=1, ncol=2),
             arrangeGrob(p12g,p12h,p12, nrow=1, ncol=3), 
             nrow=4,
             top=textGrob("Effect of false positive sequence thresholds on expected Vs observed co-occurence",
                          gp=gpar(cex=2)))



#'
#' Effect of sequence thresholds on environmental predictors of GCN
#' 


## Make new data frames
spp_NT <- no_threshold[,-c(62:66)]
spp_0.0005 <- threshold_0.0005[,-c(62:66)]
spp_0.001 <- threshold_0.001[,-c(62:66)]
spp_0.005 <- threshold_0.005[,-c(62:66)]
spp_0.01 <- threshold_0.01[,-c(62:66)]
spp_0.05 <- threshold_0.05[,-c(62:66)]
spp_0.1 <- threshold_0.1[,-c(62:66)]
spp_0.3 <- threshold_0.3[,-c(62:66)]
spp_SS <- spp_threshold[,-c(62:66)]


## Create column containing the total number of species detected in
## each sample i.e. species richness
spp_NT$Sp_richness <- rowSums(spp_NT[,2:61])
spp_0.0005$Sp_richness <- rowSums(spp_0.0005[,2:61])
spp_0.001$Sp_richness <- rowSums(spp_0.001[,2:61])
spp_0.005$Sp_richness <- rowSums(spp_0.005[,2:61])
spp_0.01$Sp_richness <- rowSums(spp_0.01[,2:61])
spp_0.05$Sp_richness <- rowSums(spp_0.05[,2:61])
spp_0.1$Sp_richness <- rowSums(spp_0.1[,2:61])
spp_0.3$Sp_richness <- rowSums(spp_0.3[,2:61])
spp_SS$Sp_richness <- rowSums(spp_SS[,2:61])


## Import habitat metadata for eDNA samples
metadata <- read.csv("../Data/HabitatMetadata.csv", header=TRUE)


## merge the habitat metadata with the metabarcoding data using the 'sample' 
## column from the metadata dataframe and the 'sample' column from the 
## metabarcoding dataframes

habitat_NT <- merge(metadata, spp_NT, all=TRUE)
habitat_0.0005 <- merge(metadata, spp_0.0005, all=TRUE)
habitat_0.001 <- merge(metadata, spp_0.001, all=TRUE)
habitat_0.005 <- merge(metadata, spp_0.005, all=TRUE)
habitat_0.01 <- merge(metadata, spp_0.01, all=TRUE)
habitat_0.05 <- merge(metadata, spp_0.05, all=TRUE)
habitat_0.1 <- merge(metadata, spp_0.1, all=TRUE)
habitat_0.3 <- merge(metadata, spp_0.3, all=TRUE)
habitat_SS <- merge(metadata, spp_SS, all=TRUE)


## Remove samples which did not have associated pond metadata
## (ADAS UK Ltd samples and 4 samples from FERA)
habitat_NT <- habitat_NT[-c(505:532),]
habitat_0.0005 <- habitat_0.0005[-c(505:532),]
habitat_0.001 <- habitat_0.001[-c(505:532),]
habitat_0.005 <- habitat_0.005[-c(505:532),]
habitat_0.01 <- habitat_0.01[-c(505:532),]
habitat_0.05 <- habitat_0.05[-c(505:532),]
habitat_0.1 <- habitat_0.1[-c(505:532),]
habitat_0.3 <- habitat_0.3[-c(505:532),]
habitat_SS <- habitat_SS[-c(505:532),]


#'
#' Now identify a suitable parsimonious model of great crested newt
#' presence-absence at difference false positive sequence thresholds.
#' 
#' 
#' Exploratory analysis that dealt with collinearity has already been 
#' performed for the species-specific false positive sequence threshold 
#' analysis. So we need only examine relative importance of explanatory 
#' variables at different thresholds. 
#' 
#' Analysis of spatial autocorrelation has also already been performed
#' in the species-specific false positive sequence threshold analysis.
#' 
#' As such, we know a Generalised Linear Mixed Model is required to analyse
#' the data in order to deal with dependencies within sites through the
#' introduction of random effects, thus controlling the influence of 
#' spatial autocorrelation on the analysis.
#' 
#' We also know that a single model containing all explanatory variables 
#' was a better fit to the data than functional models containing subsets
#' of the explanatory variables, based on AIC value.
#' 
#' Therefore, we can now go straight to model selection for the data at
#' each different false positive sequence threshold. We don't need to re-run
#' the analysis from the species-specific threshold as the results from 
#' each threshold dataframe will be plotted separately.
#' 


################
# NO THRESHOLD #
################

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bos_taurus + Lissotriton_vulgaris + Cyprinus_carpio
              + Fulica_atra + Sus_scrofa + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


GCN_tree <- rpart(f1, data = habitat_NT, method = "class",
                  minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - L. vulgaris presence
#' - pond substrate
#' - species richness
#' - terrestrial other
#' - inflow
#' - pond area
#' - common moorhen presence
#' - macrophytes
#' - fish presence
#' - pond density
#' - ruderals
#' - cow presence
#' - woodland
#' - water quality
#' - waterfowl presence
#' - outflow
#' - overhang
#' - pig presence
#' - max depth
#' - rough grass
#' - common carp presence
#' 

## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 8-10 is optimal i.e. 8-10 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - amphibian presence
#' - permanence
#' - pollution
#' - scrub/hedge
#' - overall terrestrial habitat quality
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist_NT <- data.frame(habitat_NT[,c(1,10,14:16,18,21:24,27:30,32:34,
                                       43,51,57,65,94:95,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

GCNdist_NT < scale(GCNdist_NT[,c(2:6,23)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Pond_Substrate + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Macrophytes + Fish + Pond_Density 
                + Ruderals + Bos_taurus + Woodland + Quality
                + Waterfowl + Outflow + Overhang + Sus_scrofa
                + Max_Depth + Rough_Grass + Terrestrial_Other
                + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)


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
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Macrophytes + Fish + Pond_Density 
                + Ruderals + Bos_taurus + Woodland + Quality
                + Waterfowl + Outflow + Overhang + Sus_scrofa
                + Max_Depth + Rough_Grass + Terrestrial_Other
                + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove fish presence as this variable has the lowest AIC value

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Macrophytes + Pond_Density 
                + Ruderals + Bos_taurus + Woodland + Quality
                + Waterfowl + Outflow + Overhang + Sus_scrofa
                + Max_Depth + Rough_Grass + Terrestrial_Other
                + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish presence has resulted in better model fit


drop1(model2)


## Remove water quality next as this variable has the lowest AIC value

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Macrophytes + Pond_Density 
                + Ruderals + Bos_taurus + Woodland
                + Waterfowl + Outflow + Overhang + Sus_scrofa
                + Max_Depth + Rough_Grass + Terrestrial_Other
                + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of water quality has resulted in better model fit


drop1(model3)


## Remove rough grass as this now has lowest AIC value

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Macrophytes + Pond_Density 
                + Ruderals + Bos_taurus + Woodland
                + Waterfowl + Outflow + Overhang + Sus_scrofa
                + Max_Depth + Terrestrial_Other + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of rough grass has resulted in better model fit


drop1(model4)


## Remove woodland as this now has lowest AIC value

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Macrophytes + Pond_Density 
                + Ruderals + Bos_taurus
                + Waterfowl + Outflow + Overhang + Sus_scrofa
                + Max_Depth + Terrestrial_Other + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit


drop1(model5)


## Remove pig as this now has lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Macrophytes + Pond_Density 
                + Ruderals + Bos_taurus
                + Waterfowl + Outflow + Overhang
                + Max_Depth + Terrestrial_Other + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of pig has resulted in better model fit


drop1(model6)


## Remove macrophytes as this now has lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Pond_Density + Ruderals + Bos_taurus 
                + Waterfowl + Outflow + Overhang
                + Max_Depth + Terrestrial_Other + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of macrophytes has resulted in better model fit


drop1(model7)


## Remove pond density as this now has lowest AIC value

model8 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Ruderals + Bos_taurus 
                + Waterfowl + Outflow + Overhang
                + Max_Depth + Terrestrial_Other + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond density has resulted in better model fit


drop1(model8)


## Remove waterfowl presence as this now has lowest AIC value

model9 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Pond_Area + Gallinula_chloropus
                + Inflow + Ruderals + Bos_taurus 
                + Outflow + Overhang
                + Max_Depth + Terrestrial_Other + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model9,model8)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl presence has resulted in better model fit


drop1(model9)


## Remove pond area as this now has lowest AIC value

model10 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Terrestrial_Other + Gallinula_chloropus
                + Inflow + Ruderals + Bos_taurus 
                + Outflow + Overhang
                + Max_Depth + Terrestrial_Other + Cyprinus_carpio,
                family = binomial,
                data = GCNdist_NT)

anova(model10,model9)

## p-value from LRT not significant and reduction in AIC value (<2)
## In this instance, it is better to choose model with lower AIC and
## less parameters in order to avoid overparameterisation
## Removal of pond area has resulted in better model fit


drop1(model10)


## Remove terrestrial other as this now has lowest AIC value

model11 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Gallinula_chloropus
                 + Inflow + Ruderals + Bos_taurus 
                 + Outflow + Overhang
                 + Max_Depth + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

anova(model11,model10)

## p-value from LRT not significant and reduction in AIC value (<2)
## In this instance, it is better to choose model with lower AIC and
## less parameters in order to avoid overparameterisation
## Removal of terrestrial other has resulted in better model fit


drop1(model11)


## Remove overhang as this now has lowest AIC value

model12 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Gallinula_chloropus
                 + Inflow + Ruderals + Bos_taurus 
                 + Outflow + Max_Depth + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

anova(model12,model11)


## p-value from LRT not significant and reduction in AIC value (<2)
## In this instance, it is better to choose model with lower AIC and
## less parameters in order to avoid overparameterisation
## Removal of overhang has resulted in better model fit


drop1(model12)


## Remove cow as this now has lowest AIC value

model13 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Gallinula_chloropus
                 + Inflow + Ruderals
                 + Outflow + Max_Depth + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

anova(model13,model12)

## p-value from LRT not significant and reduction in AIC value (<2)
## In this instance, it is better to choose model with lower AIC and
## less parameters in order to avoid overparameterisation
## Removal of cow presence has resulted in better model fit


drop1(model13)


## Remove moorhen as this now has lowest AIC value

model14 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness + Inflow 
                 + Ruderals + Outflow + Max_Depth + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

anova(model14,model13)

## p-value from LRT not significant but minimal reduction in AIC value 
## (<2). It is better to choose model with lower AIC and
## less parameters in order to avoid overparameterisation
## Removal of moorhen presence has resulted in better model fit


drop1(model14)


## Remove outflow as this now has lowest AIC value

model15 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness + Inflow 
                 + Ruderals + Max_Depth + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

anova(model15,model14)

## p-value from LRT not significant but minimal increase in AIC value 
## (<2). It is better to choose model with less parameters in order 
## to avoid overparameterisation
## Removal of outflow has resulted in better model fit


drop1(model15)


## Remove max depth as this now has lowest AIC value

model16 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness + Inflow 
                 + Ruderals + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

anova(model16,model15)

## p-value from LRT not significant but minimal increase in AIC value 
## (<2). It is better to choose model with less parameters in order 
## to avoid overparameterisation
## Removal of max depth has resulted in better model fit


drop1(model16)


## Remove ruderals as this now has lowest AIC value

model17 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness + Inflow 
                 + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

anova(model17,model16)

## p-value from LRT significant and increase in AIC value 
## (>2). Removal of ruderals has reduced model fit.


## Final model: 

GCNm_NT <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness + Inflow 
                 + Ruderals + Cyprinus_carpio,
                 family = binomial,
                 data = GCNdist_NT)

summary(GCNm_NT)
anova(GCNm_NT)
drop1(GCNm_NT, test = "Chi")   # for binomial/integer y models, 
                               # statistics for significance of each term

display(GCNm_NT)
se.ranef(GCNm_NT)              # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(GCNm_NT)

## Overdispersion test (chi-square)

1-pchisq(518.488, df=496)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(GCNm_NT)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(GCNm_NT, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(GCNm_NT)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(GCNm_NT)


## Plot the residuals against each x vairable
plot(sresid ~ GCNdist_NT$Lissotriton_vulgaris)
plot(sresid ~ GCNdist_NT$Sp_richness)
plot(sresid ~ GCNdist_NT$Inflow)
plot(sresid ~ GCNdist_NT$Ruderals)
plot(sresid ~ GCNdist_NT$Cyprinus_carpio)


## Check model residuals
chkres(GCNm_NT)

## Plot the fitted data against the observed data
plot(GCNdist_NT$Triturus_cristatus ~ fitted(GCNm_NT))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(GCNm_NT)


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

hoslem.test(GCNdist_NT$Triturus_cristatus, fitted(GCNm_NT))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(GCNm_NT, newdata=GCNdist_NT, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(GCNdist_NT, fit)


## PLOT 13A: Smooth newt presence and GCN presence/absence 

p13a <- ggplot(box.dat) + ggtitle('(a)')
p13a <- p13a + geom_jitter(aes(x=factor(Lissotriton_vulgaris), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p13a <- p13a + geom_boxplot(aes(x=factor(Lissotriton_vulgaris), y=fit), alpha=0.7, outlier.colour = "blue")
p13a <- p13a + scale_colour_gradient(low="grey48", high="orange")
p13a <- p13a + scale_x_discrete(breaks=c("0", "1"),
                              labels=c("Absent", "Present"))
p13a <- p13a + labs(y="Probability of GCN presence", x = "Smooth newt") 
p13a <- p13a + theme(panel.background = element_rect(fill = 'grey98'),
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position="none")
p13a



## PLOT 13B: Species richness and GCN presence/absence

## Create a range of species richness values which increase by 0.02385
## to predict GCN occupancy as species richness increases
range <- seq(from=min(GCNdist_NT$Sp_richness), to=max(GCNdist_NT$Sp_richness), by=0.02385)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Lissotriton_vulgaris=rep(0, length(range)),
                 Sp_richness=range,
                 Inflow=rep('Present', length(range)),
                 Ruderals=rep('Important', length(range)),
                 Cyprinus_carpio=rep(0, length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_NT, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.spp <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p13b <- ggplot() + ggtitle('(b)')
p13b <- p13b + geom_ribbon(aes(x=dat.spp$Sp_richness, ymin = dat.spp$cil, ymax = dat.spp$ciu), alpha = 0.25)
p13b <- p13b + geom_line(aes(x=dat.spp$Sp_richness, y=dat.spp$fit), size = 1)
p13b <- p13b + geom_jitter(aes(x=GCNdist_NT$Sp_richness, y=GCNdist_NT$Triturus_cristatus, colour=GCNdist_NT$Triturus_cristatus), cex=1, height=0.7)
p13b <- p13b + scale_colour_gradient(low="grey48", high="orange")
p13b <- p13b + labs(y = "Probability of GCN presence", x = "Species richness")
p13b <- p13b + coord_cartesian(xlim=c(0,12))
p13b <- p13b + theme(panel.background = element_rect(fill = 'grey98'), 
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position="none")
p13b



## PLOT 13C: inflow and GCN presence/absence

p13c <- ggplot(box.dat) + ggtitle('(c)')
p13c <- p13c + geom_jitter(aes(x=Inflow, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p13c <- p13c + geom_boxplot(aes(x=Inflow, y=fit), alpha=0.7, outlier.colour = "blue")
p13c <- p13c + scale_colour_gradient(low="grey48", high="orange")
p13c <- p13c + labs(y="Probability of GCN presence", x = "Pond inflow") 
p13c <- p13c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p13c



## PLOT 13D: ruderals and GCN presence/absence

p13d <- ggplot(box.dat) + ggtitle('(d)')
p13d <- p13d + geom_jitter(aes(x=Ruderals, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p13d <- p13d + geom_boxplot(aes(x=Ruderals, y=fit), alpha=0.7, outlier.colour = "blue")
p13d <- p13d + scale_colour_gradient(low="grey48", high="orange")
p13d <- p13d + labs(y="Probability of GCN presence", x = "Ruderals") 
p13d <- p13d + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p13d



## PLOT 13E: common carp and GCN presence/absence

p13e <- ggplot(box.dat) + ggtitle('(e)')
p13e <- p13e + geom_jitter(aes(x=factor(Cyprinus_carpio), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p13e <- p13e + geom_boxplot(aes(x=factor(Cyprinus_carpio), y=fit), alpha=0.7, outlier.colour = "blue")
p13e <- p13e + scale_colour_gradient(low="grey48", high="orange")
p13e <- p13e + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p13e <- p13e + labs(y="Probability of GCN presence", x = "Common carp") 
p13e <- p13e + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p13e


grid.arrange(arrangeGrob(p13a,p13b,p13c, ncol=3),
             arrangeGrob(p13d,p13e, ncol=2), nrow=2,
             top=textGrob("Predictors of great crested newt (GCN) presence with no sequence threshold",
                          gp=gpar(cex=3)))




#########
# 0.05% #
#########

## Classification tree

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris 
              + Sciurus_carolinensis + Gasterosteus_aculeatus 
              + Fulica_atra + Sus_scrofa + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Terrestrial_Habitat)

GCN_tree <- rpart(f1, data = habitat_0.0005, method = "class",
                  minsplit=5, cp = 0.001)

par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - L. vulgaris presence
#' - pond substrate
#' - species richness
#' - waterfowl presence
#' - B. bufo presence
#' - pond density
#' - woodland
#' - S. carolinensis presence
#' - G. aculeatus presence
#' - inflow
#' - macrophytes
#' - terrestrial other
#' - water quality
#' - pond area
#' - permanence
#' - outflow
#' - scrub/hedge
#' - max depth
#' - ruderals
#' - rough grass
#' - overhang
#' - F. atra presence
#' 

## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 7 is optimal i.e. 7 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - amphibian presence
#' - fish presence
#' - pollution
#' - overall terrestrial habitat quality
#' 


## Create new dataframe containing only explanatory variables to be
## modelled
GCNdist_0.0005 <- data.frame(habitat_0.0005[,c(1,10,14:16,18,20:24,27,
                                               29:33,44,56,59,65,88,
                                               95,97)])

## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.
GCNdist_0.0005 < scale(GCNdist_0.0005[,c(2:6,24)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Pond_Substrate + Sp_richness
                + Waterfowl + Bufo_bufo + Pond_Density + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Quality + Pond_Area + Permanance + Outflow
                + Scrub_Hedge + Max_Depth + Ruderals + Rough_Grass
                + Overhang + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)


## drop1 function:

drop1(model0)


## Remove pond substrate as this variable has the lowest AIC value

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Pond_Density + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Quality + Pond_Area + Permanance + Outflow
                + Scrub_Hedge + Max_Depth + Ruderals + Rough_Grass
                + Overhang + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove rough grass as this variable now has the lowest AIC value

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Pond_Density + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Quality + Pond_Area + Permanance + Outflow
                + Scrub_Hedge + Max_Depth + Ruderals
                + Overhang + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of rough grass  has resulted in better model fit


drop1(model2)


## Remove water quality as this variable now has the lowest AIC value

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Pond_Density + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Pond_Area + Permanance + Outflow
                + Scrub_Hedge + Max_Depth + Ruderals
                + Overhang + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of water quality has resulted in better model fit


drop1(model3)


## Remove scrub/hedge as this variable now has the lowest AIC value

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Pond_Density + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Pond_Area + Permanance + Outflow
                + Max_Depth + Ruderals
                + Overhang + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit


drop1(model4)


## Remove pond density as this variable now has the lowest AIC value

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Pond_Area + Permanance + Outflow
                + Max_Depth + Ruderals
                + Overhang + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond density has resulted in better model fit


drop1(model5)


## Remove outflow as this variable now has the lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Pond_Area + Permanance
                + Max_Depth + Ruderals
                + Overhang + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of outflow has resulted in better model fit


drop1(model6)


## Remove overhang as this variable now has the lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Pond_Area + Permanance
                + Max_Depth + Ruderals + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model7,model6)

## p-value from LRT not significant and minimual reduction in AIC value
## (<2) therefore it is better to choose model with lower AIC and less
## parameters to avoid overparameterisation
## removal of overhang has resulted in better model fit


drop1(model7)


## Remove max depth as this variable now has the lowest AIC value

model8 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Macrophytes + Terrestrial_Other
                + Pond_Area + Permanance
                + Ruderals + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model8,model7)

## p-value from LRT not significant and minimal difference in AIC value
## (<2) therefore it is better to choose model with lower AIC and less
## parameters to avoid overparameterisation
## removal of max depth has resulted in better model fit


drop1(model8)


## Remove macrophytes as this variable now has the lowest AIC value

model9 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Terrestrial_Other
                + Pond_Area + Permanance
                + Ruderals + Fulica_atra,
                family = binomial,
                data = GCNdist_0.0005)

anova(model9,model8)

## p-value from LRT not significant and minimal difference in AIC value
## (<2) therefore it is better to choose model with lower AIC and less
## parameters to avoid overparameterisation
## removal of macrophytes has resulted in better model fit


drop1(model9)


## Remove coot presence as this variable now has the lowest AIC value

model10 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness
                + Waterfowl + Bufo_bufo + Woodland
                + Sciurus_carolinensis + Gasterosteus_aculeatus
                + Inflow + Terrestrial_Other
                + Pond_Area + Permanance + Ruderals,
                family = binomial,
                data = GCNdist_0.0005)

anova(model10,model9)

## p-value from LRT not significant and minimal difference in AIC value
## (<2) therefore it is better to choose model with lower AIC and less
## parameters to avoid overparameterisation
## removal of coot presence has resulted in better model fit


drop1(model10)


## Remove waterfowl as this variable now has the lowest AIC value

model11 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness
                 + Bufo_bufo + Woodland
                 + Sciurus_carolinensis + Gasterosteus_aculeatus
                 + Inflow + Terrestrial_Other
                 + Pond_Area + Permanance + Ruderals,
                 family = binomial,
                 data = GCNdist_0.0005)

anova(model11,model10)

## p-value from LRT not significant and minimal difference in AIC value
## (<2) therefore it is better to choose model with lower AIC and less
## parameters to avoid overparameterisation
## removal of waterfowl has resulted in better model fit


drop1(model11)


## Remove terrestrial other as this variable now has the lowest AIC value

model12 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness
                 + Bufo_bufo + Woodland
                 + Sciurus_carolinensis + Gasterosteus_aculeatus
                 + Inflow + Pond_Area + Permanance + Ruderals,
                 family = binomial,
                 data = GCNdist_0.0005)

anova(model12,model11)

## p-value from LRT not significant and minimal increase in AIC value
## (<2) therefore it is better to choose model with less parameters 
## to avoid overparameterisation
## removal of terrestrial other has resulted in better model fit


drop1(model12)


## Remove woodland as this variable now has the lowest AIC value

model13 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness
                 + Bufo_bufo + Sciurus_carolinensis + Gasterosteus_aculeatus
                 + Inflow + Pond_Area + Permanance + Ruderals,
                 family = binomial,
                 data = GCNdist_0.0005)

anova(model13,model12)

## p-value from LRT not significant and minimal increase in AIC value
## (<2) therefore it is better to choose model with less parameters 
## to avoid overparameterisation
## removal of woodland has resulted in better model fit


drop1(model13)


## Remove permanence as this variable now has the lowest AIC value

model14 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness
                 + Bufo_bufo + Sciurus_carolinensis + Gasterosteus_aculeatus
                 + Inflow + Pond_Area + Ruderals,
                 family = binomial,
                 data = GCNdist_0.0005)

anova(model14,model13)

## p-value from LRT significant and increase in AIC value (>2) 
## removal of permanence has reduced model fit


## Final model: 

GCNm_0.0005 <- glmer(Triturus_cristatus ~ (1|sample)
                       + Lissotriton_vulgaris + Sp_richness
                       + Bufo_bufo + Sciurus_carolinensis + Gasterosteus_aculeatus
                       + Inflow + Pond_Area + Permanance + Ruderals,
                       family = binomial,
                       data = GCNdist_0.0005)

summary(GCNm_0.0005)
anova(GCNm_0.0005)
drop1(GCNm_0.0005, test = "Chi")   # for binomial/integer y models, 
                                   # statistics for significance of each term

display(GCNm_0.0005)
se.ranef(GCNm_0.0005)              # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(GCNm_0.0005)

## Overdispersion test (chi-square)

1-pchisq(422.829, df=490)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(GCNm_0.0005)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(GCNm_0.0005, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(GCNm_0.0005)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(GCNm_0.0005)


## Plot the residuals against each x vairable
plot(sresid ~ GCNdist_0.0005$Lissotriton_vulgaris)
plot(sresid ~ GCNdist_0.0005$Sp_richness)
plot(sresid ~ GCNdist_0.0005$Bufo_bufo)
plot(sresid ~ GCNdist_0.0005$Sciurus_carolinensis)
plot(sresid ~ GCNdist_0.0005$Gasterosteus_aculeatus)
plot(sresid ~ GCNdist_0.0005$Inflow)
plot(sresid ~ GCNdist_0.0005$Pond_Area)
plot(sresid ~ GCNdist_0.0005$Permanance)
plot(sresid ~ GCNdist_0.0005$Ruderals)


## Check model residuals
chkres(GCNm_0.0005)

## Plot the fitted data against the observed data
plot(GCNdist_0.0005$Triturus_cristatus ~ fitted(GCNm_0.0005))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(GCNm_0.0005)


## Hosmer and Lemeshow Goodness of Fit Test

hoslem.test(GCNdist_0.0005$Triturus_cristatus, fitted(GCNm_0.0005))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(GCNm_0.0005, newdata=GCNdist_0.0005, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(GCNdist_0.0005, fit)



## PLOT 14A: Smooth newt presence and GCN presence/absence 

p14a <- ggplot(box.dat) + ggtitle('(a)')
p14a <- p14a + geom_jitter(aes(x=factor(Lissotriton_vulgaris), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p14a <- p14a + geom_boxplot(aes(x=factor(Lissotriton_vulgaris), y=fit), alpha=0.7, outlier.colour = "blue")
p14a <- p14a + scale_colour_gradient(low="grey48", high="orange")
p14a <- p14a + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p14a <- p14a + labs(y="Probability of GCN presence", x = "Smooth newt") 
p14a <- p14a + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14a



## PLOT 14B: Species richness and GCN presence/absence

## Create a range of species richness values which increase by 0.02185
## to predict GCN occupancy as species richness increases
range <- seq(from=min(GCNdist_0.0005$Sp_richness), to=max(GCNdist_0.0005$Sp_richness), by=0.02185)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Lissotriton_vulgaris=rep(0, length(range)),
                 Sp_richness=range,
                 Bufo_bufo=rep(0, length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Permanance=rep('Never dries', length(range)),
                 Ruderals=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.0005, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.spp <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p14b <- ggplot() + ggtitle('(b)')
p14b <- p14b + geom_ribbon(aes(x=dat.spp$Sp_richness, ymin = dat.spp$cil, ymax = dat.spp$ciu), alpha = 0.25)
p14b <- p14b + geom_line(aes(x=dat.spp$Sp_richness, y=dat.spp$fit), size = 1)
p14b <- p14b + geom_jitter(aes(x=GCNdist_0.0005$Sp_richness, y=GCNdist_0.0005$Triturus_cristatus, colour=GCNdist_0.0005$Triturus_cristatus), cex=1, height=0.7)
p14b <- p14b + scale_colour_gradient(low="grey48", high="orange")
p14b <- p14b + labs(y = "Probability of GCN presence", x = "Species richness")
p14b <- p14b + coord_cartesian(xlim=c(0,12))
p14b <- p14b + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14b



## PLOT 14C: common toad presence and GCN presence/absence 

p14c <- ggplot(box.dat) + ggtitle('(c)')
p14c <- p14c + geom_jitter(aes(x=factor(Bufo_bufo), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p14c <- p14c + geom_boxplot(aes(x=factor(Bufo_bufo), y=fit), alpha=0.7, outlier.colour = "blue")
p14c <- p14c + scale_colour_gradient(low="grey48", high="orange")
p14c <- p14c + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p14c <- p14c + labs(y="Probability of GCN presence", x = "Common toad") 
p14c <- p14c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14c



## PLOT 14D: grey squirrel presence and GCN presence/absence 

p14d <- ggplot(box.dat) + ggtitle('(d)')
p14d <- p14d + geom_jitter(aes(x=factor(Sciurus_carolinensis), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p14d <- p14d + geom_boxplot(aes(x=factor(Sciurus_carolinensis), y=fit), alpha=0.7, outlier.colour = "blue")
p14d <- p14d + scale_colour_gradient(low="grey48", high="orange")
p14d <- p14d + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p14d <- p14d + labs(y="Probability of GCN presence", x = "Grey squirrel") 
p14d <- p14d + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14d



## PLOT 14E: three-spined stickleback presence and GCN presence/absence 

p14e <- ggplot(box.dat) + ggtitle('(e)')
p14e <- p14e + geom_jitter(aes(x=factor(Gasterosteus_aculeatus), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p14e <- p14e + geom_boxplot(aes(x=factor(Gasterosteus_aculeatus), y=fit), alpha=0.7, outlier.colour = "blue")
p14e <- p14e + scale_colour_gradient(low="grey48", high="orange")
p14e <- p14e + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p14e <- p14e + labs(y="Probability of GCN presence", x = "Three-spined stickleback") 
p14e <- p14e + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14e



## PLOT 14F: inflow and GCN presence/absence

p14f <- ggplot(box.dat) + ggtitle('(f)')
p14f <- p14f + geom_jitter(aes(x=Inflow, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p14f <- p14f + geom_boxplot(aes(x=Inflow, y=fit), alpha=0.7, outlier.colour = "blue")
p14f <- p14f + scale_colour_gradient(low="grey48", high="orange")
p14f <- p14f + labs(y="Probability of GCN presence", x = "Pond inflow") 
p14f <- p14f + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14f



## PLOT 14G: Pond area and GCN presence/absence

## Create a range of pond area values which increase by 18.6012
## to predict GCN occupancy as pond area increases
range <- seq(from=min(GCNdist_0.0005$Pond_Area), to=max(GCNdist_0.0005$Pond_Area), by=18.6012)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Lissotriton_vulgaris=rep(0, length(range)),
                 Sp_richness=rep(3.551587, length(range)),
                 Bufo_bufo=rep(0, length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Pond_Area=range,
                 Permanance=rep('Never dries', length(range)),
                 Ruderals=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.0005, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.PA <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p14g <- ggplot() + ggtitle('(g)')
p14g <- p14g + geom_ribbon(aes(x=dat.PA$Pond_Area, ymin = dat.PA$cil, ymax = dat.PA$ciu), alpha = 0.25)
p14g <- p14g + geom_line(aes(x=dat.PA$Pond_Area, y=dat.PA$fit), size = 1)
p14g <- p14g + geom_jitter(aes(x=GCNdist_0.0005$Pond_Area, y=GCNdist_0.0005$Triturus_cristatus, colour=GCNdist_0.0005$Triturus_cristatus), cex=1, height=0.7)
p14g <- p14g + scale_colour_gradient(low="grey48", high="orange")
p14g <- p14g + labs(y = "Probability of GCN presence", x = "Pond area ("~m^2~")")
p14g <- p14g + coord_cartesian(xlim=c(0,10000))
p14g <- p14g + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14g



## PLOT 14H: permanence and GCN presence/absence

f.Permanence <- ordered(GCNdist_0.0005$Permanance, levels = c("Never dries", "Rarely dries","Sometimes dries","Dries annually"))
f.Permanence

p14h <- ggplot(box.dat) + ggtitle('(h)')
p14h <- p14h + geom_jitter(aes(x=f.Permanence, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p14h <- p14h + geom_boxplot(aes(x=f.Permanence, y=fit), alpha=0.7, outlier.colour = "blue")
p14h <- p14h + scale_colour_gradient(low="grey48", high="orange")
p14h <- p14h + labs(y="Probability of GCN presence", x = "Pond permanence") 
p14h <- p14h + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14h



## PLOT 14I: ruderals and GCN presence/absence

f.Ruderals <- ordered(GCNdist_0.0005$Ruderals, levels = c("None", "Some","Important"))
f.Ruderals

p14i <- ggplot(box.dat) + ggtitle('(i)')
p14i <- p14i + geom_jitter(aes(x=f.Ruderals, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p14i <- p14i + geom_boxplot(aes(x=f.Ruderals, y=fit), alpha=0.7, outlier.colour = "blue")
p14i <- p14i + scale_colour_gradient(low="grey48", high="orange")
p14i <- p14i + labs(y="Probability of GCN presence", x = "Ruderals") 
p14i <- p14i + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p14i


grid.arrange(arrangeGrob(p14a,p14b,p14c, ncol=3),
             arrangeGrob(p14d,p14e,p14f, ncol=3),
             arrangeGrob(p14g,p14h,p14i, ncol=3), nrow=3,
             top=textGrob("Predictors of great crested newt (GCN) presence with 0.05% sequence threshold",
                          gp=gpar(cex=3)))



#########
# 0.1%  #
#########

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris + Sciurus_carolinensis
              + Gasterosteus_aculeatus + Fulica_atra + Sus_scrofa + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


GCN_tree <- rpart(f1, data = habitat_0.001, method = "class",
                  minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - L. vulgaris presence
#' - pond substrate
#' - species richness
#' - amphibian presence
#' - pond density
#' - B. bufo presence
#' - waterfowl presence
#' - inflow
#' - S.carolinensis presence
#' - max depth
#' - G. aculeatus presence
#' - permanence
#' - macrophytes
#' - overhang
#' - pond area
#' - ruderals
#' - woodland
#' - rough grass
#' - water quality
#' - G. chloropus presence
#' - terrestrial other
#' - overall terrestrial habitat quality
#' - scrub/hedge
#' - S. scrofa presence
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 7 or 21 is optimal i.e. 7 or 21 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - fish presence
#' - pollution
#' - outflow
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist_0.001 <- data.frame(habitat_0.001[,c(1,10,14:16,18,20:23,26:27,
                                             29:34,44,57,59,65,88,94:95,
                                             97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

GCNdist_0.001 < scale(GCNdist_0.001[,c(2:6,26)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Pond_Substrate + Sp_richness 
                + Amphibians + Pond_Density + Bufo_bufo + Waterfowl
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Rough_Grass + Quality + Gallinula_chloropus
                + Terrestrial_Other + Scrub_Hedge 
                + Sus_scrofa + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.001)


## drop1 function:

drop1(model0)


## Remove pond substrate  as this variable has the lowest AIC value

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Pond_Density + Bufo_bufo + Waterfowl
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Rough_Grass + Quality + Gallinula_chloropus
                + Terrestrial_Other + Scrub_Hedge 
                + Sus_scrofa + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.001)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove amphibian presence as this variable has the lowest AIC value

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Pond_Density + Bufo_bufo + Waterfowl
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Rough_Grass + Quality + Gallinula_chloropus
                + Terrestrial_Other + Scrub_Hedge 
                + Sus_scrofa + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.001)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibian presence has resulted in better model fit


drop1(model2)


## Remove overall terrestrial habitat quality as this variable has the 
## lowest AIC value

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Pond_Density + Bufo_bufo + Waterfowl
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Rough_Grass + Quality + Gallinula_chloropus
                + Terrestrial_Other + Scrub_Hedge + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish presence has resulted in better model fit


drop1(model3)


## Remove water quality as this variable has the lowest AIC value

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Pond_Density + Bufo_bufo + Waterfowl
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Rough_Grass + Gallinula_chloropus
                + Terrestrial_Other + Scrub_Hedge + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of water quality has resulted in better model fit


drop1(model4)


## Remove rough grass as this variable has the lowest AIC value

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Pond_Density + Bufo_bufo + Waterfowl
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Gallinula_chloropus + Terrestrial_Other 
                + Scrub_Hedge + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of rough grass has resulted in better model fit


drop1(model5)


## Remove scrub/hedge as this variable has the lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Pond_Density + Bufo_bufo + Waterfowl
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Gallinula_chloropus + Terrestrial_Other 
                + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit


drop1(model6)


## Remove waterfowl presence as this variable has the lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Pond_Density + Bufo_bufo
                + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Gallinula_chloropus + Terrestrial_Other 
                + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl presence has resulted in better model fit


drop1(model7)


## Remove pond density as this variable has the lowest AIC value

model8 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Bufo_bufo + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Gallinula_chloropus + Terrestrial_Other 
                + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond density has resulted in better model fit


drop1(model8)


## Remove moorhen presence as this variable has the lowest AIC value

model9 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Bufo_bufo + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Overhang + Pond_Area + Ruderals + Woodland
                + Terrestrial_Other + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model9,model8)

## p-value from LRT not significant but minimal reduction in AIC value
## (<2). It is better to choose model with lower AIC value and less
## parameters to avoid overparameterisation.
## removal of moorhen presence has resulted in better model fit


drop1(model9)


## Remove overhang as this variable has the lowest AIC value

model10 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Bufo_bufo + Inflow + Sciurus_carolinensis + Max_Depth
                + Gasterosteus_aculeatus + Permanance + Macrophytes
                + Pond_Area + Ruderals + Woodland
                + Terrestrial_Other + Sus_scrofa,
                family = binomial,
                data = GCNdist_0.001)

anova(model10,model9)

## p-value from LRT not significant but minimal reduction in AIC value
## (<2). It is better to choose model with lower AIC value and less
## parameters to avoid overparameterisation.
## removal of overhang has resulted in better model fit


drop1(model10)


## Remove permanence as this variable has the lowest AIC value

model11 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Lissotriton_vulgaris + Sp_richness 
                 + Bufo_bufo + Inflow + Sciurus_carolinensis + Max_Depth
                 + Gasterosteus_aculeatus + Macrophytes
                 + Pond_Area + Ruderals + Woodland
                 + Terrestrial_Other + Sus_scrofa,
                 family = binomial,
                 data = GCNdist_0.001)

anova(model11,model10)

## p-value from LRT not significant but minimal reduction in AIC value
## (<2). It is better to choose model with lower AIC value and less
## parameters to avoid overparameterisation.
## removal of permanence has resulted in better model fit


drop1(model11)


## Remove smooth newt as this variable has the lowest AIC value

model12 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo + Inflow 
                 + Sciurus_carolinensis + Max_Depth
                 + Gasterosteus_aculeatus + Macrophytes
                 + Pond_Area + Ruderals + Woodland
                 + Terrestrial_Other + Sus_scrofa,
                 family = binomial,
                 data = GCNdist_0.001)

anova(model12,model11)

## p-value from LRT not significant but minimal increase in AIC value
## (<2). It is better to choose model with less parameters to avoid 
## overparameterisation.
## removal of smooth newt has resulted in better model fit


drop1(model12)


## Remove pig presence as this variable has the lowest AIC value

model13 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo + Inflow 
                 + Sciurus_carolinensis + Max_Depth
                 + Gasterosteus_aculeatus + Macrophytes
                 + Pond_Area + Ruderals + Woodland
                 + Terrestrial_Other,
                 family = binomial,
                 data = GCNdist_0.001)

anova(model13,model12)

## p-value from LRT not significant but minimal increase in AIC value
## (<2). It is better to choose model with less parameters to avoid 
## overparameterisation.
## removal of pig presence has resulted in better model fit


drop1(model13)


## Remove macrophytes as this variable has the lowest AIC value

model14 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo + Inflow 
                 + Sciurus_carolinensis + Max_Depth
                 + Gasterosteus_aculeatus
                 + Pond_Area + Ruderals + Woodland
                 + Terrestrial_Other,
                 family = binomial,
                 data = GCNdist_0.001)

anova(model14,model13)

## p-value from LRT not significant and increase in AIC value
## (>2). Removal of macrophytes has reduced model fit.



## Final model: 

GCNm_0.001 <- glmer(Triturus_cristatus ~ (1|sample)
                    + Sp_richness + Bufo_bufo + Inflow 
                    + Sciurus_carolinensis + Max_Depth
                    + Gasterosteus_aculeatus + Macrophytes
                    + Pond_Area + Ruderals + Woodland
                    + Terrestrial_Other,
                    family = binomial,
                    data = GCNdist_0.001)

summary(GCNm_0.001)
anova(GCNm_0.001)
drop1(GCNm_0.001, test = "Chi")   # for binomial/integer y models, 
                                  # statistics for significance of each term

display(GCNm_0.001)
se.ranef(GCNm_0.001)              # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(GCNm_0.001)

## Overdispersion test (chi-square)

1-pchisq(421.517, df=488)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(GCNm_0.001)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(GCNm_0.001, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(GCNm_0.001)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(GCNm_0.001)


## Plot the residuals against each x vairable
plot(sresid ~ GCNdist_0.001$Sp_richness)
plot(sresid ~ GCNdist_0.001$Bufo_bufo)
plot(sresid ~ GCNdist_0.001$Inflow)
plot(sresid ~ GCNdist_0.001$Sciurus_carolinensis)
plot(sresid ~ GCNdist_0.001$Max_Depth)
plot(sresid ~ GCNdist_0.001$Gasterosteus_aculeatus)
plot(sresid ~ GCNdist_0.001$Macrophytes)
plot(sresid ~ GCNdist_0.001$Pond_Area)
plot(sresid ~ GCNdist_0.001$Ruderals)
plot(sresid ~ GCNdist_0.001$Woodland)
plot(sresid ~ GCNdist_0.001$Terrestrial_Other)


## Check model residuals
chkres(GCNm_0.001)

## Plot the fitted data against the observed data
plot(GCNdist_0.001$Triturus_cristatus ~ fitted(GCNm_0.001))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(GCNm_0.001)


## Hosmer and Lemeshow Goodness of Fit Test

hoslem.test(GCNdist_0.001$Triturus_cristatus, fitted(GCNm_0.001))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(GCNm_0.001, newdata=GCNdist_0.001, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(GCNdist_0.001, fit)



## PLOT 15A: Species richness and GCN presence/absence

## Create a range of species richness values which increase by 0.02185
## to predict GCN occupancy as species richness increases
range <- seq(from=min(GCNdist_0.001$Sp_richness), to=max(GCNdist_0.001$Sp_richness), by=0.02185)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=range,
                 Bufo_bufo=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Max_Depth=rep(1.344147, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Macrophytes=rep(25.87302, length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Ruderals=rep('Important', length(range)),
                 Woodland=rep('Important', length(range)),
                 Terrestrial_Other=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.001, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.spp <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p15a <- ggplot() + ggtitle('(a)')
p15a <- p15a + geom_ribbon(aes(x=dat.spp$Sp_richness, ymin = dat.spp$cil, ymax = dat.spp$ciu), alpha = 0.25)
p15a <- p15a + geom_line(aes(x=dat.spp$Sp_richness, y=dat.spp$fit), size = 1)
p15a <- p15a + geom_jitter(aes(x=GCNdist_0.001$Sp_richness, y=GCNdist_0.001$Triturus_cristatus, colour=GCNdist_0.001$Triturus_cristatus), cex=1, height=0.7)
p15a <- p15a + scale_colour_gradient(low="grey48", high="orange")
p15a <- p15a + labs(y = "Probability of GCN presence", x = "Species richness")
p15a <- p15a + coord_cartesian(xlim=c(0,12))
p15a <- p15a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15a



## PLOT 15B: common toad presence and GCN presence/absence 

p15b <- ggplot(box.dat) + ggtitle('(b)')
p15b <- p15b + geom_jitter(aes(x=factor(Bufo_bufo), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15b <- p15b + geom_boxplot(aes(x=factor(Bufo_bufo), y=fit), alpha=0.7, outlier.colour = "blue")
p15b <- p15b + scale_colour_gradient(low="grey48", high="orange")
p15b <- p15b + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p15b <- p15b + labs(y="Probability of GCN presence", x = "Common toad") 
p15b <- p15b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15b



## PLOT 15C: inflow and GCN presence/absence

p15c <- ggplot(box.dat) + ggtitle('(c)')
p15c <- p15c + geom_jitter(aes(x=Inflow, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15c <- p15c + geom_boxplot(aes(x=Inflow, y=fit), alpha=0.7, outlier.colour = "blue")
p15c <- p15c + scale_colour_gradient(low="grey48", high="orange")
p15c <- p15c + labs(y="Probability of GCN presence", x = "Pond inflow") 
p15c <- p15c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15c



## PLOT 15D: grey squirrel presence and GCN presence/absence 

p15d <- ggplot(box.dat) + ggtitle('(d)')
p15d <- p15d + geom_jitter(aes(x=factor(Sciurus_carolinensis), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15d <- p15d + geom_boxplot(aes(x=factor(Sciurus_carolinensis), y=fit), alpha=0.7, outlier.colour = "blue")
p15d <- p15d + scale_colour_gradient(low="grey48", high="orange")
p15d <- p15d + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p15d <- p15d + labs(y="Probability of GCN presence", x = "Grey squirrel") 
p15d <- p15d + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15d



## PLOT 15E: Max depth and GCN presence/absence

## Create a range of depth values which increase by 0.01587302
## to predict GCN occupancy as depth increases
range <- seq(from=min(GCNdist_0.001$Max_Depth), to=max(GCNdist_0.001$Max_Depth), by=0.01587302)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=rep(3.418651, length(range)),
                 Bufo_bufo=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Max_Depth=range,
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Macrophytes=rep(25.87302, length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Ruderals=rep('Important', length(range)),
                 Woodland=rep('Important', length(range)),
                 Terrestrial_Other=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.001, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.MD <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p15e <- ggplot() + ggtitle('(e)')
p15e <- p15e + geom_ribbon(aes(x=dat.MD$Max_Depth, ymin = dat.MD$cil, ymax = dat.MD$ciu), alpha = 0.25)
p15e <- p15e + geom_line(aes(x=dat.MD$Max_Depth, y=dat.MD$fit), size = 1)
p15e <- p15e + geom_jitter(aes(x=GCNdist_0.001$Max_Depth, y=GCNdist_0.001$Triturus_cristatus, colour=GCNdist_0.001$Triturus_cristatus), cex=1, height=0.7)
p15e <- p15e + scale_colour_gradient(low="grey48", high="orange")
p15e <- p15e + labs(y = "Probability of GCN presence", x = "Maximum pond depth (m)")
p15e <- p15e + coord_cartesian(xlim=c(0,10))
p15e <- p15e + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15e



## PLOT 15F: three-spined stickleback presence and GCN presence/absence 

p15f <- ggplot(box.dat) + ggtitle('(f)')
p15f <- p15f + geom_jitter(aes(x=factor(Gasterosteus_aculeatus), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15f <- p15f + geom_boxplot(aes(x=factor(Gasterosteus_aculeatus), y=fit), alpha=0.7, outlier.colour = "blue")
p15f <- p15f + scale_colour_gradient(low="grey48", high="orange")
p15f <- p15f + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p15f <- p15f + labs(y="Probability of GCN presence", x = "Three-spined stickleback") 
p15f <- p15f + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15f



## PLOT 15G: Macrophytes and GCN presence/absence

## Create a range of macrophyte cover values which increase by 0.1984127)
## to predict GCN occupancy as macrophyte cover increases
range <- seq(from=min(GCNdist_0.001$Macrophytes), to=max(GCNdist_0.001$Macrophytes), by=0.1984127)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=rep(3.418651, length(range)),
                 Bufo_bufo=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Max_Depth=rep(1.344147, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Macrophytes=range,
                 Pond_Area=rep(599.9762, length(range)),
                 Ruderals=rep('Important', length(range)),
                 Woodland=rep('Important', length(range)),
                 Terrestrial_Other=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.001, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.MC <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p15g <- ggplot() + ggtitle('(g)')
p15g <- p15g + geom_ribbon(aes(x=dat.MC$Macrophytes, ymin = dat.MC$cil, ymax = dat.MC$ciu), alpha = 0.25)
p15g <- p15g + geom_line(aes(x=dat.MC$Macrophytes, y=dat.MC$fit), size = 1)
p15g <- p15g + geom_jitter(aes(x=GCNdist_0.001$Macrophytes, y=GCNdist_0.001$Triturus_cristatus, colour=GCNdist_0.001$Triturus_cristatus), cex=1, height=0.7)
p15g <- p15g + scale_colour_gradient(low="grey48", high="orange")
p15g <- p15g + labs(y = "Probability of GCN presence", x = "Macrophyte cover (%)")
p15g <- p15g + coord_cartesian(xlim=c(0,100))
p15g <- p15g + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15g



## PLOT 15H: Pond area and GCN presence/absence

## Create a range of pond area values which increase by 18.6012
## to predict GCN occupancy as pond area increases
range <- seq(from=min(GCNdist_0.001$Pond_Area), to=max(GCNdist_0.001$Pond_Area), by=18.6012)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=rep(3.418651, length(range)),
                 Bufo_bufo=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Max_Depth=rep(1.344147, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Macrophytes=rep(25.87302, length(range)),
                 Pond_Area=range,
                 Ruderals=rep('Important', length(range)),
                 Woodland=rep('Important', length(range)),
                 Terrestrial_Other=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.001, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.PA <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p15h <- ggplot() + ggtitle('(g)')
p15h <- p15h + geom_ribbon(aes(x=dat.PA$Pond_Area, ymin = dat.PA$cil, ymax = dat.PA$ciu), alpha = 0.25)
p15h <- p15h + geom_line(aes(x=dat.PA$Pond_Area, y=dat.PA$fit), size = 1)
p15h <- p15h + geom_jitter(aes(x=GCNdist_0.001$Pond_Area, y=GCNdist_0.001$Triturus_cristatus, colour=GCNdist_0.001$Triturus_cristatus), cex=1, height=0.7)
p15h <- p15h + scale_colour_gradient(low="grey48", high="orange")
p15h <- p15h + labs(y = "Probability of GCN presence", x = "Pond area ("~m^2~")")
p15h <- p15h + coord_cartesian(xlim=c(0,10000))
p15h <- p15h + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15h



## PLOT 15I: ruderals and GCN presence/absence

f.Ruderals <- ordered(GCNdist_0.001$Ruderals, levels = c("None", "Some","Important"))
f.Ruderals

p15i <- ggplot(box.dat) + ggtitle('(i)')
p15i <- p15i + geom_jitter(aes(x=f.Ruderals, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15i <- p15i + geom_boxplot(aes(x=f.Ruderals, y=fit), alpha=0.7, outlier.colour = "blue")
p15i <- p15i + scale_colour_gradient(low="grey48", high="orange")
p15i <- p15i + labs(y="Probability of GCN presence", x = "Ruderals") 
p15i <- p15i + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15i



## PLOT 15J: woodland and GCN presence/absence

f.Woodland<- ordered(GCNdist_0.001$Woodland, levels = c("None", "Some","Important"))
f.Woodland

p15j <- ggplot(box.dat) + ggtitle('(i)')
p15j <- p15j + geom_jitter(aes(x=f.Woodland, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15j <- p15j + geom_boxplot(aes(x=f.Woodland, y=fit), alpha=0.7, outlier.colour = "blue")
p15j <- p15j + scale_colour_gradient(low="grey48", high="orange")
p15j <- p15j + labs(y="Probability of GCN presence", x = "Woodland") 
p15j <- p15j + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15j



## PLOT 15K: terrestrial other and GCN presence/absence

f.Terrestrial_Other <- ordered(GCNdist_0.001$Terrestrial_Other, levels = c("None", "Some","Important"))
f.Terrestrial_Other

p15k <- ggplot(box.dat) + ggtitle('(k)')
p15k <- p15k + geom_jitter(aes(x=f.Terrestrial_Other, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15k <- p15k + geom_boxplot(aes(x=f.Terrestrial_Other, y=fit), alpha=0.7, outlier.colour = "blue")
p15k <- p15k + scale_colour_gradient(low="grey48", high="orange")
p15k <- p15k + labs(y="Probability of GCN presence", x = "Terrestrial other") 
p15k <- p15k + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15k


grid.arrange(arrangeGrob(p15a,p15b,p15c,p15d, ncol=4),
             arrangeGrob(p15e,p15f,p15g,p15h, ncol=4),
             arrangeGrob(p15i,p15j,p15k, ncol=3), nrow=3,
             top=textGrob("Predictors of great crested newt (GCN) presence with 0.1% sequence threshold",
                          gp=gpar(cex=3)))




#########
# 0.5%  #
#########

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris + Sciurus_carolinensis
              + Gasterosteus_aculeatus + Fulica_atra + Sus_scrofa 
              + Gallinula_chloropus + Phasianus_colchicus + Esox_lucius
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


GCN_tree <- rpart(f1, data = habitat_0.005, method = "class",
                  minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - L. vulgaris presence
#' - pond substrate
#' - species richness
#' - amphibian presence
#' - B. bufo presence
#' - Terrestrial other
#' - max depth
#' - pond density
#' - inflow
#' - overall terrestrial habitat quality
#' - waterfowl presence
#' - pond area
#' - fish presence
#' - G. aculeatus presence
#' - permanence
#' - macrophytes
#' - G. chloropus presence
#' - ruderals
#' - rough grass
#' - S.carolinensis presence
#' - water quality
#' - scrub/hedge
#' - woodland
#' - overhang
#' 

## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Optimal tree size is high. The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - pollution
#' - outflow
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist_0.005 <- data.frame(habitat_0.005[,c(1,10,14:16,18,20:23,26:34,
                                             44,57,59,65,88,95,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

GCNdist_0.005 < scale(GCNdist_0.005[,c(2:6,26)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Pond_Substrate + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow + Waterfowl
                + Pond_Area + Fish + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Rough_Grass + Sciurus_carolinensis 
                + Quality + Scrub_Hedge + Woodland + Overhang
                + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.005)


## drop1 function:

drop1(model0)


## Remove pond substrate  as this variable has the lowest AIC value

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow + Waterfowl
                + Pond_Area + Fish + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Rough_Grass + Sciurus_carolinensis 
                + Quality + Scrub_Hedge + Woodland + Overhang
                + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.005)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove water quality as this variable has the lowest AIC value

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow + Waterfowl
                + Pond_Area + Fish + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Rough_Grass + Sciurus_carolinensis 
                + Scrub_Hedge + Woodland + Overhang
                + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.005)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of water quality has resulted in better model fit


drop1(model2)


## Remove waterfowl as this variable has the lowest AIC value

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow
                + Pond_Area + Fish + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Rough_Grass + Sciurus_carolinensis 
                + Scrub_Hedge + Woodland + Overhang
                + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.005)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl has resulted in better model fit


drop1(model3)


## Remove overall terrestrial habitat quality as this variable has the 
## lowest AIC value

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow
                + Pond_Area + Fish + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Rough_Grass + Sciurus_carolinensis 
                + Scrub_Hedge + Woodland + Overhang,
                family = binomial,
                data = GCNdist_0.005)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of overall terrestrial habitat quality has resulted in 
## better model fit


drop1(model4)


## Remove rough grass as this variable has the lowest AIC value

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow
                + Pond_Area + Fish + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Sciurus_carolinensis 
                + Scrub_Hedge + Woodland + Overhang,
                family = binomial,
                data = GCNdist_0.005)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of rough grass has resulted in better model fit


drop1(model5)


## Remove scrub/hedge as this variable has the lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow
                + Pond_Area + Fish + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Sciurus_carolinensis 
                + Woodland + Overhang,
                family = binomial,
                data = GCNdist_0.005)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit


drop1(model6)


## Remove fish as this variable has the lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Pond_Density + Inflow
                + Pond_Area + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Sciurus_carolinensis 
                + Woodland + Overhang,
                family = binomial,
                data = GCNdist_0.005)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish has resulted in better model fit


drop1(model7)


## Remove pond density as this variable has the lowest AIC value

model8 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Inflow
                + Pond_Area + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Sciurus_carolinensis 
                + Woodland + Overhang,
                family = binomial,
                data = GCNdist_0.005)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond density has resulted in better model fit


drop1(model8)


## Remove overhang as this variable has the lowest AIC value

model9 <- glmer(Triturus_cristatus ~ (1|sample)
                + Lissotriton_vulgaris + Sp_richness 
                + Amphibians + Bufo_bufo + Terrestrial_Other
                + Max_Depth + Inflow
                + Pond_Area + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Sciurus_carolinensis 
                + Woodland,
                family = binomial,
                data = GCNdist_0.005)

anova(model9,model8)


## p-value from LRT not significant and reduction in AIC value
## so removal of overhang has resulted in better model fit


drop1(model9)


## Remove smooth newt as this variable has the lowest AIC value

model10 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Amphibians + Bufo_bufo 
                + Terrestrial_Other + Max_Depth + Inflow
                + Pond_Area + Gasterosteus_aculeatus 
                + Permanance + Macrophytes + Gallinula_chloropus
                + Ruderals + Sciurus_carolinensis 
                + Woodland,
                family = binomial,
                data = GCNdist_0.005)

anova(model10,model9)

## p-value from LRT significant but reduction in AIC value
## in this instance, it is better to choose the model with less
## parameters to avoid overparameterisation
## so removal of smooth newt has improved model fit


drop1(model10)


## Remove amphibians as this variable has the lowest AIC value

model11 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo 
                 + Terrestrial_Other + Max_Depth + Inflow
                 + Pond_Area + Gasterosteus_aculeatus 
                 + Permanance + Macrophytes + Gallinula_chloropus
                 + Ruderals + Sciurus_carolinensis 
                 + Woodland,
                 family = binomial,
                 data = GCNdist_0.005)

anova(model11,model10)

## p-value from LRT significant but reduction in AIC value (<2)
## in this instance, it is better to choose the model with less
## parameters to avoid overparameterisation
## so removal of amphibian presence has improved model fit


drop1(model11)


## Remove terrestrial other as this variable has the lowest AIC value

model12 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo 
                 + Max_Depth + Inflow
                 + Pond_Area + Gasterosteus_aculeatus 
                 + Permanance + Macrophytes + Gallinula_chloropus
                 + Ruderals + Sciurus_carolinensis + Woodland,
                 family = binomial,
                 data = GCNdist_0.005)

anova(model12,model11)

## p-value from LRT not significant and minimal reduction in AIC value
## in this instance, it is better to choose the model with less
## parameters and lower AIC value to avoid overparameterisation
## so removal of terrestrial other has improved model fit


drop1(model12)


## Remove ruderals as this variable has the lowest AIC value

model13 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo 
                 + Max_Depth + Inflow
                 + Pond_Area + Gasterosteus_aculeatus 
                 + Permanance + Macrophytes + Gallinula_chloropus
                 + Sciurus_carolinensis + Woodland,
                 family = binomial,
                 data = GCNdist_0.005)

anova(model13,model12)

## p-value from LRT not significant and minimal reduction in AIC value
## in this instance, it is better to choose the model with less
## parameters and lower AIC value to avoid overparameterisation
## so removal of ruderals has improved model fit


drop1(model13)


## Remove moorhen presence as this variable has the lowest AIC value

model14 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo 
                 + Max_Depth + Inflow
                 + Pond_Area + Gasterosteus_aculeatus 
                 + Permanance + Macrophytes
                 + Sciurus_carolinensis + Woodland,
                 family = binomial,
                 data = GCNdist_0.005)

anova(model14,model13)

## p-value from LRT not significant and minimal reduction in AIC value
## in this instance, it is better to choose the model with less
## parameters and lower AIC value to avoid overparameterisation
## so removal of moorhen presence has improved model fit


drop1(model14)


## Remove macrophytes this variable has the lowest AIC value

model15 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo 
                 + Max_Depth + Inflow
                 + Pond_Area + Gasterosteus_aculeatus 
                 + Permanance + Sciurus_carolinensis + Woodland,
                 family = binomial,
                 data = GCNdist_0.005)

anova(model15,model14)

## p-value from LRT not significant despite minimal increase in AIC value (<2)
## in this instance, it is better to choose the model with less
## parameters and lower AIC value to avoid overparameterisation
## so removal of macrophytes has improved model fit


drop1(model15)


## Remove max depth as this variable has the lowest AIC value

model16 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo + Inflow
                 + Pond_Area + Gasterosteus_aculeatus 
                 + Permanance + Sciurus_carolinensis + Woodland,
                 family = binomial,
                 data = GCNdist_0.005)

anova(model16,model15)

## p-value from LRT not significant despite minimal increase in AIC value (<2)
## in this instance, it is better to choose the model with less
## parameters and lower AIC value to avoid overparameterisation
## so removal of max depth has improved model fit


drop1(model16)


## Remove pond area as this variable has the lowest AIC value

model17 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Bufo_bufo + Inflow + Gasterosteus_aculeatus 
                 + Permanance + Sciurus_carolinensis + Woodland,
                 family = binomial,
                 data = GCNdist_0.005)

anova(model17,model16)

## p-value from LRT highly significant and increase in AIC value (>2)
## so removal of pond area has reduced model fit



## Final model: 

GCNm_0.005 <- glmer(Triturus_cristatus ~ (1|sample)
                    + Sp_richness + Bufo_bufo + Inflow
                    + Pond_Area + Gasterosteus_aculeatus 
                    + Permanance + Sciurus_carolinensis + Woodland,
                    family = binomial, 
                    data = GCNdist_0.005)

summary(GCNm_0.005)
anova(GCNm_0.005)
drop1(GCNm_0.005, test = "Chi")   # for binomial/integer y models, 
                                  # statistics for significance of each term

display(GCNm_0.005)
se.ranef(GCNm_0.005)              # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(GCNm_0.005)

## Overdispersion test (chi-square)

1-pchisq(364.931, df=491)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(GCNm_0.005)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(GCNm_0.005, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(GCNm_0.005)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(GCNm_0.005)


## Plot the residuals against each x vairable
plot(sresid ~ GCNdist_0.005$Sp_richness)
plot(sresid ~ GCNdist_0.005$Bufo_bufo)
plot(sresid ~ GCNdist_0.005$Inflow)
plot(sresid ~ GCNdist_0.005$Pond_Area)
plot(sresid ~ GCNdist_0.005$Gasterosteus_aculeatus)
plot(sresid ~ GCNdist_0.005$Permanance)
plot(sresid ~ GCNdist_0.005$Sciurus_carolinensis)
plot(sresid ~ GCNdist_0.005$Woodland)


## Check model residuals
chkres(GCNm_0.005)

## Plot the fitted data against the observed data
plot(GCNdist_0.005$Triturus_cristatus ~ fitted(GCNm_0.005))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(GCNm_0.005)


## Hosmer and Lemeshow Goodness of Fit Test

hoslem.test(GCNdist_0.005$Triturus_cristatus, fitted(GCNm_0.005))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(GCNm_0.005, newdata=GCNdist_0.005, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(GCNdist_0.005, fit)



## PLOT 16A: Species richness and GCN presence/absence

## Create a range of species richness values which increase by 0.01786
## to predict GCN occupancy as species richness increases
range <- seq(from=min(GCNdist_0.005$Sp_richness), to=max(GCNdist_0.005$Sp_richness), by=0.01786)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=range,
                 Bufo_bufo=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Permanance=rep('Never dries', length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Woodland=rep('Important', length(range)))

## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.005, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.spp <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p16a <- ggplot() + ggtitle('(a)')
p16a <- p16a + geom_ribbon(aes(x=dat.spp$Sp_richness, ymin = dat.spp$cil, ymax = dat.spp$ciu), alpha = 0.25)
p16a <- p16a + geom_line(aes(x=dat.spp$Sp_richness, y=dat.spp$fit), size = 1)
p16a <- p16a + geom_jitter(aes(x=GCNdist_0.001$Sp_richness, y=GCNdist_0.001$Triturus_cristatus, colour=GCNdist_0.001$Triturus_cristatus), cex=1, height=0.7)
p16a <- p16a + scale_colour_gradient(low="grey48", high="orange")
p16a <- p16a + labs(y = "Probability of GCN presence", x = "Species richness")
p16a <- p16a + coord_cartesian(xlim=c(0,10))
p16a <- p16a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16a



## PLOT 16B: common toad presence and GCN presence/absence 

p16b <- ggplot(box.dat) + ggtitle('(b)')
p16b <- p16b + geom_jitter(aes(x=factor(Bufo_bufo), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p16b <- p16b + geom_boxplot(aes(x=factor(Bufo_bufo), y=fit), alpha=0.7, outlier.colour = "blue")
p16b <- p16b + scale_colour_gradient(low="grey48", high="orange")
p16b <- p16b + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p16b <- p16b + labs(y="Probability of GCN presence", x = "Common toad") 
p16b <- p16b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16b



## PLOT 15C: inflow and GCN presence/absence

p16c <- ggplot(box.dat) + ggtitle('(c)')
p16c <- p16c + geom_jitter(aes(x=Inflow, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p16c <- p16c + geom_boxplot(aes(x=Inflow, y=fit), alpha=0.7, outlier.colour = "blue")
p16c <- p16c + scale_colour_gradient(low="grey48", high="orange")
p16c <- p16c + labs(y="Probability of GCN presence", x = "Pond inflow") 
p16c <- p16c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16c



## PLOT 16D: Pond area and GCN presence/absence

## Create a range of pond area values which increase by 18.6012
## to predict GCN occupancy as pond area increases
range <- seq(from=min(GCNdist_0.005$Pond_Area), to=max(GCNdist_0.005$Pond_Area), by=18.6012)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=rep(2.912698, length(range)),
                 Bufo_bufo=rep(0, length(range)),
                 Inflow=rep('Present', length(range)),
                 Pond_Area=range,
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Permanance=rep('Never dries', length(range)),
                 Sciurus_carolinensis=rep(0, length(range)),
                 Woodland=rep('Important', length(range)))

## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.005, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.PA <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between pond area and GCN occupancy
## with predicted values from model

p16d <- ggplot() + ggtitle('(d)')
p16d <- p16d + geom_ribbon(aes(x=dat.PA$Pond_Area, ymin = dat.PA$cil, ymax = dat.PA$ciu), alpha = 0.25)
p16d <- p16d + geom_line(aes(x=dat.PA$Pond_Area, y=dat.PA$fit), size = 1)
p16d <- p16d + geom_jitter(aes(x=GCNdist_0.005$Pond_Area, y=GCNdist_0.005$Triturus_cristatus, colour=GCNdist_0.005$Triturus_cristatus), cex=1, height=0.7)
p16d <- p16d + scale_colour_gradient(low="grey48", high="orange")
p16d <- p16d + labs(y = "Probability of GCN presence", x = "Pond area ("~m^2~")")
p16d <- p16d + coord_cartesian(xlim=c(0,10000))
p16d <- p16d + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16d



## PLOT 16E: three-spined stickleback presence and GCN presence/absence 

p16e <- ggplot(box.dat) + ggtitle('(e)')
p16e <- p16e + geom_jitter(aes(x=factor(Gasterosteus_aculeatus), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p16e <- p16e + geom_boxplot(aes(x=factor(Gasterosteus_aculeatus), y=fit), alpha=0.7, outlier.colour = "blue")
p16e <- p16e + scale_colour_gradient(low="grey48", high="orange")
p16e <- p16e + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p16e <- p16e + labs(y="Probability of GCN presence", x = "Three-spined stickleback") 
p16e <- p16e + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16e



## PLOT 16F: permanence and GCN presence/absence

f.Permanence <- ordered(GCNdist_0.005$Permanance, levels = c("Never dries", "Rarely dries","Sometimes dries","Dries annually"))
f.Permanence

p16f <- ggplot(box.dat) + ggtitle('(h)')
p16f <- p16f + geom_jitter(aes(x=f.Permanence, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p16f <- p16f + geom_boxplot(aes(x=f.Permanence, y=fit), alpha=0.7, outlier.colour = "blue")
p16f <- p16f + scale_colour_gradient(low="grey48", high="orange")
p16f <- p16f + labs(y="Probability of GCN presence", x = "Pond permanence") 
p16f <- p16f + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16f



## PLOT 16G: grey squirrel presence and GCN presence/absence 

p16g <- ggplot(box.dat) + ggtitle('(g)')
p16g <- p16g + geom_jitter(aes(x=factor(Sciurus_carolinensis), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p16g <- p16g + geom_boxplot(aes(x=factor(Sciurus_carolinensis), y=fit), alpha=0.7, outlier.colour = "blue")
p16g <- p16g + scale_colour_gradient(low="grey48", high="orange")
p16g <- p16g + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p16g <- p16g + labs(y="Probability of GCN presence", x = "Grey squirrel") 
p16g <- p16g + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16g



## PLOT 16H: woodland and GCN presence/absence

f.Woodland<- ordered(GCNdist_0.005$Woodland, levels = c("None", "Some","Important"))
f.Woodland

p16h <- ggplot(box.dat) + ggtitle('(h)')
p16h <- p16h + geom_jitter(aes(x=f.Woodland, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p16h <- p16h + geom_boxplot(aes(x=f.Woodland, y=fit), alpha=0.7, outlier.colour = "blue")
p16h <- p16h + scale_colour_gradient(low="grey48", high="orange")
p16h <- p16h + labs(y="Probability of GCN presence", x = "Woodland") 
p16h <- p16h + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p16h



## PLOT 15K: terrestrial other and GCN presence/absence

f.Terrestrial_Other <- ordered(GCNdist_0.001$Terrestrial_Other, levels = c("None", "Some","Important"))
f.Terrestrial_Other

p15k <- ggplot(box.dat) + ggtitle('(k)')
p15k <- p15k + geom_jitter(aes(x=f.Terrestrial_Other, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p15k <- p15k + geom_boxplot(aes(x=f.Terrestrial_Other, y=fit), alpha=0.7, outlier.colour = "blue")
p15k <- p15k + scale_colour_gradient(low="grey48", high="orange")
p15k <- p15k + labs(y="Probability of GCN presence", x = "Terrestrial other") 
p15k <- p15k + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p15k


grid.arrange(arrangeGrob(p16a,p16b,p16c, ncol=3),
             arrangeGrob(p16d,p16e,p16f, ncol=3),
             arrangeGrob(p16g,p16h, ncol=2),nrow=3,
             top=textGrob("Predictors of great crested newt (GCN) presence with 0.5% sequence threshold",
                          gp=gpar(cex=3)))



#######
# 1%  #
#######

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris + Sciurus_carolinensis
              + Gasterosteus_aculeatus + Sus_scrofa 
              + Gallinula_chloropus + Phasianus_colchicus + Esox_lucius
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


GCN_tree <- rpart(f1, data = habitat_0.01, method = "class",
                  minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - species richness
#' - pond substrate
#' - L. vulgaris presence
#' - pond density
#' - overhang
#' - fish presence
#' - terrestrial other
#' - G. aculeatus presence
#' - permanance
#' - pond area
#' - inflow
#' - scrub/hedge
#' - overall terrestrial habitat quality
#' - max depth
#' - woodland
#' - macrophytes
#' - ruderals
#' - quality
#' - rough grass
#' - G. chloropus presence
#'


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 11 is optimal i.e. 11 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - amphibian presence
#' - waterfowl presence
#' - pollution
#' - outflow
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist_0.01 <- data.frame(habitat_0.01[,c(1,10,14:16,18,20:23,28:34,
                                           57,59,65,95,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

GCNdist_0.01 < scale(GCNdist_0.01[,c(2:6,22)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Substrate + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish + Terrestrial_Other
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Scrub_Hedge + Max_Depth + Woodland
                + Macrophytes + Ruderals + Quality + Rough_Grass
                + Gallinula_chloropus + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.01)

## drop1 function:

drop1(model0)


## Remove pond substrate  as this variable has the lowest AIC value

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish + Terrestrial_Other
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Scrub_Hedge + Max_Depth + Woodland
                + Macrophytes + Ruderals + Quality + Rough_Grass
                + Gallinula_chloropus + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.01)

anova(model1,model0)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of pond substrate has improved model fit


drop1(model1)


## Remove water quality as this variable has the lowest AIC value

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish + Terrestrial_Other
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Scrub_Hedge + Max_Depth + Woodland
                + Macrophytes + Ruderals + Rough_Grass
                + Gallinula_chloropus + Overall_Habitat_Score,
                family = binomial,
                data = GCNdist_0.01)

anova(model2,model1)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of water quality has improved model fit


drop1(model2)


## Remove overall terrestrial habitat quality as this variable has the lowest 
## AIC value

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish + Terrestrial_Other
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Scrub_Hedge + Max_Depth + Woodland
                + Macrophytes + Ruderals + Rough_Grass + Gallinula_chloropus,
                family = binomial,
                data = GCNdist_0.01)

anova(model3,model2)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of  overall terrestrial habitat quality has 
## improved model fit


drop1(model3)


## Remove ruderals as this variable has the lowest AIC value

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish + Terrestrial_Other
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Scrub_Hedge + Max_Depth + Woodland
                + Macrophytes + Rough_Grass + Gallinula_chloropus,
                family = binomial,
                data = GCNdist_0.01)

anova(model4,model3)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of ruderals has improved model fit


drop1(model4)


## Remove scrub/hedge as this variable has the lowest AIC value

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish + Terrestrial_Other
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Max_Depth + Woodland
                + Macrophytes + Rough_Grass + Gallinula_chloropus,
                family = binomial,
                data = GCNdist_0.01)

anova(model5,model4)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of scrub/hedge has improved model fit


drop1(model5)


## Remove rough grass as this variable has the lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish + Terrestrial_Other
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Max_Depth + Woodland
                + Macrophytes + Gallinula_chloropus,
                family = binomial,
                data = GCNdist_0.01)

anova(model6,model5)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of rough grass has improved model fit


drop1(model6)


## Remove terrestrial other as this variable has the lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang + Fish
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Max_Depth + Woodland
                + Macrophytes + Gallinula_chloropus,
                family = binomial,
                data = GCNdist_0.01)

anova(model7,model6)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of terrestrial other has improved model fit


drop1(model7)


## Remove fish presence as this variable has the lowest AIC value

model8 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris
                + Pond_Density + Overhang
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Max_Depth + Woodland
                + Macrophytes + Gallinula_chloropus,
                family = binomial,
                data = GCNdist_0.01)

anova(model8,model7)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of fish presence has improved model fit


drop1(model8)


## Remove pond density as this variable has the lowest AIC value

model9 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris + Overhang
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Max_Depth + Woodland
                + Macrophytes + Gallinula_chloropus,
                family = binomial,
                data = GCNdist_0.01)

anova(model9,model8)

## p-value from LRT not significant and AIC value of model reduced
## therefore removal of pond density has improved model fit


drop1(model9)


## Remove moorhen presence as this variable has the lowest AIC value

model10 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris + Overhang
                + Gasterosteus_aculeatus + Permanance + Pond_Area
                + Inflow + Max_Depth + Woodland + Macrophytes,
                family = binomial,
                data = GCNdist_0.01)

anova(model10,model9)

## p-value from LRT not significant yet minimal difference in AIC value 
## (<2). In these instances, it is better to choose model with less
## parameters to avoid overparameterisation.
## Therefore removal of moorhen presence has improved model fit


drop1(model9)


## Remove smooth newt presence as this variable has the lowest AIC value

model11 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Overhang
                 + Gasterosteus_aculeatus + Permanance + Pond_Area
                 + Inflow + Max_Depth + Woodland + Macrophytes,
                 family = binomial,
                 data = GCNdist_0.01)

anova(model11,model10)

## p-value from LRT not significant yet minimal difference in AIC value 
## (<2). In these instances, it is better to choose model with less
## parameters to avoid overparameterisation.
## Therefore removal of smooth newt presence has improved model fit


drop1(model11)


## Remove permanence as this variable has the lowest AIC value

model12 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Overhang + Gasterosteus_aculeatus 
                 + Pond_Area + Inflow + Max_Depth + Woodland 
                 + Macrophytes,
                 family = binomial,
                 data = GCNdist_0.01)

anova(model12,model11)

## p-value from LRT not significant yet minimal difference in AIC value 
## (<2). In these instances, it is better to choose model with less
## parameters to avoid overparameterisation.
## Therefore removal of permanence has improved model fit


drop1(model12)


## Remove macrophytes as this variable has the lowest AIC value

model13 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Overhang + Gasterosteus_aculeatus 
                 + Pond_Area + Inflow + Max_Depth + Woodland,
                 family = binomial,
                 data = GCNdist_0.01)

anova(model13,model12)

## p-value from LRT not significant yet minimal difference in AIC value 
## (<2). In these instances, it is better to choose model with less
## parameters to avoid overparameterisation.
## Therefore removal of macrophytes has improved model fit


drop1(model13)


## Remove woodland as this variable has the lowest AIC value

model14 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Overhang + Gasterosteus_aculeatus 
                 + Pond_Area + Inflow + Max_Depth,
                 family = binomial,
                 data = GCNdist_0.01)

anova(model14,model13)

## p-value from LRT not significant yet minimal increase in AIC value 
## (<2). In these instances, it is better to choose model with less
## parameters to avoid overparameterisation.
## Therefore removal of woodland has improved model fit


drop1(model14)


## Remove max. depth as this variable has the lowest AIC value

model15 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Overhang + Gasterosteus_aculeatus 
                 + Pond_Area + Inflow,
                 family = binomial,
                 data = GCNdist_0.01)

anova(model15,model14)

## p-value from LRT significant and increase in AIC value (>2). 
## Therefore removal of max. depth has reduced model fit


## Final model: 

GCNm_0.01 <- glmer(Triturus_cristatus ~ (1|sample)
                    + Sp_richness + Overhang + Gasterosteus_aculeatus 
                    + Pond_Area + Inflow + Max_Depth,
                    family = binomial,
                    data = GCNdist_0.01)

summary(GCNm_0.01)
anova(GCNm_0.01)
drop1(GCNm_0.01, test = "Chi")   # for binomial/integer y models, 
                                 # statistics for significance of each term

display(GCNm_0.01)
se.ranef(GCNm_0.01)              # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(GCNm_0.01)

## Overdispersion test (chi-square)

1-pchisq(421.025, df=496)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(GCNm_0.01)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(GCNm_0.01, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(GCNm_0.01)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(GCNm_0.01)


## Plot the residuals against each x vairable
plot(sresid ~ GCNdist_0.01$Sp_richness)
plot(sresid ~ GCNdist_0.01$Overhang)
plot(sresid ~ GCNdist_0.01$Gasterosteus_aculeatus)
plot(sresid ~ GCNdist_0.01$Pond_Area)
plot(sresid ~ GCNdist_0.01$Inflow)
plot(sresid ~ GCNdist_0.01$Max_Depth)


## Check model residuals
chkres(GCNm_0.01)

## Plot the fitted data against the observed data
plot(GCNdist_0.01$Triturus_cristatus ~ fitted(GCNm_0.01))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(GCNm_0.01)


## Hosmer and Lemeshow Goodness of Fit Test

hoslem.test(GCNdist_0.01$Triturus_cristatus, fitted(GCNm_0.01))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(GCNm_0.01, newdata=GCNdist_0.01, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(GCNdist_0.01, fit)



## PLOT 17A: Species richness and GCN presence/absence

## Create a range of species richness values which increase by 0.0159
## to predict GCN occupancy as species richness increases
range <- seq(from=min(GCNdist_0.01$Sp_richness), to=max(GCNdist_0.01$Sp_richness), by=0.0159)


# Create new data frame where only species richness changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=range,
                 Overhang=rep(39.36508, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Inflow=rep('Present', length(range)),
                 Max_Depth=rep(1.344147, length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.01, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.spp <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p17a <- ggplot() + ggtitle('(a)')
p17a <- p17a + geom_ribbon(aes(x=dat.spp$Sp_richness, ymin = dat.spp$cil, ymax = dat.spp$ciu), alpha = 0.25)
p17a <- p17a + geom_line(aes(x=dat.spp$Sp_richness, y=dat.spp$fit), size = 1)
p17a <- p17a + geom_jitter(aes(x=GCNdist_0.01$Sp_richness, y=GCNdist_0.01$Triturus_cristatus, colour=GCNdist_0.01$Triturus_cristatus), cex=1, height=0.7)
p17a <- p17a + scale_colour_gradient(low="grey48", high="orange")
p17a <- p17a + labs(y = "Probability of GCN presence", x = "Species richness")
p17a <- p17a + coord_cartesian(xlim=c(0,8))
p17a <- p17a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p17a



## PLOT 17B: Overhang and GCN presence/absence

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict GCN occupancy as terrestrial overhang increases
range <- seq(from=min(GCNdist_0.01$Overhang), to=max(GCNdist_0.01$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=rep(2.490079, length(range)),
                 Overhang=range,
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Inflow=rep('Present', length(range)),
                 Max_Depth=rep(1.344147, length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.01, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and GCN occupancy
## with predicted values from model

p17b <- ggplot() + ggtitle('(b)')
p17b <- p17b + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p17b <- p17b + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p17b <- p17b + geom_jitter(aes(x=GCNdist_0.01$Overhang, y=GCNdist_0.01$Triturus_cristatus, colour=GCNdist_0.01$Triturus_cristatus), cex=1, height=0.7)
p17b <- p17b + scale_colour_gradient(low="grey48", high="orange")
p17b <- p17b + labs(y = "Probability of GCN presence", x = "Terrestrial overhang (%)")
p17b <- p17b + coord_cartesian(xlim=c(0,100))
p17b <- p17b + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p17b



## PLOT 17C: three-spined stickleback presence and GCN presence/absence 

p17c <- ggplot(box.dat) + ggtitle('(c)')
p17c <- p17c + geom_jitter(aes(x=factor(Gasterosteus_aculeatus), y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p17c <- p17c + geom_boxplot(aes(x=factor(Gasterosteus_aculeatus), y=fit), alpha=0.7, outlier.colour = "blue")
p17c <- p17c + scale_colour_gradient(low="grey48", high="orange")
p17c <- p17c + scale_x_discrete(breaks=c("0", "1"),
                                labels=c("Absent", "Present"))
p17c <- p17c + labs(y="Probability of GCN presence", x = "Three-spined stickleback") 
p17c <- p17c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p17c



## PLOT 17D: Pond area and GCN presence/absence

## Create a range of pond area values which increase by 18.6012
## to predict GCN occupancy as pond area increases
range <- seq(from=min(GCNdist_0.01$Pond_Area), to=max(GCNdist_0.01$Pond_Area), by=18.6012)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=rep(2.490079, length(range)),
                 Overhang=rep(39.36508, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Pond_Area=range,
                 Inflow=rep('Present', length(range)),
                 Max_Depth=rep(1.344147, length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.01, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.PA <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between pond area and GCN occupancy
## with predicted values from model

p17d <- ggplot() + ggtitle('(d)')
p17d <- p17d + geom_ribbon(aes(x=dat.PA$Pond_Area, ymin = dat.PA$cil, ymax = dat.PA$ciu), alpha = 0.25)
p17d <- p17d + geom_line(aes(x=dat.PA$Pond_Area, y=dat.PA$fit), size = 1)
p17d <- p17d + geom_jitter(aes(x=GCNdist_0.01$Pond_Area, y=GCNdist_0.01$Triturus_cristatus, colour=GCNdist_0.01$Triturus_cristatus), cex=1, height=0.7)
p17d <- p17d + scale_colour_gradient(low="grey48", high="orange")
p17d <- p17d + labs(y = "Probability of GCN presence", x = "Pond area ("~m^2~")")
p17d <- p17d + coord_cartesian(xlim=c(0,10000))
p17d <- p17d + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p17d



## PLOT 17E: inflow and GCN presence/absence

p17e <- ggplot(box.dat) + ggtitle('(e)')
p17e <- p17e + geom_jitter(aes(x=Inflow, y=Triturus_cristatus, colour=Triturus_cristatus), width=0.5, height=0.7, cex=1)
p17e <- p17e + geom_boxplot(aes(x=Inflow, y=fit), alpha=0.7, outlier.colour = "blue")
p17e <- p17e + scale_colour_gradient(low="grey48", high="orange")
p17e <- p17e + labs(y="Probability of GCN presence", x = "Pond inflow") 
p17e <- p17e + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p17e



## PLOT 17F: Max depth and GCN presence/absence

## Create a range of depth values which increase by 0.01587302
## to predict GCN occupancy as depth increases
range <- seq(from=min(GCNdist_0.01$Max_Depth), to=max(GCNdist_0.01$Max_Depth), by=0.01587302)


# Create new data frame where only depth changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 Sp_richness=rep(2.490079, length(range)),
                 Overhang=rep(39.36508, length(range)),
                 Gasterosteus_aculeatus=rep(0, length(range)),
                 Pond_Area=rep(599.9762, length(range)),
                 Inflow=rep('Present', length(range)),
                 Max_Depth=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(GCNm_0.01, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.MD <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between species richness and GCN occupancy
## with predicted values from model

p17f <- ggplot() + ggtitle('(f)')
p17f <- p17f + geom_ribbon(aes(x=dat.MD$Max_Depth, ymin = dat.MD$cil, ymax = dat.MD$ciu), alpha = 0.25)
p17f <- p17f + geom_line(aes(x=dat.MD$Max_Depth, y=dat.MD$fit), size = 1)
p17f <- p17f + geom_jitter(aes(x=GCNdist_0.01$Max_Depth, y=GCNdist_0.01$Triturus_cristatus, colour=GCNdist_0.01$Triturus_cristatus), cex=1, height=0.7)
p17f <- p17f + scale_colour_gradient(low="grey48", high="orange")
p17f <- p17f + labs(y = "Probability of GCN presence", x = "Maximum pond depth (m)")
p17f <- p17f + coord_cartesian(xlim=c(0,8))
p17f <- p17f + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p17f


grid.arrange(arrangeGrob(p17a,p17b,p17c, ncol=3),
             arrangeGrob(p17d,p17e,p17f, ncol=3), nrow=2,
             top=textGrob("Predictors of great crested newt (GCN) presence with 1% sequence threshold",
                          gp=gpar(cex=3)))




#######
# 5%  #
#######

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Bufo_bufo + Lissotriton_vulgaris + Cyprinus_carpio
              + Gasterosteus_aculeatus + Gallinula_chloropus 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


GCN_tree <- rpart(f1, data = habitat_0.05, method = "class",
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
#' - pond substrate
#' - terrestrial other
#' - inflow
#' - amphibian presence
#' - ruderals
#' - permanence
#' - B. bufo presence
#' - max depth
#' - macrophytes
#' - fish presence
#' - G. aculeatus presence
#' - pond area
#' - pond density
#' - rough grass
#' - overhang
#' - woodland
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 9 is optimal i.e. 9 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - waterfowl presence
#' - water quality
#' - pollution
#' - outflow
#' - scrub/hedge
#' - overall terrestrial habitat quality
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist_0.05 <- data.frame(habitat_0.05[,c(1,10,14:16,18,20,22:23,26,28:30,
                                           32:33,44,59,65,95,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

GCNdist_0.05 < scale(GCNdist_0.05[,c(2:6,20)])


## Apply step-wise down nested model selection using glmmML

model0 <- glmer(Triturus_cristatus ~ (1|sample) 
                + Lissotriton_vulgaris + Sp_richness + Pond_Substrate
                + Terrestrial_Other + Inflow + Amphibians + Ruderals
                + Permanance + Bufo_bufo + Max_Depth + Macrophytes
                + Fish + Gasterosteus_aculeatus + Pond_Area + Pond_Density
                + Rough_Grass + Overhang + Woodland,
                family = binomial,
                data = GCNdist_0.05)

## glmer could not fit to the data



########
# 10%  #
########

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Fulica_atra + Gallinula_chloropus + Lissotriton_vulgaris
              + Sus_scrofa + Bufo_bufo + Gasterosteus_aculeatus
              + Sciurus_carolinensis + Permanance + Quality 
              + Pond_Substrate + Inflow + Outflow + Pollution 
              + Amphibians + Waterfowl + Fish + Woodland + Rough_Grass 
              + Scrub_Hedge + Ruderals + Terrestrial_Other
              + Overall_Habitat_Score)


GCN_tree <- rpart(f1, data = habitat_0.1, method = "class",
                  minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - species richness
#' - pond substrate
#' - outflow
#' - pond area
#' - pond density
#' - fish presence
#' - ruderals
#' - macrophytes
#' - waterfowl presence
#' - scrub/hedge
#' - overhang
#' - woodland
#' - inflow
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 18 is optimal i.e. 18 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of species that have interaction with GCN
#' - amphibian presence
#' - water quality
#' - pollution
#' - terrestrial other
#' - rough grass
#' - max depth
#' - permanence
#' - overall terrestrial habitat quality
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist_0.1 <- data.frame(habitat_0.1[,c(1,14:16,18,22:24,27:29,31:32,95,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

GCNdist_0.1 < scale(GCNdist_0.1[,c(2:5,15)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Substrate + Outflow + Pond_Area
                + Pond_Density + Fish + Ruderals + Macrophytes + Waterfowl
                + Scrub_Hedge + Overhang + Woodland + Inflow,
                family = binomial,
                data = GCNdist_0.1)


## drop1 function:

drop1(model0)


## Remove pond substrate as this variable had lowest AIC value in previous
## analyses

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Area
                + Pond_Density + Fish + Ruderals + Macrophytes + Waterfowl
                + Scrub_Hedge + Overhang + Woodland + Inflow,
                family = binomial,
                data = GCNdist_0.1)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove scrub/hedge as this variable had lowest AIC value in previous
## analyses

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Area
                + Pond_Density + Fish + Ruderals + Macrophytes + Waterfowl
                + Overhang + Woodland + Inflow,
                family = binomial,
                data = GCNdist_0.1)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit


drop1(model2)


## Remove ruderals as this variable had lowest AIC value in previous
## analyses

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Area
                + Pond_Density + Fish + Macrophytes + Waterfowl
                + Overhang + Woodland + Inflow,
                family = binomial,
                data = GCNdist_0.1)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has resulted in better model fit


drop1(model3)


## Remove fish as this variable had lowest AIC value in previous
## analyses

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Area
                + Pond_Density + Macrophytes + Waterfowl
                + Overhang + Woodland + Inflow,
                family = binomial,
                data = GCNdist_0.1)

anova(model4,model3)

## p-value from LRT not significant and slight reduction in AIC value
## better to choose model with lower AIC value and less parameters to avoid
## overparameterisation
## so removal of fish has resulted in better model fit


drop1(model4)


## Remove waterfowl as this variable had lowest AIC value in previous
## analyses

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Area
                + Pond_Density + Macrophytes
                + Overhang + Woodland + Inflow,
                family = binomial,
                data = GCNdist_0.1)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl presence has resulted in better model fit


drop1(model5)


## Remove pond area as this variable has lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Density + Macrophytes
                + Overhang + Woodland + Inflow,
                family = binomial,
                data = GCNdist_0.1)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond area has resulted in better model fit


drop1(model6)


## Remove inflow as this variable has lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Density + Macrophytes
                + Overhang + Woodland,
                family = binomial,
                data = GCNdist_0.1)

anova(model7,model6)

## p-value from LRT not significant yet minimal reduction in AIC value
## better to choose model with less parameters to avoid overparameterisation
## so removal of inflow has resulted in better model fit


drop1(model7)


## Remove woodland as this variable has lowest AIC value

model8 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Outflow + Pond_Density + Macrophytes
                + Overhang,
                family = binomial,
                data = GCNdist_0.1)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit


drop1(model8)


## Remove outflow as this variable has lowest AIC value

model9 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Macrophytes + Overhang,
                family = binomial,
                data = GCNdist_0.1)

anova(model9,model8)

## p-value from LRT not significant and reduction in AIC value
## so removal of outflow has resulted in better model fit


drop1(model9)


## Remove overhang as this variable has lowest AIC value

model10 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Macrophytes,
                family = binomial,
                data = GCNdist_0.1)

anova(model10,model9)

## p-value from LRT not significant and reduction in AIC value
## so removal of overhang has resulted in better model fit


drop1(model10)


## Remove macrophytes as this variable has lowest AIC value

model11 <- glmer(Triturus_cristatus ~ (1|sample)
                 + Sp_richness + Pond_Density,
                 family = binomial,
                 data = GCNdist_0.1)

anova(model11,model10)

## p-value from LRT not significant and reduction in AIC value
## so removal of macrophytes has resulted in better model fit


drop1(model11)



## Remove pond density as this variable has lowest AIC value

model12 <- glmer(Triturus_cristatus ~ (1|sample) + Sp_richness,
                 family = binomial,
                 data = GCNdist_0.1)

anova(model12,model11)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond density has resulted in better model fit


## Compare to null model to see if model with species richness is better

model13 <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = GCNdist_0.1)

anova(model13,model12)

## p-value from LRT borderline significant and increase in AIC value (<2)
## Removal of species richness has not significantly affected model fit



#######
# 30% #
#######

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Triturus_cristatus ~ Sp_richness + Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Fulica_atra + Gallinula_chloropus + Lissotriton_vulgaris
              + Sus_scrofa + Bufo_bufo + Gasterosteus_aculeatus
              + Sciurus_carolinensis + Permanance + Quality 
              + Pond_Substrate + Inflow + Outflow + Pollution 
              + Amphibians + Waterfowl + Fish + Woodland + Rough_Grass 
              + Scrub_Hedge + Ruderals + Terrestrial_Other
              + Overall_Habitat_Score)


GCN_tree <- rpart(f1, data = habitat_0.3, method = "class",
                  minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(GCN_tree,uniform=TRUE, margin=0.1)
text(GCN_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - species richness
#' - pond density
#' - max depth
#' - pond area
#' - ruderals
#' - L. vulgaris
#' - water quality
#' - macrophytes
#' - permanence
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(GCN_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 4 is optimal i.e. 4 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - presence of other species that have interaction with GCN
#' - amphibian presence
#' - fish presence
#' - waterfowl presence
#' - inflow
#' - outflow 
#' - overhang
#' - rough grass
#' - woodland
#' - scrub/hedge
#' - terrestrial other
#' - pond substrate
#' - pollution
#' - overall terrestrial habitat quality
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

GCNdist_0.3 <- data.frame(habitat_0.1[,c(1,10,14:15,18,20:21,32,65,95,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

GCNdist_0.3 < scale(GCNdist_0.3[,c(2:5,11)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Max_Depth + Pond_Area
                + Ruderals + Lissotriton_vulgaris + Quality + Macrophytes 
                + Permanance,
                family = binomial,
                data = GCNdist_0.3)


## drop1 function:

drop1(model0)


## Remove macrophytes from model as this has lowest AIC value

model1 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Max_Depth + Pond_Area
                + Ruderals + Lissotriton_vulgaris + Quality + Permanance,
                family = binomial,
                data = GCNdist_0.3)

anova(model1,model0)

## p-value from LRT not significant and reduction in model AIC values
## therefore removal of macrophytes has improved model fit


drop1(model1)     # ERROR


## Remove permanence from model as this relationship has been inconsistant
## in previous analysis

model2 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Max_Depth + Pond_Area
                + Ruderals + Lissotriton_vulgaris + Quality,
                family = binomial,
                data = GCNdist_0.3)

anova(model2,model1)

## p-value from LRT not significant and reduction in model AIC values
## therefore removal of permanence has improved model fit


drop1(model2)


## Remove ruderals from model as this relationship has been inconsistant
## in previous analysis

model3 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Max_Depth + Pond_Area
                + Lissotriton_vulgaris + Quality,
                family = binomial,
                data = GCNdist_0.3)

anova(model3,model2)

## p-value from LRT not significant and reduction in model AIC values
## therefore removal of ruderals has improved model fit


drop1(model3)


## Remove water quality from model as this has lowest AIC value

model4 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Max_Depth + Pond_Area
                + Lissotriton_vulgaris,
                family = binomial,
                data = GCNdist_0.3)

anova(model4,model3)

## p-value from LRT not significant and reduction in model AIC values
## therefore removal of water quality has improved model fit


drop1(model4)


## Remove pond area from model as this has lowest AIC value

model5 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Pond_Density + Max_Depth + Lissotriton_vulgaris,
                family = binomial,
                data = GCNdist_0.3)

anova(model5,model4)

## p-value from LRT not significant and reduction in model AIC values
## therefore removal of pond area has improved model fit


drop1(model5)


## Remove pond density from model as this has lowest AIC value

model6 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Max_Depth + Lissotriton_vulgaris,
                family = binomial,
                data = GCNdist_0.3)

anova(model6,model5)

## p-value from LRT not significant and reduction in model AIC values
## therefore removal of pond density has improved model fit


drop1(model6)



## Remove max depth from model as this has lowest AIC value

model7 <- glmer(Triturus_cristatus ~ (1|sample)
                + Sp_richness + Lissotriton_vulgaris,
                family = binomial,
                data = GCNdist_0.3)

anova(model7,model6)

## p-value from LRT not significant and reduction in model AIC values
## therefore removal of max depth has improved model fit


drop1(model7)


## Remove smooth newt presence from model as this has lowest AIC value

model8 <- glmer(Triturus_cristatus ~ (1|sample) + Sp_richness,
                family = binomial,
                data = GCNdist_0.3)

anova(model8,model7)


## p-value from LRT not significant and reduction in model AIC values
## therefore removal of smooth newt presence has improved model fit


## Compare final model from model selection to null model

model9 <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                family = binomial,
                data = GCNdist_0.3)

anova(model9,model8)


## p-value from LRT not significant despite slight increase in model AIC 
## value therefore removal of species richness has not significantly 
## affected model fit

## MODEL NOT BETTER THAN NULL MODEL THEREFORE NO EXPLANATORY VARIABLES
## ADEQUATELY EXPLAIN GREAT CRESTED NEWT OCCURRENCE




#' ----
#' 
#' Now, test HSI score for correlation with GCN presence in different
#' threshold datasets.
#' 


################
# No threshold #
################

HSI_NT <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
                family = binomial,
                data = habitat_NT)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_NT)

anova(HSI_NT, HSInull)    # HSI_NT better


## Model output:
summary(HSI_NT)
anova(HSI_NT)
drop1(HSI_NT, test="Chi")   # for binomial/integer y models, 
                            # statistics for significance of each term

display(HSI_NT)
se.ranef(HSI_NT)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_NT)

## Overdispersion test (chi-square)

1-pchisq(630.552, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_NT)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_NT, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_NT)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_NT)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_NT$HSI)


## Check model residuals
chkres(HSI_NT)

## Plot the fitted data against the observed data
plot(habitat_NT$Triturus_cristatus ~ fitted(HSI_NT))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_NT)



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

hoslem.test(habitat_NT$Triturus_cristatus, fitted(HSI_NT))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 18A: HSI score and GCN presence/absence (no threshold)

## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as pond area increases
range <- seq(from=min(habitat_NT$HSI), to=max(habitat_NT$HSI), by=0.00142)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d2 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


# Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_NT, newdata=d2, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d2, fit, ciu, cil) # This is now the data frame


# Plot showing relationship between pond area and GCN presence/absence 
# with predicted values from model

p18a <- ggplot() + ggtitle('(a) No threshold')
p18a <- p18a + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p18a <- p18a + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p18a <- p18a + geom_jitter(aes(x=habitat_NT$HSI, y=habitat_NT$Triturus_cristatus, colour=habitat_NT$Triturus_cristatus), height=0.7, cex=1)
p18a <- p18a + scale_colour_gradient(low="grey48", high="orange")
p18a <- p18a + xlim(0,1)
p18a <- p18a + labs(y = "Probability of GCN presence", x = "Habitat Suitability Index (HSI) score")
p18a <- p18a + theme(panel.background = element_rect(fill = 'grey98'), 
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position="none")
p18a



#########
# 0.05% #
#########

HSI_0.05 <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
                family = binomial,
                data = habitat_0.0005)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_0.0005)

anova(HSI_0.05, HSInull)    # HSI_0.05 better


## Model output:
summary(HSI_0.05)
anova(HSI_0.05)
drop1(HSI_0.05, test="Chi")   # for binomial/integer y models, 
                              # statistics for significance of each term

display(HSI_0.05)
se.ranef(HSI_0.05)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_0.05)

## Overdispersion test (chi-square)

1-pchisq(579.455, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_0.05)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_0.05, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_0.05)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_0.05)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.0005$HSI)


## Check model residuals
chkres(HSI_0.05)

## Plot the fitted data against the observed data
plot(habitat_0.0005$Triturus_cristatus ~ fitted(HSI_0.05))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_0.05)



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

hoslem.test(habitat_0.0005$Triturus_cristatus, fitted(HSI_0.05))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 18B: HSI score and GCN presence/absence (0.05% threshold)

## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as pond area increases
range <- seq(from=min(habitat_0.0005$HSI), to=max(habitat_0.0005$HSI), by=0.00142)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d2 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


# Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_0.05, newdata=d2, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d2, fit, ciu, cil) # This is now the data frame


# Plot showing relationship between pond area and GCN presence/absence 
# with predicted values from model

p18b <- ggplot() + ggtitle('(b) 0.05%')
p18b <- p18b + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p18b <- p18b + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p18b <- p18b + geom_jitter(aes(x=habitat_0.0005$HSI, y=habitat_0.0005$Triturus_cristatus, colour=habitat_0.0005$Triturus_cristatus), height=0.7, cex=1)
p18b <- p18b + scale_colour_gradient(low="grey48", high="orange")
p18b <- p18b + xlim(0,1)
p18b <- p18b + labs(y = "Probability of GCN presence", x = "Habitat Suitability Index (HSI) score")
p18b <- p18b + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p18b



########
# 0.1% #
########

HSI_0.1 <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
                  family = binomial,
                  data = habitat_0.001)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_0.001)

anova(HSI_0.1, HSInull)    # HSI_0.1 better


## Model output:
summary(HSI_0.1)
anova(HSI_0.1)
drop1(HSI_0.1, test="Chi")   # for binomial/integer y models, 
                             # statistics for significance of each term

display(HSI_0.1)
se.ranef(HSI_0.1)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_0.1)

## Overdispersion test (chi-square)

1-pchisq(569.934, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_0.1)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_0.1, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_0.1)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_0.1)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.001$HSI)


## Check model residuals
chkres(HSI_0.1)

## Plot the fitted data against the observed data
plot(habitat_0.001$Triturus_cristatus ~ fitted(HSI_0.1))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_0.1)



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

hoslem.test(habitat_0.001$Triturus_cristatus, fitted(HSI_0.1))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 18C: HSI score and GCN presence/absence (0.1% threshold)

## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as pond area increases
range <- seq(from=min(habitat_0.001$HSI), to=max(habitat_0.001$HSI), by=0.00142)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d2 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


# Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_0.1, newdata=d2, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d2, fit, ciu, cil) # This is now the data frame


# Plot showing relationship between pond area and GCN presence/absence 
# with predicted values from model

p18c <- ggplot() + ggtitle('(c) 0.1%')
p18c <- p18c + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p18c <- p18c + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p18c <- p18c + geom_jitter(aes(x=habitat_0.001$HSI, y=habitat_0.001$Triturus_cristatus, colour=habitat_0.001$Triturus_cristatus), height=0.7, cex=1)
p18c <- p18c + scale_colour_gradient(low="grey48", high="orange")
p18c <- p18c + xlim(0,1)
p18c <- p18c + labs(y = "Probability of GCN presence", x = "Habitat Suitability Index (HSI) score")
p18c <- p18c + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p18c



########
# 0.5% #
########

HSI_0.5 <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
                 family = binomial,
                 data = habitat_0.005)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_0.005)

anova(HSI_0.5, HSInull)    # HSI_0.5 better


## Model output:
summary(HSI_0.5)
anova(HSI_0.5)
drop1(HSI_0.5, test="Chi")    # for binomial/integer y models, 
                              # statistics for significance of each term

display(HSI_0.1)
se.ranef(HSI_0.1)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_0.5)

## Overdispersion test (chi-square)

1-pchisq(544.658, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_0.5)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_0.5, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_0.5)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_0.5)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.005$HSI)


## Check model residuals
chkres(HSI_0.5)

## Plot the fitted data against the observed data
plot(habitat_0.005$Triturus_cristatus ~ fitted(HSI_0.5))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_0.5)



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

hoslem.test(habitat_0.005$Triturus_cristatus, fitted(HSI_0.5))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 18D: HSI score and GCN presence/absence (0.5% threshold)

## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as pond area increases
range <- seq(from=min(habitat_0.005$HSI), to=max(habitat_0.005$HSI), by=0.00142)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d2 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


# Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_0.5, newdata=d2, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d2, fit, ciu, cil) # This is now the data frame


# Plot showing relationship between pond area and GCN presence/absence 
# with predicted values from model

p18d <- ggplot() + ggtitle('(d) 0.5%')
p18d <- p18d + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p18d <- p18d + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p18d <- p18d + geom_jitter(aes(x=habitat_0.005$HSI, y=habitat_0.005$Triturus_cristatus, colour=habitat_0.005$Triturus_cristatus), height=0.7, cex=1)
p18d <- p18d + scale_colour_gradient(low="grey48", high="orange")
p18d <- p18d + xlim(0,1)
p18d <- p18d + labs(y = "Probability of GCN presence", x = "Habitat Suitability Index (HSI) score")
p18d <- p18d + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p18d



######
# 1% #
######

HSI_1 <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
               family = binomial,
               data = habitat_0.01)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_0.01)

anova(HSI_1, HSInull)    # HSI_1 better


## Model output:
summary(HSI_1)
anova(HSI_1)
drop1(HSI_1, test="Chi")    # for binomial/integer y models, 
                            # statistics for significance of each term

display(HSI_1)
se.ranef(HSI_1)              # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_1)

## Overdispersion test (chi-square)

1-pchisq(516.275, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_1)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_1, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_1)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_1)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.01$HSI)


## Check model residuals
chkres(HSI_1)

## Plot the fitted data against the observed data
plot(habitat_0.005$Triturus_cristatus ~ fitted(HSI_0.5))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_1)



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

hoslem.test(habitat_0.01$Triturus_cristatus, fitted(HSI_1))


## Model does not fit well as p-value is significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 18E: HSI score and GCN presence/absence (1% threshold)

## Create a range of HSI score values which increase by 0.00142
## to predict GCN occupancy as pond area increases
range <- seq(from=min(habitat_0.01$HSI), to=max(habitat_0.01$HSI), by=0.00142)


# Create new data frame where only pond area changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d2 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


# Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_1, newdata=d2, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d2, fit, ciu, cil) # This is now the data frame


# Plot showing relationship between pond area and GCN presence/absence 
# with predicted values from model

p18e <- ggplot() + ggtitle('(e) 1%')
p18e <- p18e + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p18e <- p18e + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p18e <- p18e + geom_jitter(aes(x=habitat_0.01$HSI, y=habitat_0.01$Triturus_cristatus, colour=habitat_0.01$Triturus_cristatus), height=0.7, cex=1)
p18e <- p18e + scale_colour_gradient(low="grey48", high="orange")
p18e <- p18e + xlim(0,1)
p18e <- p18e + labs(y = "Probability of GCN presence", x = "Habitat Suitability Index (HSI) score")
p18e <- p18e + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p18e




######
# 5% #
######

HSI_5 <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
               family = binomial,
               data = habitat_0.05)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_0.05)

anova(HSI_5, HSInull)    # null model better



#######
# 10% #
#######

HSI_10 <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
               family = binomial,
               data = habitat_0.1)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_0.1)

anova(HSI_10, HSInull)    # null model better



#######
# 30% #
#######

HSI_30 <- glmer(Triturus_cristatus ~ (1|sample) + HSI,
                family = binomial,
                data = habitat_0.3)


## Compare to null model

HSInull <- glmer(Triturus_cristatus ~ (1|sample) + 1,
                 family = binomial,
                 data = habitat_0.3)

anova(HSI_30, HSInull)    # null model better



## SHOW ALL PLOTS

grid.arrange(arrangeGrob(p18a,p18b,p18c, ncol=3),
             arrangeGrob(p18d,p18e, ncol=2), nrow=2,
             top=textGrob("Relationship between habitat suitability and great crested newt (GCN) presence at different sequence thresholds",
                          gp=gpar(cex=3)))





#' ---
#'
#' ## Effect of sequence thresholds on environmental predictors of species 
#' richness
#' 
#' Follow a similar approach as predictors of GCN occupancy to identify
#' predictors of pond species richness.
#' 


################
# NO THRESHOLD #
################

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_NT, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - pond substrate
#' - amphibians
#' - overall terrestrial habitat quality
#' - pond area
#' - macrophytes
#' - pond density
#' - permanance
#' - waterfowl
#' - woodland
#' - rough grass
#' - fish
#' - water quality
#' - overhang
#' - max depth
#' - inflow
#' - terrestrial other
#' - permanance
#' - permanence
#' - ruderals
#' - rough grass
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 6 is optimal i.e. 6 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Variables not included in the classification tree were:
#' 
#' - pollution
#' - scrub/hedge
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_NT <- data.frame(habitat_NT[,c(1,10,14:16,18,20:24,26:30,32:34,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_NT < scale(SRdist_NT[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Substrate + Amphibians + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth + Woodland
                + Waterfowl + Quality + Permanance
                + Ruderals + Fish + Terrestrial_Other + Rough_Grass 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

## drop1() function:

drop1(model0)


## Remove amphibians as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Substrate + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth + Woodland
                + Waterfowl + Quality + Permanance
                + Ruderals + Fish + Terrestrial_Other + Rough_Grass 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model1)


## Remove pond substrate as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth + Woodland
                + Waterfowl + Quality + Permanance
                + Ruderals + Fish + Terrestrial_Other + Rough_Grass 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model2)


## Remove fish as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth + Woodland
                + Waterfowl + Quality + Permanance
                + Ruderals + Terrestrial_Other + Rough_Grass 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish has resulted in better model fit


drop1(model3)


## Remove waterfowl as this variable has the lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth + Woodland
                + Quality + Permanance
                + Ruderals + Terrestrial_Other + Rough_Grass 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl has resulted in better model fit


drop1(model4)


## Remove water quality as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth + Woodland
                + Permanance + Ruderals + Terrestrial_Other 
                + Rough_Grass + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of water quality has resulted in better model fit


drop1(model5)


## Remove woodland as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth
                + Permanance + Ruderals + Terrestrial_Other 
                + Rough_Grass + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit


drop1(model6)


## Remove ruderals as this variable has the lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang + Max_Depth
                + Permanance + Terrestrial_Other 
                + Rough_Grass + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has resulted in better model fit


drop1(model7)


## Remove max depth as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang
                + Permanance + Terrestrial_Other 
                + Rough_Grass + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_NT)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of max depth has resulted in better model fit


drop1(model8)


## Remove overall terrestrial habitat as this variable has the 
## lowest AIC value (in conjunction with pond area)

model9 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Pond_Area + Overhang
                + Permanance + Terrestrial_Other 
                + Rough_Grass + Outflow,
                family = poisson,
                data = SRdist_NT)

anova(model9,model8)

## p-value from LRT not significant and reduction in AIC value (<2) so
## better to remove variable to avoid overparameterisation in this
## instance. 
## Removal of overall terrestrial habitat has resulted in better model 
## fit


drop1(model9)


## Remove pond area as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Macrophytes
                + Inflow + Overhang
                + Permanance + Terrestrial_Other 
                + Rough_Grass + Outflow,
                family = poisson,
                data = SRdist_NT)

anova(model10,model9)

## p-value from LRT not significant and reduction in AIC value (<2) so
## better to remove variable to avoid overparameterisation in this
## instance. 
## Removal of pond area has resulted in better model fit


drop1(model10)


## Remove terrestrial other as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Macrophytes
                 + Inflow + Overhang
                 + Permanance + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_NT)

anova(model11,model10)

## p-value from LRT not significant and reduction in AIC value (<2) so
## better to remove variable to avoid overparameterisation in this
## instance. 
## Removal of terrestrial other has resulted in better model fit


drop1(model11)


## Remove inflow as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Macrophytes + Overhang
                 + Permanance + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_NT)

anova(model12,model11)

## p-value from LRT not significant and minimal reduction in AIC value 
## (<2) so better to remove variable to avoid overparameterisation in 
## this instance. 
## Removal of inflow has resulted in better model fit


drop1(model12)


## Remove macrophyte cover as this variable has the lowest AIC value

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang
                 + Permanance + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_NT)

anova(model13,model12)

## p-value from LRT not significant and minimal increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of macrophytes has resulted in better model fit


drop1(model13)


## Remove pond density as this variable has the lowest AIC value

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Permanance + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_NT)

anova(model14,model13)

## p-value from LRT not significant yet slight increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has resulted in better model fit


drop1(model14)


## Remove permanance as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_NT)

anova(model15,model14)

## p-value from LRT significant and increase in AIC value (>2)
## so removal of permanance has reduced model fit



## Final model: 

SRm_NT <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Rough_Grass + Outflow,
                family = poisson,
                data = SRdist_NT)

summary(SRm_NT)
anova(SRm_NT)
drop1(SRm_NT, test = "Chi")   # for binomial/integer y models, 
                              # statistics for significance of each term

display(SRm_NT)
se.ranef(SRm_NT)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_NT)

## Overdispersion test (chi-square)

1-pchisq(461.863, df=498)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_NT)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_NT, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_NT)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_NT)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_NT$Overhang)
plot(sresid ~ SRdist_NT$Scrub_Hedge)
plot(sresid ~ SRdist_NT$Permanance)
plot(sresid ~ SRdist_NT$Rough_Grass)
plot(sresid ~ SRdist_NT$Outflow)


## Check model residuals
chkres(SRm_NT)

## Plot the fitted data against the observed data
plot(SRdist_NT$Sp_richness ~ fitted(SRm_NT))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_NT)


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

hoslem.test(SRdist_NT$Sp_richness, fitted(SRm_NT))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: GCNdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_NT, newdata=SRdist_NT, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_NT, fit)



## PLOT 19A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_NT$Overhang), to=max(SRdist_NT$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range,
                 Scrub_Hedge=rep('None', length(range)),
                 Permanance=rep('Never dries', length(range)),
                 Rough_Grass=rep('Some', length(range)),
                 Outflow=rep('Present', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_NT, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p19a <- ggplot() + ggtitle('(a)')
p19a <- p19a + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p19a <- p19a + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p19a <- p19a + geom_jitter(aes(x=SRdist_NT$Overhang, y=SRdist_NT$Sp_richness, colour=SRdist_NT$Sp_richness), cex=1, height=0.7)
p19a <- p19a + scale_colour_gradient(low="grey48", high="grey60")
p19a <- p19a + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p19a <- p19a + coord_cartesian(xlim=c(0,100))
p19a <- p19a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p19a



## PLOT 19B: rough grass and species richness

f.Rough_Grass <- ordered(SRdist_NT$Rough_Grass, levels = c("None","Some","Important"))
f.Rough_Grass

p19b <- ggplot(box.dat) + ggtitle('(b)')
p19b <- p19b + geom_jitter(aes(x=f.Rough_Grass, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p19b <- p19b + geom_boxplot(aes(x=f.Rough_Grass, y=fit), alpha=0.7, outlier.colour = "red")
p19b <- p19b + scale_colour_gradient(low="grey48", high="grey60")
p19b <- p19b + labs(y="Species richness", x = "Rough grass habitat") 
p19b <- p19b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p19b



## PLOT 19C: outflow and species richness

p19c <- ggplot(box.dat) + ggtitle('(c)')
p19c <- p19c + geom_jitter(aes(x=Outflow, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p19c <- p19c + geom_boxplot(aes(x=Outflow, y=fit), alpha=0.7, outlier.colour = "red")
p19c <- p19c + scale_colour_gradient(low="grey48", high="grey60")
p19c <- p19c + labs(y="Species richness", x = "Pond outflow") 
p19c <- p19c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p19c



grid.arrange(arrangeGrob(p19a,p19b,p19c, ncol=3),
             top=textGrob("Predictors of species richness with no sequence threshold",
                          gp=gpar(cex=3)))




#########
# 0.05% #
#########

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_0.0005, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - pond substrate
#' - woodland
#' - amphibians
#' - waterfowl
#' - pond density
#' - permanence
#' - ruderals
#' - rough grass
#' - max depth
#' - terrestrial other
#' - pond area
#' - scrub/hedge
#' - fish
#' - overall terrestrial habitat
#' - water quality
#' - inflow
#' - macrophytes
#' - outflow
#' - overhang
#' - woodland
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 10 is optimal i.e. 10 explanatory variables. 
#' The best tree will be subjective.
#' 
#' Pollutions was not included in the classification tree.
#' 



#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_0.0005 <- data.frame(habitat_0.0005[,c(1,10,14:16,18,20:24,
                                              26:34,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_0.0005 < scale(SRdist_0.0005[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality + Pond_Substrate 
                + Inflow + Outflow + Amphibians + Waterfowl 
                + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

## drop1() function:

drop1(model0)


## Remove amphibians as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality + Pond_Substrate 
                + Inflow + Outflow + Waterfowl 
                + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model1)


## Remove pond substrate as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality 
                + Inflow + Outflow + Waterfowl 
                + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model2)


## Remove woodland as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                 + Max_Depth + Pond_Area + Pond_Density + Overhang 
                 + Macrophytes + Permanance + Quality 
                 + Inflow + Outflow + Waterfowl 
                 + Fish + Rough_Grass + Scrub_Hedge + Ruderals 
                 + Terrestrial_Other + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.0005)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit


drop1(model3)


## Remove waterfowl as this variable has the lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality 
                + Inflow + Outflow 
                + Fish + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl has resulted in better model fit


drop1(model4)


## Remove ruderals as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality 
                + Inflow + Outflow + Fish + Rough_Grass + Scrub_Hedge
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals  has resulted in better model fit


drop1(model5)


## Remove permanence as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Inflow + Outflow + Fish
                +  Rough_Grass + Scrub_Hedge + Terrestrial_Other 
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value
## so removal of permanence has resulted in better model fit


drop1(model6)


## Remove max depth as this variable has the lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Inflow + Outflow + Fish
                +  Rough_Grass + Scrub_Hedge + Terrestrial_Other 
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of max depth has resulted in better model fit


drop1(model7)


## Remove pond area as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Overhang + Macrophytes + Quality 
                + Inflow + Outflow + Fish +  Rough_Grass + Scrub_Hedge 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model8,model7)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of pond area has improved model fit


drop1(model8)


## Remove fish as this variable has the lowest AIC value

model9 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Overhang + Macrophytes + Quality 
                + Inflow + Outflow +  Rough_Grass + Scrub_Hedge 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model9,model8)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of fish has improved model fit


drop1(model9)


## Remove inflow as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Overhang + Macrophytes + Quality 
                + Outflow +  Rough_Grass + Scrub_Hedge 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.0005)

anova(model10,model9)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of inflow has improved model fit


drop1(model10)


## Remove terrestrial other as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Macrophytes + Quality 
                 + Outflow +  Rough_Grass + Scrub_Hedge 
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.0005)

anova(model11,model10)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of terrestrial other has improved model fit


drop1(model11)


## Remove water quality as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Macrophytes
                 + Outflow +  Rough_Grass + Scrub_Hedge 
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.0005)

anova(model12,model11)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of water quality has improved model fit


drop1(model12)


## Remove macrophytes as this variable has the lowest AIC value

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Outflow 
                 +  Rough_Grass + Scrub_Hedge + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.0005)

anova(model13,model12)

## p-value from LRT not significant although slight increase in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of macrophytes has improved model fit


drop1(model13)


## Remove pond density as this variable has the lowest AIC value

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Rough_Grass + Scrub_Hedge 
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.0005)

anova(model14,model13)

## p-value from LRT not significant although increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has improved model fit


drop1(model14)


## Remove overall terrestrial habitat as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Rough_Grass + Scrub_Hedge,
                 family = poisson,
                 data = SRdist_0.0005)

anova(model15,model14)

## p-value from LRT not significant and minimal increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of overall terrestrial habitat has improved model fit


drop1(model15)


## Remove scrub/hedge as this variable has the lowest AIC value

model16 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.0005)

anova(model16,model15)

## p-value from LRT significant and increase in AIC value (>2)
## so removal of scrub/hedge has reduced model fit




## Final model: 

SRm_0.0005 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Outflow + Rough_Grass + Scrub_Hedge,
                family = poisson,
                data = SRdist_0.0005)

summary(SRm_0.0005)
anova(SRm_0.0005)
drop1(SRm_0.0005, test = "Chi")   # for binomial/integer y models, 
                                  # statistics for significance of each term

display(SRm_0.0005)
se.ranef(SRm_0.0005)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_0.0005)

## Overdispersion test (chi-square)

1-pchisq(482.463, df=496)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_0.0005)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_0.0005, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_0.0005)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_0.0005)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_0.0005$Overhang)
plot(sresid ~ SRdist_0.0005$Outflow)
plot(sresid ~ SRdist_0.0005$Rough_Grass)
plot(sresid ~ SRdist_0.0005$Scrub_Hedge)



## Check model residuals
chkres(SRm_0.0005)

## Plot the fitted data against the observed data
plot(SRdist_0.0005$Sp_richness ~ fitted(SRm_0.0005))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_0.0005)


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

hoslem.test(SRdist_0.0005$Sp_richness, fitted(SRm_0.0005))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: SRdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_0.0005, newdata=SRdist_0.0005, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_0.0005, fit)



## PLOT 20A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_0.0005$Overhang), to=max(SRdist_0.0005$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range,
                 Outflow=rep('Present', length(range)),
                 Rough_Grass=rep('Some', length(range)),
                 Scrub_Hedge=rep('None', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_0.0005, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p20a <- ggplot() + ggtitle('(a)')
p20a <- p20a + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p20a <- p20a + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p20a <- p20a + geom_jitter(aes(x=SRdist_0.0005$Overhang, y=SRdist_0.0005$Sp_richness, colour=SRdist_0.0005$Sp_richness), cex=1, height=0.7)
p20a <- p20a + scale_colour_gradient(low="grey48", high="grey60")
p20a <- p20a + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p20a <- p20a + coord_cartesian(xlim=c(0,100))
p20a <- p20a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p20a



## PLOT 20B: outflow and species richness

p20b <- ggplot(box.dat) + ggtitle('(b)')
p20b <- p20b + geom_jitter(aes(x=Outflow, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p20b <- p20b + geom_boxplot(aes(x=Outflow, y=fit), alpha=0.7, outlier.colour = "red")
p20b <- p20b + scale_colour_gradient(low="grey48", high="grey60")
p20b <- p20b + labs(y="Species richness", x = "Pond outflow") 
p20b <- p20b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p20b



## PLOT 20C: rough grass and species richness

f.Rough_Grass <- ordered(SRdist_0.0005$Rough_Grass, levels = c("None","Some","Important"))
f.Rough_Grass

p20c <- ggplot(box.dat) + ggtitle('(c)')
p20c <- p20c + geom_jitter(aes(x=f.Rough_Grass, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p20c <- p20c + geom_boxplot(aes(x=f.Rough_Grass, y=fit), alpha=0.7, outlier.colour = "red")
p20c <- p20c + scale_colour_gradient(low="grey48", high="grey60")
p20c <- p20c + labs(y="Species richness", x = "Rough grass habitat") 
p20c <- p20c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p20c



## PLOT 20D: scrub/hedge and species richness

f.Scrub_Hedge <- ordered(SRdist_0.0005$Scrub_Hedge, levels = c("None","Some","Important"))
f.Scrub_Hedge

p20d <- ggplot(box.dat) + ggtitle('(d)')
p20d <- p20d + geom_jitter(aes(x=f.Scrub_Hedge, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p20d <- p20d + geom_boxplot(aes(x=f.Scrub_Hedge, y=fit), alpha=0.7, outlier.colour = "red")
p20d <- p20d + scale_colour_gradient(low="grey48", high="grey60")
p20d <- p20d + labs(y="Species richness", x = "Scrub/hedge habitat") 
p20d <- p20d + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p20d



grid.arrange(arrangeGrob(p20a,p20b, ncol=2),
             arrangeGrob(p20c,p20d, ncol=2), nrow=2,
             top=textGrob("Predictors of species richness with 0.05% sequence threshold",
                          gp=gpar(cex=3)))




########
# 0.1% #
########

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_0.001, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - pond substrate
#' - woodland
#' - amphibians
#' - pond density
#' - pollution
#' - rough grass
#' - pond area
#' - permanence
#' - waterfowl
#' - overhang
#' - terrestrial other
#' - scrub/hedge
#' - fish
#' - water quality
#' - macrophytes
#' - overall terrestrial habitat quality
#' - max depth
#' - ruderals
#' - outflow
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 8 is optimal i.e. 8 explanatory variables. 
#' The best tree will be subjective.
#' 
#' All variables excluding inflow were included in the classification tree.
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_0.001 <- data.frame(habitat_0.001[,c(1,10,14:16,18,20:22,24:34,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_0.001 < scale(SRdist_0.001[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality + Pond_Substrate 
                + Outflow + Pollution + Amphibians + Waterfowl 
                + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

## drop1() function:

drop1(model0)


## Remove amphibians as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality + Pond_Substrate 
                + Outflow + Pollution + Waterfowl 
                + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model1)


## Remove pond substrate as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality 
                + Outflow + Pollution + Waterfowl 
                + Fish + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model2)


## Remove woodland as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality 
                + Outflow + Pollution + Waterfowl 
                + Fish + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit


drop1(model3)


## Remove waterfowl as this variable has the lowest AIC value (alongside fish)

model4 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality 
                + Outflow + Pollution + Fish + Rough_Grass + Scrub_Hedge 
                + Ruderals + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl has resulted in better model fit


drop1(model4)


## Remove fish as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality 
                + Outflow + Pollution + Rough_Grass + Scrub_Hedge 
                + Ruderals + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model5,model4)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC/less parameters to avoid 
## overparameterisation
## so removal of fish has resulted in better model fit


drop1(model5)


## Remove permanence as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Outflow + Pollution 
                + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model6,model5)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of permanence has resulted in better model fit


drop1(model6)


## Remove ruderals as this variable has the lowest AIC value (alongside
## max depth)

model7 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Outflow + Pollution 
                + Rough_Grass + Scrub_Hedge
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model7,model6)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of ruderals has resulted in better model fit


drop1(model7)


## Remove max depth as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Area + Pond_Density + Overhang + Macrophytes 
                + Quality + Outflow + Pollution + Rough_Grass + Scrub_Hedge
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model8,model7)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of max depth has improved model fit


drop1(model8)


## Remove pond area as this variable has the lowest AIC value

model9 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Overhang + Macrophytes 
                + Quality + Outflow + Pollution + Rough_Grass + Scrub_Hedge
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model9,model8)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of pond area has improved model fit


drop1(model9)


## Remove pollution as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Overhang + Macrophytes 
                + Quality + Outflow + Rough_Grass + Scrub_Hedge
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.001)

anova(model10,model9)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of pollution has improved model fit


drop1(model10)


## Remove water quality as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Macrophytes 
                 + Outflow + Rough_Grass + Scrub_Hedge
                 + Terrestrial_Other + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.001)

anova(model11,model10)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of water quality has improved model fit


drop1(model11)


## Remove terrestrial other as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Macrophytes + Outflow 
                 + Rough_Grass + Scrub_Hedge + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.001)

anova(model12,model11)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of terrestrial other has improved model fit


drop1(model12)


## Remove pond density as this variable has the lowest AIC value

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Macrophytes + Outflow 
                 + Rough_Grass + Scrub_Hedge + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.001)


anova(model13,model12)

## p-value from LRT not significant although slight increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has improved model fit


drop1(model13)


## Remove macrophytes as this variable has the lowest AIC value

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Rough_Grass + Scrub_Hedge 
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.001)

anova(model14,model13)

## p-value from LRT not significant although increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of macrophytes has improved model fit


drop1(model14)


## Remove overall terrestrial habitat as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Rough_Grass + Scrub_Hedge,
                 family = poisson,
                 data = SRdist_0.001)

anova(model15,model14)

## p-value from LRT not significant although increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of overall terrestrial habitat has improved model fit


drop1(model15)


## Remove rough grass as this variable has the lowest AIC value

model16 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Scrub_Hedge,
                 family = poisson,
                 data = SRdist_0.001)

anova(model16,model15)

## p-value from LRT significant and increase in AIC value (>2)
## so removal of rough grass has reduced model fit


## Final model: 

SRm_0.001 <- glmer(Sp_richness ~ (1|sample)
                    + Overhang + Outflow + Rough_Grass + Scrub_Hedge,
                    family = poisson,
                    data = SRdist_0.001)

summary(SRm_0.001)
anova(SRm_0.001)
drop1(SRm_0.001, test = "Chi")   # for binomial/integer y models, 
                                 # statistics for significance of each term

display(SRm_0.001)
se.ranef(SRm_0.001)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_0.001)

## Overdispersion test (chi-square)

1-pchisq(484.798, df=496)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_0.001)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_0.001, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_0.001)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_0.001)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_0.001$Overhang)
plot(sresid ~ SRdist_0.001$Outflow)
plot(sresid ~ SRdist_0.001$Rough_Grass)
plot(sresid ~ SRdist_0.001$Scrub_Hedge)



## Check model residuals
chkres(SRm_0.001)

## Plot the fitted data against the observed data
plot(SRdist_0.001$Sp_richness ~ fitted(SRm_0.001))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_0.001)


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

hoslem.test(SRdist_0.001$Sp_richness, fitted(SRm_0.001))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: SRdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_0.001, newdata=SRdist_0.001, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_0.001, fit)



## PLOT 21A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_0.001$Overhang), to=max(SRdist_0.001$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range,
                 Outflow=rep('Present', length(range)),
                 Rough_Grass=rep('Some', length(range)),
                 Scrub_Hedge=rep('None', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_0.001, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p21a <- ggplot() + ggtitle('(a)')
p21a <- p21a + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p21a <- p21a + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p21a <- p21a + geom_jitter(aes(x=SRdist_0.001$Overhang, y=SRdist_0.001$Sp_richness, colour=SRdist_0.001$Sp_richness), cex=1, height=0.7)
p21a <- p21a + scale_colour_gradient(low="grey48", high="grey60")
p21a <- p21a + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p21a <- p21a + coord_cartesian(xlim=c(0,100))
p21a <- p21a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p21a



## PLOT 21B: outflow and species richness

p21b <- ggplot(box.dat) + ggtitle('(b)')
p21b <- p21b + geom_jitter(aes(x=Outflow, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p21b <- p21b + geom_boxplot(aes(x=Outflow, y=fit), alpha=0.7, outlier.colour = "red")
p21b <- p21b + scale_colour_gradient(low="grey48", high="grey60")
p21b <- p21b + labs(y="Species richness", x = "Pond outflow") 
p21b <- p21b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p21b



## PLOT 21C: rough grass and species richness

f.Rough_Grass <- ordered(SRdist_0.001$Rough_Grass, levels = c("None","Some","Important"))
f.Rough_Grass

p21c <- ggplot(box.dat) + ggtitle('(c)')
p21c <- p21c + geom_jitter(aes(x=f.Rough_Grass, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p21c <- p21c + geom_boxplot(aes(x=f.Rough_Grass, y=fit), alpha=0.7, outlier.colour = "red")
p21c <- p21c + scale_colour_gradient(low="grey48", high="grey60")
p21c <- p21c + labs(y="Species richness", x = "Rough grass habitat") 
p21c <- p21c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p21c



## PLOT 21D: scrub/hedge and species richness

f.Scrub_Hedge <- ordered(SRdist_0.001$Scrub_Hedge, levels = c("None","Some","Important"))
f.Scrub_Hedge

p21d <- ggplot(box.dat) + ggtitle('(d)')
p21d <- p21d + geom_jitter(aes(x=f.Scrub_Hedge, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p21d <- p21d + geom_boxplot(aes(x=f.Scrub_Hedge, y=fit), alpha=0.7, outlier.colour = "red")
p21d <- p21d + scale_colour_gradient(low="grey48", high="grey60")
p21d <- p21d + labs(y="Species richness", x = "Scrub/hedge habitat") 
p21d <- p21d + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p21d



grid.arrange(arrangeGrob(p21a,p21b, ncol=2),
             arrangeGrob(p21c,p21d, ncol=2), nrow=2,
             top=textGrob("Predictors of species richness with 0.1% sequence threshold",
                          gp=gpar(cex=3)))




########
# 0.5% #
########

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_0.005, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - pond substrate
#' - rough grass
#' - macrophytes
#' - amphibians
#' - max depth
#' - water quality
#' - pond area
#' - permanence
#' - overhang
#' - pond density
#' - terrestrial other
#' - woodland
#' - fish
#' - waterfowl
#' - ruderals
#' - inflow
#' - overall terrestrial habitat
#' - outflow
#' - scrub/hedge
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 4 is optimal i.e. 4 explanatory variables. 
#' The best tree will be subjective.
#' 
#' All variables excluding pollution were included in the classification 
#' tree.
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_0.005 <- data.frame(habitat_0.005[,c(1,10,14:16,18,20:24,
                                            26:34,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_0.005 < scale(SRdist_0.005[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality + Pond_Substrate 
                + Outflow + Amphibians + Fish + Waterfowl + Woodland 
                + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

## drop1() function:

drop1(model0)


## Remove amphibians as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality + Pond_Substrate 
                + Outflow + Fish + Waterfowl + Woodland 
                + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model1)


## Remove pond substrate as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality
                + Outflow + Fish + Waterfowl + Woodland 
                + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model2)


## Remove woodland as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality
                + Outflow + Fish + Waterfowl  
                + Rough_Grass + Scrub_Hedge + Ruderals 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit


drop1(model3)


## Remove ruderals as this variable has the lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Permanance + Quality + Outflow 
                + Fish + Waterfowl + Rough_Grass + Scrub_Hedge 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has resulted in better model fit


drop1(model4)


## Remove permanence as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Outflow 
                + Fish + Waterfowl + Rough_Grass + Scrub_Hedge 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value 
## so removal of permanence has resulted in better model fit


drop1(model5)


## Remove waterfowl as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Outflow 
                + Fish + Rough_Grass + Scrub_Hedge 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model6,model5)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of waterfowl has resulted in better model fit


drop1(model6)


## Remove fish as this variable has the lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Outflow + Rough_Grass 
                + Scrub_Hedge + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish has resulted in better model fit


drop1(model7)


## Remove scrub/hedge as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Area + Pond_Density + Overhang 
                + Macrophytes + Quality + Outflow + Rough_Grass 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has improved model fit


drop1(model8)


## Remove pond area as this variable has the lowest AIC value

model9 <- glmer(Sp_richness ~ (1|sample)
                + Max_Depth + Pond_Density + Overhang 
                + Macrophytes + Quality + Outflow + Rough_Grass 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model9,model8)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of pond area has improved model fit


drop1(model9)


## Remove max depth as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Pond_Density + Overhang + Macrophytes 
                + Quality + Outflow + Rough_Grass 
                + Terrestrial_Other + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.005)

anova(model10,model9)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of max depth has improved model fit


drop1(model10)


## Remove macrophytes as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Quality + Outflow 
                 + Rough_Grass + Terrestrial_Other + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.005)

anova(model11,model10)

## p-value from LRT not significant and slight decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of macrophytes has improved model fit


drop1(model11)


## Remove terrestrial other as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Quality + Outflow 
                 + Rough_Grass + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.005)

anova(model12,model11)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of terrestrial other has improved model fit


drop1(model12)


## Remove water quality as this variable has the lowest AIC value

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Pond_Density + Overhang + Outflow 
                 + Rough_Grass + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.005)

anova(model13,model12)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of water quality has improved model fit


drop1(model13)


## Remove pond density as this variable has the lowest AIC value

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Rough_Grass + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.005)

anova(model14,model13)

## p-value from LRT not significant and slight decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has improved model fit


drop1(model14)


## Remove overall terrestrial habitat as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Outflow + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.005)

anova(model15,model14)

## p-value from LRT significant and increase in AIC value (>2)
## so removal of overall terrestrial habitat has reduced model fit



## Final model: 

SRm_0.005 <- glmer(Sp_richness ~ (1|sample)
                   + Overhang + Outflow + Rough_Grass + Overall_Habitat_Score,
                   family = poisson,
                   data = SRdist_0.005)

summary(SRm_0.005)
anova(SRm_0.005)
drop1(SRm_0.005, test = "Chi")   # for binomial/integer y models, 
                                 # statistics for significance of each term

display(SRm_0.005)
se.ranef(SRm_0.005)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_0.005)

## Overdispersion test (chi-square)

1-pchisq(585.732, df=496)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_0.005)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_0.005, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_0.005)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_0.005)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_0.005$Overhang)
plot(sresid ~ SRdist_0.005$Outflow)


## Check model residuals
chkres(SRm_0.005)

## Plot the fitted data against the observed data
plot(SRdist_0.005$Sp_richness ~ fitted(SRm_0.005))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_0.005)


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

hoslem.test(SRdist_0.005$Sp_richness, fitted(SRm_0.005))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: SRdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_0.005, newdata=SRdist_0.005, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_0.005, fit)



## PLOT 22A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_0.005$Overhang), to=max(SRdist_0.005$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range,
                 Outflow=rep('Present', length(range)),
                 Rough_Grass=rep('Some', length(range)),
                 Overall_Habitat_Score=rep('Moderate', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_0.005, newdata=d3, se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p22a <- ggplot() + ggtitle('(a)')
p22a <- p22a + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p22a <- p22a + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p22a <- p22a + geom_jitter(aes(x=SRdist_0.005$Overhang, y=SRdist_0.005$Sp_richness, colour=SRdist_0.005$Sp_richness), cex=1, height=0.7)
p22a <- p22a + scale_colour_gradient(low="grey48", high="grey60")
p22a <- p22a + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p22a <- p22a + coord_cartesian(xlim=c(0,100))
p22a <- p22a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p22a


## PLOT 22B: outflow and species richness

p22b <- ggplot(box.dat) + ggtitle('(b)')
p22b <- p22b + geom_jitter(aes(x=Outflow, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p22b <- p22b + geom_boxplot(aes(x=Outflow, y=fit), alpha=0.7, outlier.colour = "red")
p22b <- p22b + scale_colour_gradient(low="grey48", high="grey60")
p22b <- p22b + labs(y="Species richness", x = "Pond outflow") 
p22b <- p22b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p22b


## PLOT 22C: rough grass and species richness

f.Rough_Grass <- ordered(SRdist$Rough_Grass, levels = c("None", "Some","Important"))
f.Rough_Grass

p22c <- ggplot(box.dat) + ggtitle('(c)')
p22c <- p22c + geom_jitter(aes(x=f.Rough_Grass, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p22c <- p22c + geom_boxplot(aes(x=f.Rough_Grass, y=fit), alpha=0.7, outlier.colour = "red")
p22c <- p22c + scale_colour_gradient(low="grey48", high="grey60")
p22c <- p22c + labs(y="Species richness", x = "Pond outflow") 
p22c <- p22c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p22c


## PLOT 22D: overall terrestrial habitat and species richness

f.Overall_Habitat_Score <- ordered(SRdist_0.005$Overall_Habitat_Score, 
                                   levels = c("Poor", "Moderate","Good"))
f.Overall_Habitat_Score

p22d <- ggplot(box.dat) + ggtitle('(d)')
p22d <- p22d + geom_jitter(aes(x=f.Overall_Habitat_Score, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p22d <- p22d + geom_boxplot(aes(x=f.Overall_Habitat_Score, y=fit), alpha=0.7, outlier.colour = "red")
p22d <- p22d + scale_colour_gradient(low="grey48", high="grey60")
p22d <- p22d + labs(y="Species richness", x = "Pond outflow") 
p22d <- p22d + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p22d


## ALL PLOTS

grid.arrange(arrangeGrob(p22a,p22b, ncol=2),
             arrangeGrob(p22c,p22d, ncol=2), nrow=2,
             top=textGrob("Predictors of species richness with 0.5% sequence threshold",
                          gp=gpar(cex=3)))




######
# 1% #
######

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_0.01, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - overhang
#' - amphibians
#' - pond substrate
#' - pond density
#' - fish
#' - woodland
#' - macrophytes
#' - pond area
#' - max depth
#' - outflow
#' - water quality
#' - waterfowl
#' - ruderals
#' - terrestrial other
#' - scrub/hedge
#' - permanence
#' - rough grass
#' - inflow
#' - overall terrestrial habitat
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 14 is optimal i.e. 14 explanatory variables. 
#' The best tree will be subjective.
#' 
#' All variables excluding pollution were included in the classification 
#' tree.
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_0.01 <- data.frame(habitat_0.01[,c(1,10,14:16,18,20:24,26:34,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_0.01 < scale(SRdist_0.01[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Amphibians + Pond_Substrate + Pond_Density 
                + Fish + Woodland + Macrophytes + Pond_Area + Max_Depth
                + Outflow + Quality + Waterfowl + Ruderals 
                + Terrestrial_Other + Scrub_Hedge + Permanance 
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)


## drop1() function:

drop1(model0)


## Remove amphibians as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Substrate + Pond_Density 
                + Fish + Woodland + Macrophytes + Pond_Area + Max_Depth
                + Outflow + Quality + Waterfowl + Ruderals 
                + Terrestrial_Other + Scrub_Hedge + Permanance 
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model1)


## Remove pond substrate as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Fish + Woodland 
                + Macrophytes + Pond_Area + Max_Depth
                + Outflow + Quality + Waterfowl + Ruderals 
                + Terrestrial_Other + Scrub_Hedge + Permanance 
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model2)


## Remove permanence as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Fish + Woodland 
                + Macrophytes + Pond_Area + Max_Depth
                + Outflow + Quality + Waterfowl + Ruderals 
                + Terrestrial_Other + Scrub_Hedge
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of permanance has resulted in better model fit


drop1(model3)


## Remove ruderals as this variable has the lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Fish + Woodland 
                + Macrophytes + Pond_Area + Max_Depth
                + Outflow + Quality + Waterfowl 
                + Terrestrial_Other + Scrub_Hedge
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has resulted in better model fit


drop1(model4)


## Remove woodland as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Fish + Macrophytes 
                + Pond_Area + Max_Depth + Outflow + Quality 
                + Waterfowl + Terrestrial_Other + Scrub_Hedge
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value 
## so removal of woodland has resulted in better model fit


drop1(model5)


## Remove fish as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Macrophytes 
                + Pond_Area + Max_Depth + Outflow + Quality 
                + Waterfowl + Terrestrial_Other + Scrub_Hedge
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value 
## so removal of fish has resulted in better model fit


drop1(model6)


## Remove scrub/hedge as this variable has the lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Macrophytes + Pond_Area 
                + Max_Depth + Outflow + Quality + Waterfowl 
                + Terrestrial_Other + Rough_Grass + Inflow 
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit


drop1(model7)


## Remove waterfowl presence as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Macrophytes + Pond_Area 
                + Max_Depth + Outflow + Quality + Terrestrial_Other 
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl presence has improved model fit


drop1(model8)


## Remove water quality as this variable has the lowest AIC value

model9 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Macrophytes + Pond_Area 
                + Max_Depth + Outflow + Terrestrial_Other 
                + Rough_Grass + Inflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model9,model8)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of water quality has improved model fit


drop1(model9)


## Remove max depth as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Density + Macrophytes + Pond_Area 
                + Outflow + Terrestrial_Other + Rough_Grass + Inflow 
                + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.01)

anova(model10,model9)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of max depth has improved model fit


drop1(model10)


## Remove pond area as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Density + Macrophytes + Outflow 
                 + Terrestrial_Other + Rough_Grass + Inflow 
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.01)

anova(model11,model10)

## p-value from LRT not significant and slight decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond area has improved model fit


drop1(model11)


## Remove terrestrial other as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Density + Macrophytes + Outflow 
                 + Rough_Grass + Inflow + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.01)

anova(model12,model11)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of macrophytes has improved model fit


drop1(model12)


## Remove macrophytes as this variable has the lowest AIC value

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Density + Outflow + Rough_Grass 
                 + Inflow + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.01)

anova(model13,model12)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of macrophytes has improved model fit


drop1(model13)


## Remove inflow as this variable has the lowest AIC value (alongside
## pond density)

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Density + Outflow + Rough_Grass 
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.01)

anova(model14,model13)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of inflow has improved model fit


drop1(model14)


## Remove outflow as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Density + Rough_Grass 
                 + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.01)

anova(model15,model14)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of outflow has improved model fit


drop1(model15)


## Remove pond density as this variable has the lowest AIC value

model16 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Rough_Grass + Overall_Habitat_Score,
                 family = poisson,
                 data = SRdist_0.01)

anova(model16,model15)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has improved model fit


drop1(model16)


## Remove overall terrestrial habitat as this variable has the lowest AIC value

model17 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.01)

anova(model17,model16)

## p-value from LRT not significant yet slight increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of overall terrestrial habitat has improved model fit


drop1(model17)


## Remove rough grass as this variable has the lowest AIC value

model18 <- glmer(Sp_richness ~ (1|sample) + Overhang,
                 family = poisson,
                 data = SRdist_0.01)

anova(model18,model17)


## p-value from LRT not significant yet slight increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of overall terrestrial habitat has improved model fit


## Compare null model to model with overhang 

model19 <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = SRdist_0.01)

anova(model19,model18)

## p-value from LRT significant and increase in AIC value (>2)
## model with overhang better fit to data than null model



## Final model: 

SRm_0.01 <- glmer(Sp_richness ~ (1|sample) + Overhang,
                   family = poisson,
                   data = SRdist_0.01)

summary(SRm_0.01)
anova(SRm_0.01)
drop1(SRm_0.01, test = "Chi")   # for binomial/integer y models, 
                                # statistics for significance of each term

display(SRm_0.01)
se.ranef(SRm_0.01)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_0.01)

## Overdispersion test (chi-square)

1-pchisq(543.303, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_0.01)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_0.01, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_0.01)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_0.01)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_0.01$Overhang)


## Check model residuals
chkres(SRm_0.01)

## Plot the fitted data against the observed data
plot(SRdist_0.01$Sp_richness ~ fitted(SRm_0.01))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_0.01)


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

hoslem.test(SRdist_0.01$Sp_richness, fitted(SRm_0.01))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: SRdist_NT
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_0.01, newdata=SRdist_0.01, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_0.01, fit)



## PLOT 23A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_0.01$Overhang), to=max(SRdist_0.01$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_0.01, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p23a <- ggplot() + ggtitle('Predictors of species richness with 1% sequence threshold')
p23a <- p23a + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p23a <- p23a + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p23a <- p23a + geom_jitter(aes(x=SRdist_0.01$Overhang, y=SRdist_0.01$Sp_richness, colour=SRdist_0.01$Sp_richness), cex=1, height=0.7)
p23a <- p23a + scale_colour_gradient(low="grey48", high="grey60")
p23a <- p23a + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p23a <- p23a + coord_cartesian(xlim=c(0,100))
p23a <- p23a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold"),
                     text = element_text(size=24),
                     legend.position="none")
p23a




######
# 5% #
######

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_0.05, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - overhang
#' - permanance
#' - pond substrate
#' - amphibians
#' - pond area
#' - waterfowl
#' - max depth
#' - scrub hedge
#' - pond density
#' - inflow
#' - macrophytes
#' - fish
#' - overall terrestrial habitat
#' - terrestrial other
#' - ruderals
#' - rough grass
#' - water quality
#' - outflow
#' - woodland
#' 
 

## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 4 is optimal i.e. 4 explanatory variables. 
#' The best tree will be subjective.
#' 
#' All variables excluding pollution were included in the classification 
#' tree.
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_0.05 <- data.frame(habitat_0.05[,c(1,10,14:16,18,20:24,26:34,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_0.05 < scale(SRdist_0.05[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Permanance + Pond_Substrate + Amphibians 
                + Pond_Area + Waterfowl + Max_Depth + Scrub_Hedge
                + Pond_Density + Inflow + Macrophytes + Fish + Rough_Grass
                + Terrestrial_Other + Ruderals + Quality
                + Outflow + Woodland + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.05)

## drop1() function:

drop1(model0)


## Remove pond substrate as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Permanance + Amphibians 
                + Pond_Area + Waterfowl + Max_Depth + Scrub_Hedge
                + Pond_Density + Inflow + Macrophytes + Fish + Rough_Grass
                + Terrestrial_Other + Ruderals + Quality
                + Outflow + Woodland + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.05)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove amphibians as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Permanance + Pond_Area + Waterfowl 
                + Max_Depth + Scrub_Hedge + Pond_Density + Inflow 
                + Macrophytes + Fish + Rough_Grass
                + Terrestrial_Other + Ruderals + Quality
                + Outflow + Woodland + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.05)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model2)


## Remove ruderals as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Permanance + Pond_Area + Waterfowl 
                + Max_Depth + Scrub_Hedge + Pond_Density + Inflow 
                + Macrophytes + Fish + Rough_Grass
                + Terrestrial_Other + Quality
                + Outflow + Woodland + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.05)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has resulted in better model fit


drop1(model3)


## Remove woodland as this variable has the lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Permanance + Pond_Area + Waterfowl 
                + Max_Depth + Scrub_Hedge + Pond_Density + Inflow 
                + Macrophytes + Fish + Rough_Grass + Terrestrial_Other 
                + Quality + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.05)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of woodland has resulted in better model fit


drop1(model4)


## Remove overall terrestrial habitat as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Permanance + Pond_Area + Waterfowl 
                + Max_Depth + Scrub_Hedge + Pond_Density + Inflow 
                + Macrophytes + Fish + Rough_Grass + Terrestrial_Other 
                + Quality + Outflow,
                family = poisson,
                data = SRdist_0.05)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value 
## so removal of overall terrestrial habitat has resulted in better model fit


drop1(model5)


## Remove terrestrial other as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Permanance + Pond_Area + Waterfowl 
                + Max_Depth + Scrub_Hedge + Pond_Density + Inflow 
                + Macrophytes + Fish + Rough_Grass + Quality + Outflow,
                family = poisson,
                data = SRdist_0.05)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value 
## so removal of terrestrial other has resulted in better model fit


drop1(model6)


## Remove permanance as this variable has the lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Area + Waterfowl + Max_Depth 
                + Scrub_Hedge + Pond_Density + Inflow + Macrophytes 
                + Fish + Rough_Grass + Quality + Outflow,
                family = poisson,
                data = SRdist_0.05)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl presence has resulted in better model fit


drop1(model7)


## Remove waterfowl as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Area + Max_Depth 
                + Scrub_Hedge + Pond_Density + Inflow + Macrophytes 
                + Fish + Rough_Grass + Quality + Outflow,
                family = poisson,
                data = SRdist_0.05)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of waterfowl has improved model fit


drop1(model8)


## Remove max depth as this variable has the lowest AIC value
## (alongside macrophyte cover)

model9 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Area + Scrub_Hedge 
                + Pond_Density + Inflow + Macrophytes 
                + Fish + Rough_Grass + Quality + Outflow,
                family = poisson,
                data = SRdist_0.05)

anova(model9,model8)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with lower AIC value/less parameters to avoid
## overparameterisation
## so removal of max depth has improved model fit


drop1(model9)


## Remove macrophyte as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Area + Scrub_Hedge 
                + Pond_Density + Inflow + Fish 
                + Rough_Grass + Quality + Outflow,
                family = poisson,
                data = SRdist_0.05)

anova(model10,model9)

## p-value from LRT not significant and reduction in AIC value 
## so removal of macrophyte cover has improved model fit


drop1(model10)


## Remove fish as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Area + Scrub_Hedge 
                 + Pond_Density + Inflow + Rough_Grass 
                 + Quality + Outflow,
                 family = poisson,
                 data = SRdist_0.05)

anova(model11,model10)

## p-value from LRT not significant and slight decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of fish has improved model fit


drop1(model11)


## Remove pond density as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Area + Scrub_Hedge + Inflow 
                 + Rough_Grass + Quality + Outflow,
                 family = poisson,
                 data = SRdist_0.05)

anova(model12,model11)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has improved model fit


drop1(model12)


## Remove water quality as this variable has the lowest AIC value
## (alongside pond area)

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Area + Scrub_Hedge + Inflow 
                 + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_0.05)

anova(model13,model12)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of water quality has improved model fit


drop1(model13)


## Remove scrub/hedge as this variable has the lowest AIC value
## (alongside pond area)

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Area + Inflow + Rough_Grass 
                 + Outflow,
                 family = poisson,
                 data = SRdist_0.05)

anova(model14,model13)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of scrub/hedge has improved model fit


drop1(model14)


## Remove pond area as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Inflow + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_0.05)

anova(model15,model14)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond area has improved model fit


drop1(model15)


## Remove outflow as this variable has the lowest AIC value

model16 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Inflow + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.05)

anova(model16,model15)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of outflow has improved model fit


drop1(model16)


## Remove inflow as this variable has the lowest AIC value

model17 <- glmer(Sp_richness ~ (1|sample) + Overhang + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.05)

anova(model17,model16)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of inflow has improved model fit


drop1(model17)


## Remove rough grass as this variable has the lowest AIC value

model18 <- glmer(Sp_richness ~ (1|sample) + Overhang,
                 family = poisson,
                 data = SRdist_0.05)

anova(model18,model17)

## p-value from LRT significant and increase in AIC value (>2)
## removal of rough grass has reduced model fit



## Final model: 

SRm_0.05 <- glmer(Sp_richness ~ (1|sample) + Overhang + Rough_Grass,
                  family = poisson,
                  data = SRdist_0.05)

summary(SRm_0.05)
anova(SRm_0.05)
drop1(SRm_0.05, test = "Chi")   # for binomial/integer y models, 
                                # statistics for significance of each term

display(SRm_0.05)
se.ranef(SRm_0.05)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_0.05)

## Overdispersion test (chi-square)

1-pchisq(453.047, df=499)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_0.05)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_0.05, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_0.05)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_0.05)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_0.05$Overhang)
plot(sresid ~ SRdist_0.05$Rough_Grass)


## Check model residuals
chkres(SRm_0.05)

## Plot the fitted data against the observed data
plot(SRdist_0.05$Sp_richness ~ fitted(SRm_0.05))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_0.05)


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

hoslem.test(SRdist_0.05$Sp_richness, fitted(SRm_0.05))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: SRdist_0.05
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_0.05, newdata=SRdist_0.05, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_0.05, fit)



## PLOT 24A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_0.05$Overhang), to=max(SRdist_0.05$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range,
                 Rough_Grass=rep('Some', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_0.05, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p24a <- ggplot() + ggtitle('(a)')
p24a <- p24a + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p24a <- p24a + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p24a <- p24a + geom_jitter(aes(x=SRdist_0.05$Overhang, y=SRdist_0.05$Sp_richness, colour=SRdist_0.05$Sp_richness), cex=1, height=0.7)
p24a <- p24a + scale_colour_gradient(low="grey48", high="grey60")
p24a <- p24a + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p24a <- p24a + coord_cartesian(xlim=c(0,100))
p24a <- p24a + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p24a



## PLOT 24B: rough grass and species richness

f.Rough_Grass <- ordered(SRdist_0.05$Rough_Grass, levels = c("None","Some","Important"))
f.Rough_Grass

p24b <- ggplot(box.dat) + ggtitle('(b)')
p24b <- p24b + geom_jitter(aes(x=f.Rough_Grass, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p24b <- p24b + geom_boxplot(aes(x=f.Rough_Grass, y=fit), alpha=0.7, outlier.colour = "red")
p24b <- p24b + scale_colour_gradient(low="grey48", high="grey60")
p24b <- p24b + labs(y="Species richness", x = "Rough grass habitat") 
p24b <- p24b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p24b


grid.arrange(arrangeGrob(p24a,p24b, ncol=2), nrow=1,
             top=textGrob("Predictors of species richness with 5% sequence threshold",
                          gp=gpar(cex=3)))




#######
# 10% #
#######

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_0.1, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - overhang
#' - pond substrate
#' - max depth
#' - permanence
#' - inflow
#' - fish
#' - waterfowl
#' - macrophytes
#' - woodland
#' - pond density
#' - pond area
#' - terrestrial other
#' - rough grass
#' - ruderals
#' - water quality
#' - scrub hedge
#' - amphibians
#' - outflow
#' - overall terrestrial habitat
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 9 is optimal i.e. 9 explanatory variables. 
#' The best tree will be subjective.
#' 
#' All variables excluding pollution were included in the classification 
#' tree.
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_0.1 <- data.frame(habitat_0.1[,c(1,10,14:16,18,20:24,26:34,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_0.1 < scale(SRdist_0.1[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Substrate + Max_Depth + Permanance
                + Inflow + Fish + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge 
                + Amphibians + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.1)

## drop1() function:

drop1(model0)


## Remove pond substrate as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Fish + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge 
                + Amphibians + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.1)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove amphibians as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Fish + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.1)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model2)


## Remove permanance as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Fish 
                + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.1)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of permanance has resulted in better model fit


drop1(model3)


## Remove fish presence as this variable has the lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Woodland + Pond_Density + Pond_Area 
                + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge 
                + Outflow + Overall_Habitat_Score,
                family = poisson,
                data = SRdist_0.1)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish presence has resulted in better model fit


drop1(model4)


## Remove overall terrestrial habitat as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Woodland + Pond_Density + Pond_Area 
                + Terrestrial_Other + Macrophytes + Rough_Grass 
                + Ruderals + Quality + Scrub_Hedge + Outflow,
                family = poisson,
                data = SRdist_0.1)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value 
## so removal of overall terrestrial habitat has resulted in better model fit


drop1(model5)


## Remove scrub/hedge as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Woodland + Pond_Density + Pond_Area 
                + Terrestrial_Other + Macrophytes + Rough_Grass 
                + Ruderals + Quality + Outflow,
                family = poisson,
                data = SRdist_0.1)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value 
## so removal of scrub/hedge has resulted in better model fit


drop1(model6)


## Remove terrestrial other as this variable has the lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Woodland + Pond_Density + Pond_Area 
                + Macrophytes + Rough_Grass  + Ruderals 
                + Quality + Outflow,
                family = poisson,
                data = SRdist_0.1)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of terrestrial other has resulted in better model fit


drop1(model7)


## Remove water quality as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Woodland + Pond_Density + Pond_Area 
                + Macrophytes + Rough_Grass  + Ruderals 
                + Outflow,
                family = poisson,
                data = SRdist_0.1)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of water quality has improved model fit


drop1(model8)


## Remove ruderals as this variable has the lowest AIC value

model9 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Woodland + Pond_Density + Pond_Area 
                + Macrophytes + Rough_Grass + Outflow,
                family = poisson,
                data = SRdist_0.1)

anova(model9,model8)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has improved model fit


drop1(model9)


## Remove woodland as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Pond_Density + Pond_Area + Macrophytes 
                + Rough_Grass + Outflow,
                family = poisson,
                data = SRdist_0.1)

anova(model10,model9)

## p-value from LRT not significant and reduction in AIC value 
## so removal of woodland has improved model fit


drop1(model10)


## Remove pond density as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Max_Depth + Inflow + Waterfowl 
                 + Pond_Area + Macrophytes + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_0.1)

anova(model11,model10)

## p-value from LRT not significant and slight decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has improved model fit


drop1(model11)


## Remove max depth as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Inflow + Waterfowl + Pond_Area 
                 + Macrophytes + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_0.1)

anova(model12,model11)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of max depth has improved model fit


drop1(model12)


## Remove inflow as this variable has the lowest AIC value

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Waterfowl + Pond_Area 
                 + Macrophytes + Rough_Grass + Outflow,
                 family = poisson,
                 data = SRdist_0.1)

anova(model13,model12)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of inflow has improved model fit


drop1(model13)


## Remove outflow as this variable has the lowest AIC value

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Waterfowl + Pond_Area + Macrophytes 
                 + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.1)

anova(model14,model13)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of outflow has improved model fit


drop1(model14)


## Remove waterfowl presence as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Pond_Area + Macrophytes 
                 + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.1)

anova(model15,model14)

## p-value from LRT not significant and decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of waterfowl presence has improved model fit


drop1(model15)


## Remove pond area as this variable has the lowest AIC value

model16 <- glmer(Sp_richness ~ (1|sample) 
                 + Overhang + Macrophytes + Rough_Grass,
                 family = poisson,
                 data = SRdist_0.1)

anova(model16,model15)

## p-value from LRT not significant and minimal increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond area has improved model fit


drop1(model16)


## Remove rough grass as this variable has the lowest AIC value

model17 <- glmer(Sp_richness ~ (1|sample) + Overhang + Macrophytes,
                 family = poisson,
                 data = SRdist_0.1)

anova(model17,model16)

## p-value from LRT not significant and minimal increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of rough grass has improved model fit


drop1(model17)


## Remove macrophytes this variable has the lowest AIC value

model18 <- glmer(Sp_richness ~ (1|sample) + Overhang,
                 family = poisson,
                 data = SRdist_0.1)

anova(model18,model17)

## p-value from LRT not significant and minimal increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of macrophytes has improved model fit


## Compare to null model to see if final model has been reached

model19 <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = SRdist_0.1)

anova(model19,model18)

## p-value from LRT highly significant and model AIc increased (>2)
## therefore model with overhang is better fit to data than null model



## Final model: 

SRm_0.1 <- glmer(Sp_richness ~ (1|sample) + Overhang,
                  family = poisson,
                  data = SRdist_0.1)

summary(SRm_0.1)
anova(SRm_0.1)
drop1(SRm_0.1, test = "Chi")   # for binomial/integer y models, 
                               # statistics for significance of each term

display(SRm_0.1)
se.ranef(SRm_0.1)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_0.1)

## Overdispersion test (chi-square)

1-pchisq(438.947, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_0.1)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_0.1, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_0.1)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_0.1)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_0.1$Overhang)


## Check model residuals
chkres(SRm_0.1)

## Plot the fitted data against the observed data
plot(SRdist_0.1$Sp_richness ~ fitted(SRm_0.1))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_0.1)


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

hoslem.test(SRdist_0.1$Sp_richness, fitted(SRm_0.1))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: SRdist_0.1
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_0.1, newdata=SRdist_0.1, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_0.1, fit)



## PLOT 25A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_0.1$Overhang), to=max(SRdist_0.1$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_0.1, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p25 <- ggplot() + ggtitle('Predictors of species richness with 10% sequence threshold')
p25 <- p25 + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p25 <- p25 + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p25 <- p25 + geom_jitter(aes(x=SRdist_0.1$Overhang, y=SRdist_0.1$Sp_richness, colour=SRdist_0.1$Sp_richness), cex=1, height=0.7)
p25 <- p25 + scale_colour_gradient(low="grey48", high="grey60")
p25 <- p25 + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p25 <- p25 + coord_cartesian(xlim=c(0,100))
p25 <- p25 + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold"),
                     text = element_text(size=24),
                     legend.position="none")
p25



#######
# 30% #
#######

## Classification trees allow detailed investigation into the relative
## importance of explanatory variables

f1 <- formula(Sp_richness ~ Max_Depth + Pond_Area 
              + Pond_Density + Overhang + Macrophytes 
              + Permanance + Quality + Pond_Substrate + Inflow
              + Outflow + Pollution + Amphibians + Waterfowl + Fish 
              + Woodland + Rough_Grass + Scrub_Hedge + Ruderals 
              + Terrestrial_Other + Overall_Habitat_Score)


SR_tree <- rpart(f1, data = habitat_0.3, method = "class",
                 minsplit=5, cp = 0.001)


par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(SR_tree,uniform=TRUE, margin=0.1)
text(SR_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' Most important variables are:
#' 
#' - overhang
#' - pond substrate
#' - waterfowl
#' - permanence
#' - max depth
#' - macrophytes
#' - pond area
#' - water quality
#' - ruderals
#' - pond density
#' - scrub/hedge
#' - woodland
#' - inflow
#' - fish
#' - rough grass
#' - amphibians
#' - terrestrial other
#' 


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(SR_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))


#'
#' A tree of size of 11 is optimal i.e. 11 explanatory variables. 
#' The best tree will be subjective.
#' 
#' All variables excluding pollution, outflow, and overall terrestrial
#' habitat were included in the classification tree.
#' 


#' 
#' Now, we must standardise the dataframes containing only the variables
#' which are to be modelled.
#' 

## Create new dataframe containing only explanatory variables to be
## modelled

SRdist_0.3 <- data.frame(habitat_0.3[,c(1,10,14:16,18,20:23,26:33,97)])


## Standardise explanatory variables to have a mean of 0 and standard
## deviation of 1. This will improve convergence of fitting algorithm
## and put estimated coefficients on the same scale allowing effect 
## sizes to be compared more easily.

SRdist_0.3 < scale(SRdist_0.3[,c(2:6)])


## Apply step-wise down nested model selection using glmer

model0 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Pond_Substrate + Max_Depth + Permanance
                + Inflow + Fish + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge 
                + Amphibians,
                family = poisson,
                data = SRdist_0.3)

## drop1() function:

drop1(model0)


## Remove pond substrate as this variable has the lowest AIC value

model1 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Fish + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge 
                + Amphibians,
                family = poisson,
                data = SRdist_0.3)

anova(model1,model0)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond substrate has resulted in better model fit


drop1(model1)


## Remove amphibians as this variable has the lowest AIC value

model2 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Fish + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge,
                family = poisson,
                data = SRdist_0.3)

anova(model2,model1)

## p-value from LRT not significant and reduction in AIC value
## so removal of amphibians has resulted in better model fit


drop1(model2)


## Remove fish presence as this variable has the lowest AIC value

model3 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Ruderals + Quality + Scrub_Hedge,
                family = poisson,
                data = SRdist_0.3)

anova(model3,model2)

## p-value from LRT not significant and reduction in AIC value
## so removal of fish presence has resulted in better model fit


drop1(model3)


## Remove ruderals this variable has the lowest AIC value

model4 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Rough_Grass + Quality + Scrub_Hedge,
                family = poisson,
                data = SRdist_0.3)

anova(model4,model3)

## p-value from LRT not significant and reduction in AIC value
## so removal of ruderals has resulted in better model fit


drop1(model4)


## Remove rough grass as this variable has the lowest AIC value

model5 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Terrestrial_Other + Macrophytes 
                + Quality + Scrub_Hedge,
                family = poisson,
                data = SRdist_0.3)

anova(model5,model4)

## p-value from LRT not significant and reduction in AIC value 
## so removal of rough grass has resulted in better model fit


drop1(model5)


## Remove terrestrial other as this variable has the lowest AIC value

model6 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Macrophytes + Quality + Scrub_Hedge,
                family = poisson,
                data = SRdist_0.3)

anova(model6,model5)

## p-value from LRT not significant and reduction in AIC value 
## so removal of terrestrial other has resulted in better model fit


drop1(model6)


## Remove scrub/hedge as this variable has the lowest AIC value

model7 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Waterfowl + Woodland + Pond_Density 
                + Pond_Area + Macrophytes + Quality,
                family = poisson,
                data = SRdist_0.3)

anova(model7,model6)

## p-value from LRT not significant and reduction in AIC value
## so removal of scrub/hedge has resulted in better model fit


drop1(model7)


## Remove pond area as this variable has the lowest AIC value

model8 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Permanance
                + Inflow + Waterfowl + Woodland + Pond_Density 
                + Macrophytes + Quality,
                family = poisson,
                data = SRdist_0.3)

anova(model8,model7)

## p-value from LRT not significant and reduction in AIC value
## so removal of pond area has improved model fit


drop1(model8)


## Remove permanence as this variable has the lowest AIC value

model9 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Max_Depth + Inflow + Waterfowl 
                + Woodland + Pond_Density + Macrophytes + Quality,
                family = poisson,
                data = SRdist_0.3)

anova(model9,model8)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## In this instance, it is better to choose model with lower AIC/
## less parameters to avoid overparameterisation.
## So removal of permanence has improved model fit


drop1(model9)


## Remove max depth as this variable has the lowest AIC value

model10 <- glmer(Sp_richness ~ (1|sample)
                + Overhang + Inflow + Waterfowl 
                + Woodland + Pond_Density + Macrophytes + Quality,
                family = poisson,
                data = SRdist_0.3)

anova(model10,model9)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## In this instance, it is better to choose model with lower AIC/
## less parameters to avoid overparameterisation.
## So removal of permanence has improved model fit


drop1(model10)


## Remove pond density as this variable has the lowest AIC value

model11 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Inflow + Waterfowl + Woodland 
                 + Macrophytes + Quality,
                 family = poisson,
                 data = SRdist_0.3)

anova(model11,model10)

## p-value from LRT not significant and slight decrease in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of pond density has improved model fit


drop1(model11)


## Remove water quality as this variable has the lowest AIC value

model12 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Inflow + Waterfowl + Woodland 
                 + Macrophytes,
                 family = poisson,
                 data = SRdist_0.3)

anova(model12,model11)

## p-value from LRT not significant and slight reduction in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of water quality has improved model fit


drop1(model12)


## Remove inflow as this variable has the lowest AIC value

model13 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Waterfowl + Woodland + Macrophytes,
                 family = poisson,
                 data = SRdist_0.3)

anova(model13,model12)

## p-value from LRT not significant and minimal difference in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of inflow has improved model fit


drop1(model13)


## Remove macrophytes as this variable has the lowest AIC value

model14 <- glmer(Sp_richness ~ (1|sample)
                 + Overhang + Waterfowl + Woodland,
                 family = poisson,
                 data = SRdist_0.3)

anova(model14,model13)

## p-value from LRT not significant despite slight increase in AIC value (<2)
## better to choose model with less parameters to avoid overparameterisation
## so removal of macrophytes has improved model fit


drop1(model14)


## Remove woodland as this variable has the lowest AIC value

model15 <- glmer(Sp_richness ~ (1|sample) + Overhang + Waterfowl,
                 family = poisson,
                 data = SRdist_0.3)

anova(model15,model14)

## p-value from LRT significant and increase in AIC value (>2)
## so removal of woodland has reduced model fit



## Final model: 

SRm_0.3 <- glmer(Sp_richness ~ (1|sample) 
                 + Overhang + Waterfowl + Woodland,
                 family = poisson,
                 data = SRdist_0.3)

summary(SRm_0.3)
anova(SRm_0.3)
drop1(SRm_0.3, test = "Chi")   # for binomial/integer y models, 
                               # statistics for significance of each term

display(SRm_0.3)
se.ranef(SRm_0.3)             # levels of random factor centred around 0


## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(SRm_0.3)

## Overdispersion test (chi-square)

1-pchisq(360.7, df=497)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(SRm_0.3)

## Model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(SRm_0.3, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(SRm_0.3)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(SRm_0.3)


## Plot the residuals against each x vairable
plot(sresid ~ SRdist_0.3$Overhang)
plot(sresid ~ SRdist_0.3$Waterfowl)
plot(sresid ~ SRdist_0.3$Woodland)


## Check model residuals
chkres(SRm_0.3)

## Plot the fitted data against the observed data
plot(SRdist_0.3$Sp_richness ~ fitted(SRm_0.3))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(SRm_0.3)


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

hoslem.test(SRdist_0.3$Sp_richness, fitted(SRm_0.3))

## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT MODEL
## Obtain predicted values for the full data set: SRdist_0.1
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predictSE(SRm_0.3, newdata=SRdist_0.3, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for 
## variables in model
box.dat <- cbind(SRdist_0.3, fit)



## PLOT 26A: Overhang and species richness 

## Create a range of terrestrial overhang values which increase by 0.1984127
## to predict species richness as terrestrial overhang increases
range <- seq(from=min(SRdist_0.3$Overhang), to=max(SRdist_0.3$Overhang), by=0.1984127)


# Create new data frame where only terrestrial overhang changes
# Continuous variables are set as the mean values in this study
# Factors are set to first level in dataset

d3 <- data.frame(sample=rep('F1', length(range)),
                 Overhang=range,
                 Waterfowl=rep('Minor', length(range)),
                 Woodland=rep('Important', length(range)))


## Get predictions for this new dataset
fit <- data.frame(predictSE(SRm_0.3, newdata=d3, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.OH <- cbind(d3, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between terrestrial overhang and species richness
## with predicted values from model

p26a <- ggplot() + ggtitle('(a)')
p26a <- p26a + geom_line(aes(x=dat.OH$Overhang, y=dat.OH$fit), size = 1)
p26a <- p26a + geom_ribbon(aes(x=dat.OH$Overhang, ymin = dat.OH$cil, ymax = dat.OH$ciu), alpha = 0.25)
p26a <- p26a + geom_jitter(aes(x=SRdist_0.3$Overhang, y=SRdist_0.3$Sp_richness, colour=SRdist_0.3$Sp_richness), cex=1, height=0.7)
p26a <- p26a + scale_colour_gradient(low="grey48", high="grey60")
p26a <- p26a + labs(y = "Species richness", x = "Terrestrial overhang (%)")
p26a <- p26a + coord_cartesian(xlim=c(0,100))
p26a <- p26a + theme(panel.background = element_rect(fill = 'grey98'), 
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position="none")
p26a



## PLOT 26B: waterfowl and species richness

f.Waterfowl <- ordered(SRdist_0.3$Waterfowl, levels = c("Absent","Minor","Major"))
f.Waterfowl

p26b <- ggplot(box.dat) + ggtitle('(b)')
p26b <- p26b + geom_jitter(aes(x=f.Waterfowl, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p26b <- p26b + geom_boxplot(aes(x=f.Waterfowl, y=fit), alpha=0.7, outlier.colour = "red")
p26b <- p26b + scale_colour_gradient(low="grey48", high="grey60")
p26b <- p26b + labs(y="Species richness", x = "Waterfowl presence") 
p26b <- p26b + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p26b



## PLOT 26C: woodland and species richness

f.Woodland <- ordered(SRdist_0.3$Woodland, levels = c("None","Some","Important"))
f.Woodland

p26c <- ggplot(box.dat) + ggtitle('(c)')
p26c <- p26c + geom_jitter(aes(x=f.Woodland, y=Sp_richness, colour=Sp_richness), width=0.5, height=0.7, cex=1)
p26c <- p26c + geom_boxplot(aes(x=f.Woodland, y=fit), alpha=0.7, outlier.colour = "red")
p26c <- p26c + scale_colour_gradient(low="grey48", high="grey60")
p26c <- p26c + labs(y="Species richness", x = "Woodland") 
p26c <- p26c + theme(panel.background = element_rect(fill = 'grey98'),
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p26c


grid.arrange(arrangeGrob(p26a, ncol=1), 
             arrangeGrob(p26b,p26c, ncol=2), nrow=2,
             top=textGrob("Predictors of species richness with 30% sequence threshold",
                          gp=gpar(cex=3)))



#' ----
#' 
#' Now, test HSI score for correlation with GCN presence in different
#' threshold datasets.
#' 


################
# No threshold #
################

HSI_NT <- glmer(Sp_richness ~ (1|sample) + HSI,
                family = poisson,
                data = habitat_NT)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                  family = poisson,
                  data = habitat_NT)

anova(HSI_NT, HSInull)    

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_NT)
anova(HSI_NT)
drop1(HSI_NT, test="Chi")  # for binomial/integer y models, 
                           # statistics for significance of each term

display(HSI_NT)
se.ranef(HSI_NT)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_NT)

## Overdispersion test (chi-square)

1-pchisq(446.341, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_NT)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_NT, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_NT)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_NT)


## Plot the residuals against each x vairable
plot(sresid ~ habitat.dat$HSI)


## Check model residuals
chkres(HSI_NT)

## Plot the fitted data against the observed data
plot(habitat.dat$Sp_richness ~ fitted(HSI_NT))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_NT)



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

hoslem.test(habitat_NT$Sp_richness, fitted(HSI_NT))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 27A: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_NT$HSI), to=max(habitat_NT$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_NT, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27a <- ggplot() + ggtitle('(a) No threshold')
p27a <- p27a + geom_jitter(aes(x=habitat_NT$HSI, y=habitat_NT$Sp_richness, colour=habitat_NT$Sp_richness), height=0.7, cex=1)
p27a <- p27a + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27a <- p27a + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27a <- p27a + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27a <- p27a + coord_cartesian(xlim=c(0,1))
p27a <- p27a + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27a <- p27a + theme(panel.background = element_rect(fill = 'grey98'), 
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position="none")
p27a



#########
# 0.05% #
#########

HSI_0.05 <- glmer(Sp_richness ~ (1|sample) + HSI,
                  family = poisson,
                  data = habitat_0.0005)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = habitat_0.0005)

anova(HSI_0.05, HSInull)    

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_0.05)
anova(HSI_0.05)
drop1(HSI_0.05, test="Chi")  # for binomial/integer y models, 
                             # statistics for significance of each term

display(HSI_0.05)
se.ranef(HSI_0.05)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_0.05)

## Overdispersion test (chi-square)

1-pchisq(462.615, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_0.05)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_0.05, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_0.05)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_0.05)


## Plot the residuals against each x vairable
plot(sresid ~ habitat.dat$HSI)


## Check model residuals
chkres(HSI_0.05)

## Plot the fitted data against the observed data
plot(habitat_0.0005$Sp_richness ~ fitted(HSI_0.05))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_0.05)


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

hoslem.test(habitat_0.0005$Sp_richness, fitted(HSI_0.05))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data


## PLOT 27B: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_0.0005$HSI), to=max(habitat_0.0005$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_0.05, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27b <- ggplot() + ggtitle('(b) 0.05%')
p27b <- p27b + geom_jitter(aes(x=habitat_0.0005$HSI, y=habitat_0.0005$Sp_richness, colour=habitat_0.0005$Sp_richness), height=0.7, cex=1)
p27b <- p27b + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27b <- p27b + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27b <- p27b + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27b <- p27b + coord_cartesian(xlim=c(0,1))
p27b <- p27b + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27b <- p27b + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p27b



########
# 0.1% #
########

HSI_0.1 <- glmer(Sp_richness ~ (1|sample) + HSI,
                  family = poisson,
                  data = habitat_0.001)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = habitat_0.001)

anova(HSI_0.1, HSInull)    

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_0.1)
anova(HSI_0.1)
drop1(HSI_0.1, test="Chi")  # for binomial/integer y models, 
                            # statistics for significance of each term

display(HSI_0.1)
se.ranef(HSI_0.1)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_0.1)

## Overdispersion test (chi-square)

1-pchisq(464.718, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_0.1)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_0.1, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_0.1)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_0.1)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.001$HSI)


## Check model residuals
chkres(HSI_0.1)

## Plot the fitted data against the observed data
plot(habitat_0.001$Sp_richness ~ fitted(HSI_0.1))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_0.1)



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

hoslem.test(habitat_0.001$Sp_richness, fitted(HSI_0.1))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 27C: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_0.001$HSI), to=max(habitat_0.001$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_0.1, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27c <- ggplot() + ggtitle('(c) 0.1%')
p27c <- p27c + geom_jitter(aes(x=habitat_0.001$HSI, y=habitat_0.001$Sp_richness, colour=habitat_0.001$Sp_richness), height=0.7, cex=1)
p27c <- p27c + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27c <- p27c + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27c <- p27c + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27c <- p27c + coord_cartesian(xlim=c(0,1))
p27c <- p27c + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27c <- p27c + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p27c



########
# 0.5% #
########

HSI_0.5 <- glmer(Sp_richness ~ (1|sample) + HSI,
                 family = poisson,
                 data = habitat_0.005)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = habitat_0.005)

anova(HSI_0.5, HSInull)    

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_0.5)
anova(HSI_0.5)
drop1(HSI_0.5, test="Chi")  # for binomial/integer y models, 
                            # statistics for significance of each term

display(HSI_0.5)
se.ranef(HSI_0.5)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_0.5)

## Overdispersion test (chi-square)

1-pchisq(563.455, df=501)   # final model overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_0.5)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_0.5, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_0.5)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_0.5)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.005$HSI)


## Check model residuals
chkres(HSI_0.5)

## Plot the fitted data against the observed data
plot(habitat_0.005$Sp_richness ~ fitted(HSI_0.5))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_0.5)



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

hoslem.test(habitat_0.005$Sp_richness, fitted(HSI_0.5))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 27D: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_0.005$HSI), to=max(habitat_0.005$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_0.5, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27d <- ggplot() + ggtitle('(d) 0.5%')
p27d <- p27d + geom_jitter(aes(x=habitat_0.005$HSI, y=habitat_0.005$Sp_richness, colour=habitat_0.005$Sp_richness), height=0.7, cex=1)
p27d <- p27d + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27d <- p27d + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27d <- p27d + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27d <- p27d + coord_cartesian(xlim=c(0,1))
p27d <- p27d + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27d <- p27d + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p27d



######
# 1% #
######

HSI_1 <- glmer(Sp_richness ~ (1|sample) + HSI,
                 family = poisson,
                 data = habitat_0.01)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = habitat_0.01)

anova(HSI_1, HSInull) 

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_1)
anova(HSI_1)
drop1(HSI_1, test="Chi")  # for binomial/integer y models, 
                          # statistics for significance of each term

display(HSI_1)
se.ranef(HSI_1)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_1)

## Overdispersion test (chi-square)

1-pchisq(550.847, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_1)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_1, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_1)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_1)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.01$HSI)


## Check model residuals
chkres(HSI_1)

## Plot the fitted data against the observed data
plot(habitat_0.01$Sp_richness ~ fitted(HSI_1))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_1)



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

hoslem.test(habitat_0.01$Sp_richness, fitted(HSI_1))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 27E: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_0.01$HSI), to=max(habitat_0.01$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_1, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27e <- ggplot() + ggtitle('(e) 1%')
p27e <- p27e + geom_jitter(aes(x=habitat_0.01$HSI, y=habitat_0.01$Sp_richness, colour=habitat_0.01$Sp_richness), height=0.7, cex=1)
p27e <- p27e + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27e <- p27e + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27e <- p27e + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27e <- p27e + coord_cartesian(xlim=c(0,1))
p27e <- p27e + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27e <- p27e + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p27e



######
# 5% #
######

HSI_5 <- glmer(Sp_richness ~ (1|sample) + HSI,
               family = poisson,
               data = habitat_0.05)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = habitat_0.05)

anova(HSI_5, HSInull) 

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_5)
anova(HSI_5)
drop1(HSI_5, test="Chi")  # for binomial/integer y models, 
                          # statistics for significance of each term

display(HSI_5)
se.ranef(HSI_5)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_5)

## Overdispersion test (chi-square)

1-pchisq(477.305, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_5)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_5, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_5)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_5)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.05$HSI)


## Check model residuals
chkres(HSI_5)

## Plot the fitted data against the observed data
plot(habitat_0.05$Sp_richness ~ fitted(HSI_5))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_5)



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

hoslem.test(habitat_0.05$Sp_richness, fitted(HSI_5))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 27F: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_0.05$HSI), to=max(habitat_0.05$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_5, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27f <- ggplot() + ggtitle('(f) 5%')
p27f <- p27f + geom_jitter(aes(x=habitat_0.05$HSI, y=habitat_0.05$Sp_richness, colour=habitat_0.05$Sp_richness), height=0.7, cex=1)
p27f <- p27f + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27f <- p27f + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27f <- p27f + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27f <- p27f + coord_cartesian(xlim=c(0,1))
p27f <- p27f + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27f <- p27f + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p27f



#######
# 10% #
#######

HSI_10 <- glmer(Sp_richness ~ (1|sample) + HSI,
               family = poisson,
               data = habitat_0.1)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = habitat_0.1)

anova(HSI_10, HSInull)

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_10)
anova(HSI_10)
drop1(HSI_10, test="Chi")  # for binomial/integer y models, 
                           # statistics for significance of each term

display(HSI_10)
se.ranef(HSI_10)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_10)

## Overdispersion test (chi-square)

1-pchisq(458.423, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_10)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_10, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_10)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_10)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.1$HSI)


## Check model residuals
chkres(HSI_10)

## Plot the fitted data against the observed data
plot(habitat_0.1$Sp_richness ~ fitted(HSI_10))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_10)



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

hoslem.test(habitat_0.1$Sp_richness, fitted(HSI_10))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 27G: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_0.1$HSI), to=max(habitat_0.1$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_10, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27g <- ggplot() + ggtitle('(g) 10%')
p27g <- p27g + geom_jitter(aes(x=habitat_0.1$HSI, y=habitat_0.1$Sp_richness, colour=habitat_0.1$Sp_richness), height=0.7, cex=1)
p27g <- p27g + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27g <- p27g + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27g <- p27g + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27g <- p27g + coord_cartesian(xlim=c(0,1))
p27g <- p27g + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27g <- p27g + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p27g



#######
# 30% #
#######

HSI_30 <- glmer(Sp_richness ~ (1|sample) + HSI,
                family = poisson,
                data = habitat_0.3)


## Compare to null model

HSInull <- glmer(Sp_richness ~ (1|sample) + 1,
                 family = poisson,
                 data = habitat_0.3)

anova(HSI_30, HSInull)

## LRT p-value significant and AIC value of null model increased
## Therefore, HSI score is an appropriate explanatory variable


## Model output:
summary(HSI_30)
anova(HSI_30)
drop1(HSI_30, test="Chi")  # for binomial/integer y models, 
                           # statistics for significance of each term

display(HSI_30)
se.ranef(HSI_30)           # levels of random factor centred around 0



## Check for model overdispersion 
## Obtain residual deviance and degrees of freedom
overdisp.glmer(HSI_30)

## Overdispersion test (chi-square)

1-pchisq(387.717, df=501)   # final model not overdispersed


## summary() can give inflated value for model residual deviance so usual
## methods of calculating residual deviance can be unreliable for GLMMs

## Also use customised function from Rob Thomas & the Guidebook team:
## 'Data Analysis with R Statistical Software'
overdisp_fun(HSI_30)

## Using custom function, model is not overdispersed
## Ratio reported is equivalent to overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion


## Perform model validation checks
## Fit model using REML and check normalised residuals for normal distribution
sresid <- resid(HSI_30, type = "pearson")
hist(sresid)

## Plot residuals against fitted values
fits <- fitted(HSI_30)
plot(sresid ~ fits)

## Plot QQ plot and heteroscedasticity plot
mcp.fnc(HSI_30)


## Plot the residuals against each x vairable
plot(sresid ~ habitat_0.3$HSI)


## Check model residuals
chkres(HSI_30)

## Plot the fitted data against the observed data
plot(habitat_0.3$Sp_richness ~ fitted(HSI_30))


## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(HSI_30)



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

hoslem.test(habitat_0.3$Sp_richness, fitted(HSI_30))


## Model fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed 
## data



## PLOT 27H: HSI score and species richness

## Create a range of HSI score values which increase by 0.00142
## to predict species richness as HSI score increases
range <- seq(from=min(habitat_0.3$HSI), to=max(habitat_0.3$HSI), by=0.00142)


## Create new data frame where only HSI score changes
## Continuous variables are set as the mean values in this study
## Factors are set to first level in dataset

d1 <- data.frame(sample=rep('F1', length(range)),
                 HSI=range)


## Get predictions for this new dataset
fit <- data.frame(predictSE(HSI_30, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.HSI <- cbind(d1, fit, ciu, cil) # This is now the data frame


## Plot showing relationship between HSI score and species richness
## with predicted values from model

p27h <- ggplot() + ggtitle('(h) 30%')
p27h <- p27h + geom_jitter(aes(x=habitat_0.3$HSI, y=habitat_0.3$Sp_richness, colour=habitat_0.3$Sp_richness), height=0.7, cex=1)
p27h <- p27h + geom_line(aes(x=dat.HSI$HSI, y=dat.HSI$fit), size = 1)
p27h <- p27h + geom_ribbon(aes(x=dat.HSI$HSI, ymin = dat.HSI$cil, ymax = dat.HSI$ciu), alpha = 0.25)
p27h <- p27h + scale_colour_gradient(low="goldenrod", high="lightgoldenrod")
p27h <- p27h + coord_cartesian(xlim=c(0,1))
p27h <- p27h + labs(y = "Species richness", x = "Habitat Suitability Index (HSI) score")
p27h <- p27h + theme(panel.background = element_rect(fill = 'grey98'), 
                     plot.title =element_text(face="bold", hjust=0),
                     text = element_text(size=24),
                     legend.position="none")
p27h


grid.arrange(arrangeGrob(p27a,p27b,p27c, ncol=3), 
             arrangeGrob(p27d,p27e,p27f, ncol=3),
             arrangeGrob(p27g,p27h, ncol=2),nrow=3,
             top=textGrob("Relationship between habitat suitability and vertebrate species richness at different sequence thresholds",
                          gp=gpar(cex=3)))

