#' ---
#' Title: "Great crested newt (*Triturus cristatus*) environmental DNA (eDNA) metabarcoding data analysis"
#' Author: "Lynsey Rebecca Harper"
#' Date: "5th December 2016"
#' ---
#' 
#' 
#' A total of 532 samples previously analysed by qPCR for great crested newt 
#' (GCN) were re-analysed by eDNA metabarcoding to compare method sensitivity 
#' and evaluate a method for detection of vertebrate communities in 
#' freshwater ponds.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. Then,
#' load functions for calculating kappa coefficient, plotting model 
#' residuals and testing model fit.
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
       "AICcmodavg","gtools","xlsx","reshape2","plyr","dplyr","arm",
       "RVAideMemoire","ResourceSelection","bbmle","RColorBrewer", "MuMIn")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://star-www.st-andrews.ac.uk/cran/")

lapply(p, require, character.only = TRUE)


## Load custom functions
f <- c("kappa.coefficient.function.R","CheckResidsFunction.R",
       "OverdispersalFunction.R","CheckConvergenceFunction.R")

lapply(f, source)


#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()


#'
#' ## Metabarcoding dataset
#' 
#' Next, examine the raw data output from metaBEAT.
#' 

raw.data <- read.csv("../Data/12S-trim30_crop110_min90_merge-forwonly_nonchimera_c0.97cov5_blast98-id1-by-taxonomy-readcounts.blast.csv", 
                     header=TRUE)


summary(raw.data[,c(1:5,length(raw.data))])
head(raw.data)
names(raw.data)
str(raw.data)


#'
#' The raw output from metaBEAT is not particularly suited to applying
#' statistical analysis in R.
#' 
#' The file was manipulated manually in Microsoft excel:
#' 
#' - Samples that were not from FERA or Nottingham removed (Joe's and empty wells)
#' - top row removed to make column headers top row.
#' - first column renamed "Tax_assignment"
#' - removed taxonomy column 
#' - corrected samples which were known to have wrong tags placed on them for
#'   sequencing causing environmental samples to be labelled controls and vice
#'   versa according to sample sheet.
#' 
#' 
#' A similar file was prepared for the metaBEAT output from the 
#' unassigned sequences from the original BLAST. A separate BLAST 
#' was performed for the unassigned sequences against the full 
#' NCBI nucleotide database.
#' 


## Reimport the edited file:

origBLAST <- read.csv("../Data/metabarcoding_PrepOutput.csv", header=TRUE)

summary(origBLAST[,1:6])
head(origBLAST)
names(origBLAST)
str(origBLAST)


## Import the prepared unassigned file:

unassBLAST <- read.csv("../Data/Unassigned_PrepOutput.csv", header=TRUE)

summary(unassBLAST[,1:6])
head(unassBLAST)
names(unassBLAST)
str(unassBLAST)


#'
#' I need to do some manipulation to combine these datasets and get
#' the final dataset.
#' 
#' First, create dataframes with only the taxonomic assignment and
#' reads.
#' 


## Remove '.blast' after sample ID in both data frames
colnames(unassBLAST) <- gsub('.blast', '', colnames(unassBLAST))
colnames(origBLAST) <- gsub('.blast', '', colnames(origBLAST))


## Bind data frames
merged_df <- rbind(origBLAST, unassBLAST)

## Merge read counts by taxonomic assignment
## Anything with the same name in each data frame will be merged and
## anything new in the unassigned BLAST added to the original BLAST
## data frame
mergedBLAST <- ddply(merged_df, .(Tax_assignment), numcolwise(sum))


## Remove '_.nc' after sample ID
colnames(mergedBLAST) <- gsub('_.nc', '', colnames(mergedBLAST))


## Export as excel file
write.xlsx(mergedBLAST, "../Data/MergedBLAST.xlsx", row.names=FALSE)


#' 
#' Now, we need to further refine the metabarcoding dataset.
#' 
#' 1. Any spurious species must be removed.
#' 2. Any Genus or Family assignments containing only one species in 
#'    the UK must be changed to that species.
#' 3. Likewise, for any species assignment which is the only species in 
#'    the UK and also has genus/family reads, read counts from all
#'    assignments must be merged.
#'  
#' A record of these changes will be kept as these changes are made.
#'      


## Inspect new dataframe
summary(mergedBLAST[,1:6])
summary(mergedBLAST)
names(mergedBLAST)
str(mergedBLAST)


#'
#' First, remove spurious assignments. This includes bacteria and 
#' invertebrates that may in fact be present in ponds but should not 
#' have been amplified by the vertebrate 12S rRNA primers. Spurious
#' assignments are as follows:
#' 
#' - Branchiopoda (invert)
#' - *Cloen dipterum* (invert)
#' - *Ceratocystis* (Fungi)
#' - *Daphnia* (invert)
#' - Fungi
#' - Gadidae/Gadiformes (marine cod)
#' - Loricariidae (catfish)
#' - Nostocales (Cyanobacteria)
#' - *Ovis canadensis* (Bighorn sheep, not present in UK)
#' - *Pan* (probably misassigned human reads but can't reliably be included)
#' - *Saurida* (Lizardfishes)
#' - uncultured bacterium
#' - uncultured *Oceanospirillales* bacterium
#' 

true_assign <- mergedBLAST[-c(92,95,98,100,102:104,106,109,111:112,
                               116,118:119),]


#'
#' Now, change genus/family level assignments where only one UK species is
#' contained within that genus/family to species level. Assignments to be
#' changed are as follows:
#' 
#' - *Bos* = *Bos taurus*
#' - *Canis* = *Canis lupus*
#' - *Pungitius* = *Pungitius pungitius*
#' - *Strix* = *Strix aluco*
#' - Hominidae = *Homo sapiens*
#' - Cichlidae = *Rhamphochromis esox*
#' 
#' Additionally, correct any species which have clearly been misassigned.
#' Given the restricted distribution of wildcat to Scotland and that 
#' eDNA samples were collected in Lincolnshire, Cheshire or Kent, we 
#' can safely say that this assignment is actually *Felis catus*.
#' 
#' 
#' These species level assignments already exist in the metabarcoding
#' dataset, therefore we can simply merge the genus/family read counts
#' with the species read counts.
#' 

## Rename the genus/family assignments 
true_assign$Tax_assignment <- as.character(true_assign$Tax_assignment)

true_assign[32, "Tax_assignment"] <- "Felis_catus"
true_assign[91, "Tax_assignment"] <- "Bos_taurus"
true_assign[93, "Tax_assignment"] <- "Canis_lupus"
true_assign[95, "Tax_assignment"] <- "Rhamphochromis_esox"
true_assign[98, "Tax_assignment"] <- "Homo_sapiens"
true_assign[104, "Tax_assignment"] <- "Pungitius_pungitius"
true_assign[105, "Tax_assignment"] <- "Strix_aluco"


## Now merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
vert_seq <- ddply(true_assign, .(Tax_assignment), numcolwise(sum))



#'
#' ## Clean up dataset
#' 
#' Now the data is in a form that can be manipulated easily, I must tidy the
#' metabarcoding dataset by: applying species-specific false positive thresholds 
#' i.e. the frequency required to remove a given species from the positive 
#' control as only cichlid DNA should be present in positive controls. 
#' 
#' Arguably, this is a more effective way to deal with contamination and 
#' false positives in metabarcoding as negative controls have no template 
#' DNA for contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially and paints a false picture of contamination.
#' 
#' A false positive threshold of 0.05% removes 95% contamination from my 
#' positive controls and ensures no GCN is detected in positive controls
#' (see excel spreadsheet - threshold analysis).
#' 
#' While I could apply this threshold to my entire dataset, this will not 
#' remove major contaminants, such as human. Conversely, a species-specific
#' threshold for human would ensure this contamination is also removed from
#' the positive PCR control and presumably in my eDNA samples. However, I 
#' will have to accept some contamination even with application of false 
#' positive threshold as there may have been environmental contamination 
#' during sample collection.
#' 
#' 
#' First, calculate the species specific false positive thresholds using 
#' the positive PCR contols. Anything not cichlid is a false positive. 
#' Create a dataframe containing only read counts and taxonomic assignment 
#' for positive controls (last 114 samples in data frame).
#' 

## Create data frame
positives <- subset(vert_seq[,648:761])
head(positives[,1:6])


## Calculate total number of reads in each positive control
pos_total <- rbind(positives, colSums(positives))

## Create new dataframe containing the frequency of reads in each sample
pos_freq <- pos_total/c(pos_total[100,])


## Add new column to this dataframe containing the maximum frequency for
## each assignment across all positive controls.
pos_freq[, "Species_Threshold"] <- apply(pos_freq, 1, max)

## Check maxmimum frequency column has been created properly
head(pos_freq[,110:115])


## Add taxonomic assignment to dataframe
species <- data.frame(vert_seq$Tax_assignment)
colnames(species) <- "Tax_assignment"

total <- data.frame("Total")
colnames(total) <- "Tax_assignment"

species <- rbind(species, total)
pos_freq <- cbind(species, pos_freq)


## Finally, must manually change the species threshold value for cichlid to 0
## as this will be a true contaminant in the dataset and no false positive 
## threshold required.

pos_freq$Species_Threshold[83]=0

## Check this manual edit has occurred
head(pos_freq[83,112:116])



#'
#' Now, I must use this dataframe to apply the species specific false 
#' positive thresholds to a dataframe containing the read counts for all
#' samples.
#' 
#' A similar workflow to as above will be followed to calculate the total
#' read counts, calculate frequency of each taxonomic in an assignment,
#' and then apply a condition to the dataframe to set any assignments to
#' a species below the false positive threshold to 0.
#' 


## Create new dataframe to contain only samples and read counts
reads <- subset(vert_seq[,c(2:761)])

## Calculate total number of reads for each sample
counts <- rbind(reads, colSums(reads))

## Create new dataframe containing the frequency of reads in each sample
freq <- counts/c(counts[100,])


## Combine read frequencies with species names
species <- data.frame(vert_seq$Tax_assignment)
colnames(species) <- "Tax_assignment"

total <- data.frame("Total")
colnames(total) <- "Tax_assignment"

species <- rbind(species, total)
freq <- cbind(species, freq)


## Add the column Species_Threshold from the pos_freq dataframe to the 
## freq dataframe containing sequence frequencies in all samples.

FP_df <- cbind(freq, pos_freq$Species_Threshold)
names(FP_df)[762] <- "Species_Threshold"

## Check this has worked
head(FP_df[,760:762])


## Replace all assignments that are less than or equal to the false positive 
## threshold for a given species with zero
## Don't include negative controls as any reads in these constitutes 
## contamination
FP_df[,c(2:509,624:761)][FP_df[c(2:509,624:761)] <= FP_df[,762]] <- 0

## Check this has worked
head(FP_df[,750:762])


#'
#' Now, I need to get the dataframe back into read count format with the
#' taxonomic assignment for further analyses.
#' 


## Remove last row containing frequencies and add the total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.

sample_freq <- FP_df[1:99,1:761]
total_counts <- counts[100,]

total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

conversion <- smartbind(sample_freq, total_counts)

final_df <- conversion*c(conversion[100,])


## Remove 'Total' row and empty Tax_assignment row and add metadata to final 
## dataframe
final_df <- final_df[-c(100), -c(1)]
final_df <- cbind(vert_seq[,1], final_df)
names(final_df)[1] <- "Tax_assignment"


## Replace any NA values in dataframe with 0
final_df[is.na(final_df)] <- 0

  
## copy the current version to new dataframe to preserve it for later
vert_seq_final <- final_df


## Export as excel spreadsheet to data working directory for visual inspection
## and comparison with initial input file
write.xlsx(vert_seq_final, "../Data/ThresholdsApplied.xlsx", row.names=FALSE)



#' 
#' The species specific false positive thresholds have been applied 
#' correctly, thus we now have the final metabarcoding dataset.
#' 
#' The output file was manually edited in excel to add additional 
#' information for each taxonomic assignment.
#' 
#' - changed name of first column to 'Tax_assignment'
#' - created column detailing level of final assignment (e.g. genus, family)
#' - created column to specify whether a sample is positive control or eDNA
#' - created column to specify what family each assignment belongs to 
#' ("other" for assignments higher than family)
#' 


## Import prepared output file

prep_df <- read.csv("../Data/FinalOutput_Prep.csv", header=TRUE)


## Inspect new dataframe
summary(prep_df[,1:6])
summary(prep_df)
names(prep_df)
str(prep_df)


#'
#' ## Summary Statistics
#' 
#' Create another data frame for presence/absence detection of species in
#' eDNA samples (1 = presence). 
#' 
#' Exclude cichlid assignment as this is the positive control. 
#' Exclude positive and negative PCR controls.
#' 

vert_seq_binary <- data.frame(prep_df[-c(83), c(1:512,627:650)])
vert_seq_binary[,5:536][vert_seq_binary[,5:536] > 0] <- 1


## total up number of samples each assignment was detected in
vert_seq_binary$Samples <- rowSums(vert_seq_binary[,5:536])
head(vert_seq_binary[,535:537])


## count number of assignments
vert_seq_binary$Count <- ifelse((vert_seq_binary$Samples==0),0,1)
head(vert_seq_binary[,535:538])


## print a summary table
summary <- vert_seq_binary[,c(1:4,537:538)]
summary


#'
#' ### Level of Assignment
#' 

## Examine levels of assignment
assign_levs <- ddply(vert_seq_binary[,c(1:4,537:538)], .(AssignmentLevel), numcolwise(sum))
assign_levs


#'
#' ### Families
#' 

## Vertebrate families
families <- ddply(vert_seq_binary[,c(1:4,537:538)], .(Family), numcolwise(sum))
families


#'
#' ### Species
#'

## Vertebrate species
vert_seq_binary <- vert_seq_binary[order(vert_seq_binary$AssignmentLevel),]
species <- ddply(vert_seq_binary[c(36:98),c(1:4,537:538)], .(Tax_assignment), numcolwise(sum))
species



#'
#' ## Summary plots
#' 
#' Plot sequence read distribution across all levels of taxonomic assignment
#' and the proportion of species detected in ponds.
#' 


## Heatmap of sequence reads

## First order data frame by taxonomic level
OrderedAssignment <- ordered(prep_df$AssignmentLevel, levels = c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Unassigned"))
prep_df <- prep_df[order(OrderedAssignment),]

## Created modified dataframe containing only necessary data columns
mod.df <- prep_df[,c(1:2,5:513,627:650)]

## Order row indexes
rownames(mod.df) <- 1:nrow(mod.df)

## Remove underscore from species names
mod.df$Tax_assignment <- gsub('_', ' ', mod.df$Tax_assignment)

## Create empty column in dataframe and label all assignments with
## which vertebrate group they belong to
mod.df$Group <- NA

mod.df[c(1,3,9,12,17,21:23,26:27,33,40,45,47,49,56,59:60,64,67:68,
         75:76,80:81,83,86:87,93,96),536] <- "Bird"

mod.df[c(8,30:31,44,51,62,74,78,88,91),536] <- "Amphibian"

mod.df[c(2,5:6,11,14:16,19,24:25,41,46,50,53:54,58,66,69,82,84:85,
         89:90,92,94),536] <- "Fish"

mod.df[c(4,7,10,13,18,20,28:29,32,34:39,42:43,48,52,55,57,61,63,
         65,70:73,77,79,95,97),536] <- "Mammal"

mod.df[c(98:99),536] <- "higherassign"


## Order dataframe by vertebrate group and taxonomic assignment
heatmap_order <- mod.df[order(mod.df$Group,mod.df$AssignmentLevel),]


## Create new dataframe which contains data in an appropriate format
## for a heatmap
heatmap.df <- melt(heatmap_order, id = c("Tax_assignment","AssignmentLevel","Group"))
heatmap.df$Tax_assignment <- factor(heatmap.df$Tax_assignment, 
                                    levels = heatmap.df$Tax_assignment[order(heatmap.df$AssignmentLevel)])


## Subset data by vertebrate groups
Fish <- subset(heatmap.df, Group == "Fish")
Fish$AssignmentLevel <- factor(Fish$AssignmentLevel, levels = c("Species","Genus","Family","Order"))

Amphib <- subset(heatmap.df, Group == "Amphibian")
Bird <- subset(heatmap.df, Group == "Bird")
Mammal <- subset(heatmap.df, Group == "Mammal")
high_assign <- subset(heatmap.df, Group == "higherassign")


hm1 <- ggplot(Fish, aes(variable, Tax_assignment))
hm1 <- hm1 + ggtitle("(a)")
hm1 <- hm1 + geom_tile(aes(fill = value))
hm1 <- hm1 + scale_fill_gradient("Reads", limits=c(0,100000), 
                               low = "grey100", high = "black")
hm1 <- hm1 + labs(x="", y="")
hm1 <- hm1 + theme(panel.background = element_rect(fill = 'white'),
                 panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 text = element_text(size=20),
                 legend.key.size = unit(2, 'lines'),
                 legend.position="none")
hm1


hm2 <- ggplot(Amphib, aes(variable, Tax_assignment))
hm2 <- hm2 + ggtitle("(b)")
hm2 <- hm2 + geom_tile(aes(fill = value))
hm2 <- hm2 + scale_fill_gradient("Reads", limits=c(0,100000), 
                                 low = "grey100", high = "black")
hm2 <- hm2 + labs(x="", y="")
hm2 <- hm2 + theme(panel.background = element_rect(fill = 'white'),
                   panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.ticks = element_blank(), 
                   axis.text.x = element_blank(),
                   text = element_text(size=20),
                   legend.key.size = unit(2, 'lines'),
                   legend.position="none")
hm2


hm3 <- ggplot(Bird, aes(variable, Tax_assignment))
hm3 <- hm3 + ggtitle("(c)")
hm3 <- hm3 + geom_tile(aes(fill = value))
hm3 <- hm3 + scale_fill_gradient("Reads", limits=c(0,100000), 
                                 low = "grey100", high = "black")
hm3 <- hm3 + labs(x="", y="")
hm3 <- hm3 + theme(panel.background = element_rect(fill = 'white'),
                   panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.ticks = element_blank(), 
                   axis.text.x = element_blank(),
                   text = element_text(size=20),
                   legend.key.size = unit(2, 'lines'),
                   legend.position="none")
hm3


hm4 <- ggplot(Mammal, aes(variable, Tax_assignment))
hm4 <- hm4 + ggtitle("(d)")
hm4 <- hm4 + geom_tile(aes(fill = value))
hm4 <- hm4 + scale_fill_gradient("Reads", limits=c(0,100000), 
                                 low = "grey100", high = "black")
hm4 <- hm4 + labs(x="", y="")
hm4 <- hm4 + theme(panel.background = element_rect(fill = 'white'),
                   panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.ticks = element_blank(), 
                   axis.text.x = element_blank(),
                   text = element_text(size=20),
                   legend.key.size = unit(2, 'lines',),
                   legend.position="none")
hm4


hm5 <- ggplot(high_assign, aes(variable, Tax_assignment))
hm5 <- hm5 + ggtitle("(e)")
hm5 <- hm5 + geom_tile(aes(fill = value))
hm5 <- hm5 + scale_fill_gradient("Reads", limits=c(0,100000), 
                                 low = "grey100", high = "black")
hm5 <- hm5 + labs(x="", y="")
hm5 <- hm5 + theme(panel.background = element_rect(fill = 'white'),
                   panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                   plot.title = element_text(face="bold", hjust=0),
                   axis.ticks = element_blank(), 
                   axis.text.x = element_blank(),
                   text = element_text(size=20),
                   legend.key.size = unit(2, 'lines'))
hm5


## ALL PLOTS with common legend and common axes

source("multiplot_function.R")

hm <- multiplot(hm1,hm4,hm2,hm5,hm3, cols=3, 
                labs=list("eDNA samples (N=532)", 
                          "Taxonomic Assignment"),
                cex=25)



heatmap_order$total <- rowSums(heatmap_order[,3:535])
install.packages("metacoder")
library(metacoder)
heat_tree(heatmap_order, node_size=total, 
          node_label=Tax_assignment, node_color=total)




## Plot proportion of species across ponds
## Remove underscores from species names
species$Tax_assignment <- gsub('_', ' ', species$Tax_assignment)
species$Proportion <- species$Samples/532

sb <- ggplot(species, aes(x=Tax_assignment, y=Proportion))
sb <- sb + geom_bar(stat="identity", fill="white", colour="black") + coord_flip()
sb <- sb + scale_y_continuous(labels = percent_format(),
                              expand = c(0,0), limits = c(0,1))
sb <- sb + labs(x = "Species", y = "eDNA Samples (N = 532)")
sb <- sb + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold"),
                 legend.title = element_blank(),
                 text = element_text(size=20))
sb



#'
#' ## Statistical Analysis
#'
#' ### Dataframe
#' 
#' Here, I'll do some further tidying up to the dataset imported
#' earlier, prep_df. 
#' 
#' First, I'll remove assignment levels higher than species: these don't 
#' provide detailed information on biodiversity unless a detected genus or
#' family only contains one species in the UK which instances of have 
#' already been identified and changed to species.
#' 

vert_detections <- subset(prep_df, prep_df$AssignmentLevel=="Species")
vert_reads <- vert_detections[,c(1,5:length(vert_detections))]


#'
#' Next, transpose the dataframe.
#' 

## prepare the assignment names and transpose the dataframe
rownames(vert_reads) <- vert_reads[,1]
vert_reads <- vert_reads[,2:761]
sample_reads <- data.frame(t(vert_reads))


## move the sample names into the first column
sample_reads$sample <- rownames(sample_reads)
sample_reads <- sample_reads[,c(length(sample_reads),1:(length(sample_reads)-1))]
row.names(sample_reads) <- NULL


## Create new data frame containing only eDNA samples and new data frame
## containing only PCR controls
eDNA_only <- data.frame(sample_reads[c(1:508,623:646),])
controls_only <- data.frame(sample_reads[c(509:622,647:760),])



#'
#' The eDNA_only dataframe will be used for statistical analysis comparing
#' eDNA methods for GCN detection and GCN pond occupancy.
#' 
#' The controls_only dataframe will be used to examine contamination in the
#' dataset.
#' 
#' Now, I can import information about the PCR and qPCR corresponding to 
#' each sample ID. This information is held in a file called 
#' "GCNmetabarcoding_metadata.csv".
#' 


## Import metadata for eDNA samples
eDNA_metadata <- read.csv("../Data/metabarcoding_metadata.csv", header=TRUE)

## Import metadata for PCR controls
controls_metadata <- read.csv("../Data/Controls_metadata.csv", header=TRUE)


#'
#' I want to merge the metadata and reads for controls and eDNA samples 
#' respectively together using the 'sample' column from the 
#' metadata dataframe and the 'Sample' column from the reads dataframe.
#'

## the sample columns need to match type, at the moment one is a factor
eDNA_metadata$sample <- as.character(eDNA_metadata$sample)
controls_metadata$sample <- as.character(controls_metadata$sample)

eDNA_reads <- merge(eDNA_metadata,eDNA_only, all=TRUE)
head(eDNA_reads[,1:6])

control_reads <- merge(controls_metadata,controls_only, all=TRUE)
head(control_reads[,1:6])


## Determine whether a control is contaminated or not using read counts
control_reads$contamination <- ifelse(rowSums(control_reads[,8:70]) > 0, "Yes", "No")




#' ---
#' 
#' # Statistical Analysis
#' 
#' ## 1) Controls
#' 
#' The key things to investigate here are:
#' 
#' - is GCN DNA present in negative controls after false positive threshold applied?
#' - how much cichlid contamination is present in negative controls and eDNA samples?
#' - How much human DNA is present in negative controls and eDNA samples?
#' - what other contamination remains in negative controls after thresholds applied?
#' - did the addition of mineral oil to PCR reactions included on the second sequencing run reduce contamination?
#' 


## Create another data frame for presence/absence detection of species in
## negative PCR controls (1=presence)
## Exclude eDNA samples
control_seq_binary <- data.frame(prep_df[, c(1:4,513:626)])
control_seq_binary[,5:118][control_seq_binary[,5:118] > 0] <- 1


## total up number of samples each assignment was detected in
control_seq_binary$Samples <- rowSums(control_seq_binary[,5:118])
head(control_seq_binary[,115:119])


## count number of assignments
control_seq_binary$Count <- ifelse((control_seq_binary$Samples==0),0,1)
head(control_seq_binary[,115:120])


## print a summary table
controls_summary <- control_seq_binary[,c(1:4,119:120)]
controls_summary


#'
#' ### Level of Assignment
#' 

## Examine levels of assignment
control_assign_levs <- ddply(control_seq_binary[,c(1:4,119:120)], .(AssignmentLevel), numcolwise(sum))
control_assign_levs


#'
#' ### Families
#' 

## Vertebrate families
control_families <- ddply(control_seq_binary[,c(1:4,119:120)], .(Family), numcolwise(sum))
control_families


#'
#' ### Species
#'

## Vertebrate species
control_seq_binary <- control_seq_binary[order(control_seq_binary$AssignmentLevel),]
control_species <- ddply(control_seq_binary[c(36:98),c(1:4,119:120)], .(Tax_assignment), numcolwise(sum))
control_species



#'
#' ### Plots
#' 


## Plot presence of GCN DNA amongst controls after false positive threshold
## applied

p0 <- ggplot(control_reads, aes(x=sample, y=Triturus_cristatus, fill=type)) 
p0 <- p0 + geom_bar(stat="identity", alpha=.5)
p0 <- p0 + scale_y_continuous(limits=c(0,320))
p0 <- p0 + scale_fill_manual(breaks=c("Negative","Positive"),
                               labels=c("Negative (MGW)",
                                        "Positive (cichlid DNA)"),
                               values=c("lightcoral", "cornflowerblue"))
p0 <- p0 + labs(title="Presence of GCN DNA in PCR Controls")
p0 <- p0 + labs(x="PCR Controls", y="GCN reads")
p0 <- p0 + theme(panel.background = element_rect(fill = 'grey98'),
                 plot.title = element_text(face="bold"),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p0 <- p0 + facet_wrap(~ run)
p0



## Plot of negative controls showing number of reads assigned to cichlid

negatives <- subset(control_reads, control_reads$type == "Negative")

p1 <- ggplot(negatives, aes(x=sample, y=Rhamphochromis_esox, fill="lightcoral"))
p1 <- p1 + geom_bar(stat="identity", alpha=.5)
p1 <- p1 + scale_y_continuous(expand = c(0,0), limits = c(0,80000))
p1 <- p1 + labs(title="(a) PCR Negative Controls")
p1 <- p1 + labs(x="Negative Controls (MGW)", y="Cichlid reads")
p1 <- p1 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 legend.position="none",
                 legend.title = element_blank(),
                 text = element_text(size=24))
p1 <- p1 + facet_wrap(~ run)
p1



## Plot contamination of environmental samples with cichlid DNA

## Create colour blind friendly palette:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


p2 <- ggplot(eDNA_reads, aes(x=sample, y=Rhamphochromis_esox, fill=factor(plate)))         
p2 <- p2 + geom_bar(stat="identity", alpha=.5)
p2 <- p2 + scale_y_continuous(expand = c(0,0), limits = c(0,80000))
p2 <- p2 + labs(title="(b) eDNA Samples")
p2 <- p2 + labs(x="eDNA Samples", y="Cichlid reads")
p2 <- p2 + scale_fill_manual(name="2nd PCR Plate",
                             values=cbPalette)
p2 <- p2 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 text = element_text(size=24))
p2 <- p2 + facet_wrap(~ run)
p2



## Show both plots

grid.arrange(arrangeGrob(p1, nrow=1, ncol=1),
             arrangeGrob(p2, nrow=1, ncol=1), nrow=2)



## Plot of positive and negative controls showing number of reads assigned
## to Human

p3 <- ggplot(control_reads, aes(x=sample, y=Homo_sapiens, fill=type))
p3 <- p3 + geom_bar(stat="identity", alpha=.5)
p3 <- p3 + scale_y_continuous(expand = c(0,0), limits = c(0,250000))
p3 <- p3 + scale_fill_manual(breaks=c("Negative","Positive"),
                             labels=c("Negative (MGW)",
                                      "Positive (cichlid DNA)"),
                             values=c("lightcoral", "cornflowerblue"))
p3 <- p3 + labs(title="(a) PCR Controls")
p3 <- p3 + labs(x="Controls", y="Human reads")
p3 <- p3 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 legend.title = element_blank(),
                 text = element_text(size=24))
p3 <- p3 + facet_wrap(~ run)
p3


## Plot contamination of environmental samples with Human DNA

p4 <- ggplot(eDNA_reads, aes(x=sample, y=Homo_sapiens, fill=factor(plate)))         
p4 <- p4 + geom_bar(stat="identity", alpha=.5)
p4 <- p4 + scale_y_continuous(expand = c(0,0), limits = c(0,250000))
p4 <- p4 + labs(title="(b) eDNA Samples")
p4 <- p4 + labs(x="eDNA samples", y="Human reads")
p4 <- p4 + scale_fill_manual(name="2nd PCR Plate",
                             values=cbPalette)
p4 <- p4 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 text = element_text(size=24))
p4 <- p4 + facet_wrap(~ run)
p4



## Show all plots side by side

grid.arrange(arrangeGrob(p3, nrow=1, ncol=1),
             arrangeGrob(p4, nrow=1, ncol=1), nrow=2)



## Plot contamination of controls with DNA of detected species

negatives <- data.frame(subset(control_reads, control_reads$type == "Negative"))
fish <- melt(negatives[,c(1:2,9,12:13,18,21:23,26,31:32,48,53,57,61,65)], 
             id = c("sample","run"))

amphib <- melt(negatives[,c(1:2,15,37:38,51,58,69)], 
               id = c("sample","run"))

bird <- melt(negatives[,c(1:2,8,10,16,19,24,28:30,33:34,40,47,52,54,56,63,66:67)], 
             id = c("sample","run"))

mammal <- melt(negatives[,c(1:2,11,14,17,20,25,27,36,39,41:46,49:50,55,59,62,64,68,70)], 
               id = c("sample","run"))


## Extend colour palette for number of fish species 
colourCount = length(unique(fish$variable))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))


## Plot fish species in negative controls
p5 <- ggplot(fish, aes(x=sample, y=value, fill=variable))
p5 <- p5 + geom_bar(stat="identity")
p5 <- p5 + scale_fill_manual(name = "",
                             values = getPalette(colourCount),
                             labels = c("European eel", "Stone loach",
                                        "Common barbel", "Crucian carp",
                                        "Atlantic herring", "European bullhead",
                                        "Common carp", "Northern pike",
                                        "Three-spined stickleback", "Ruffe",
                                        "Rainbow trout", "Common minnow",
                                        "Ninespine stickleback", "Common roach", 
                                        "European chub"))
p5 <- p5 + scale_y_continuous(expand = c(0,0), limits = c(0,30000))
p5 <- p5 + labs(title="(a) Fish")
p5 <- p5 + labs(x="Negative PCR controls", y="Number of reads")
p5 <- p5 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 text = element_text(size=16))
p5 <- p5 + facet_wrap(~ run)
p5



## Set colour palette for number of amphibian species 
colourCount = length(unique(amphib$variable))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))


## Plot amphibian species in negative controls
p6 <- ggplot(amphib, aes(x=sample, y=value, fill=variable))
p6 <- p6 + geom_bar(stat="identity")
p6 <- p6 + scale_fill_manual(name = "",
                             values = getPalette(colourCount),
                             labels = c("Common toad", "Palmate newt",
                                        "Smooth newt", "Marsh frog",
                                        "Common frog", "Great crested newt"))
p6 <- p6 + scale_y_continuous(expand = c(0,0), limits = c(0,30000))
p6 <- p6 + labs(title="(b) Amphibian")
p6 <- p6 + labs(x="Negative PCR controls", y="Number of reads")
p6 <- p6 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 text = element_text(size=16))
p6 <- p6 + facet_wrap(~ run)
p6



## Extend colour palette for number of bird species 
colourCount = length(unique(bird$variable))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))


## Plot bird species in negative controls
p7 <- ggplot(bird, aes(x=sample, y=value, fill=variable))
p7 <- p7 + geom_bar(stat="identity")
p7 <- p7 + scale_fill_manual(name = "",
                             values = getPalette(colourCount),
                             labels = c("Green-winged teal", "Grey heron",
                                        "Common buzzard", "European goldfinch",
                                        "Great spotted woodpecker", "Eurasian coot",
                                        "Common moorhen", "Eurasian jay",
                                        "Eurasian oystercatcher", "Melodious warbler",
                                        "Domesticated turkey", "Helmeted guineafowl",
                                        "Common pheasant", "Green woodpecker",
                                        "Dunnock", "Eurasian nuthatch",
                                        "Tawny owl", "Common starling"))
p7 <- p7 + scale_y_continuous(expand = c(0,0), limits = c(0,30000))
p7 <- p7 + labs(title="(c) Bird")
p7 <- p7 + labs(x="Negative PCR controls", y="Number of reads")
p7 <- p7 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 text = element_text(size=16))
p7 <- p7 + facet_wrap(~ run)
p7



## Extend colour palette for number of mammal species 
colourCount = length(unique(mammal$variable))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

## Plot mammal species in negative controls
p8 <- ggplot(mammal, aes(x=sample, y=value, fill=variable))
p8 <- p8 + geom_bar(stat="identity")
p8 <- p8 + scale_fill_manual(name = "",
                             values = getPalette(colourCount),
                             labels = c("European water vole", "Cow", "Dog", "Red deer",
                                        "Horse", "Cat", "European hare", "Eurasian otter",
                                        "European badger", "Muntjac deer", "House mouse",
                                        "European polecat", "Bank vole", "Eurasian water shrew",
                                        "European rabbit", "Sheep", "Common pipistrelle",
                                        "Brown rat", "Grey squirrel", "Common shrew",
                                        "Pig", "Fox"),
                             guide = guide_legend(ncol=1))
p8 <- p8 + scale_y_continuous(expand = c(0,0), limits = c(0,30000))
p8 <- p8 + labs(title="(d) Mammal")
p8 <- p8 + labs(x="Negative PCR controls", y="Number of reads")
p8 <- p8 + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill = 'white'),
                 plot.title = element_text(face="bold", hjust=0),
                 axis.ticks = element_blank(), 
                 axis.text.x = element_blank(),
                 text = element_text(size=16))
p8 <- p8 + facet_wrap(~ run)
p8



## Show all plots side by side

grid.arrange(arrangeGrob(p5, nrow=1, ncol=1),
             arrangeGrob(p6, nrow=1, ncol=1),
             arrangeGrob(p7, nrow=1, ncol=1),
             arrangeGrob(p8, nrow=1, ncol=1), nrow=4)



## Are contamination levels significantly different between sequencing
## runs due to addition of mineral oil in PCR reactions?
## Only look at PCR negative controls and eDNA samples as all contamination 
## has been removed from PCR positive controls by application of the 
## species-specific false positive thresholds.

## Model with a binomial distribution as Contamination is a binary 
## (Yes/No, 1/0) response variable
## Model plate as a random factor as PCR is stochastic and interested
## in overall contamination levels between runs and not between plates


## NEGATIVE CONTROLS

contamination_negatives <- subset(control_reads, control_reads$type=="Negative")

contamination_negatives$total <- rowSums(contamination_negatives[,8:70] > 0)
contamination_negatives$contamination <- ifelse(contamination_negatives$total > 0, 1, 0)

m0 <- glmer(contamination ~ (1|plate) + run, 
            family=binomial,
            data = contamination_negatives)

null0 <- glmer(contamination ~ (1|plate) + 1, 
            family=binomial,
            data = contamination_negatives)


anova(null0,m0)


## LRT not significant; when compared to null model, AIC value of model 
## including sequencing run is not too different but p-value is not 
## significant so run model better fit to the data


## Final model:

summary(m0)
anova(m0)
drop1(m0, test="Chi")  # not significant


## Check for overdispersion of final model
## Obtain residual deviance and degrees of freedom
overdisp.glmer(m0)


## Overdispersion test (chi-square)

1-pchisq(45.572, df=111)   # model not overdispersed


## Run overdispersion function

overdisp_fun(m0)     # not overdispersed


## Examine model fit:

## Check model residuals
chkres(m0)

## Plot the fitted data against the observed data
plot(contamination_negatives$contamination ~ fitted(m0))

## Mineral oil did not significantly reduce contamination between sequencing runs.



## eDNA SAMPLES

contamination_eDNA <- eDNA_reads[,c(1:7,40,65)]

contamination_eDNA$contamination <- ifelse(contamination_eDNA$Homo_sapiens > 0, 1,
                                    ifelse(contamination_eDNA$Rhamphochromis_esox > 0, 1,
                                           0))

contamination_eDNA$contamination_human <- ifelse(contamination_eDNA$Homo_sapiens > 0, 1, 0)
contamination_eDNA$contamination_cichlid <- ifelse(contamination_eDNA$Rhamphochromis_esox > 0, 1, 0)


## Human and cichlid contamination

m1 <- glmer(contamination ~ (1|plate) + run, 
            family=binomial,
            data = contamination_eDNA)

null1 <- glmer(contamination ~ (1|plate) + 1, 
               family=binomial,
               data = contamination_eDNA)

anova(null1,m1)

## LRT significant; when compared to null model, AIC value of model including
## sequencing run is lower and p-value significant indicating removal of
## sequencing run has reduced fit and simpler model cannot adequately 
## explain variation in data.

## Final model:

summary(m1)
anova(m1)
drop1(m1, test="Chi")  # significant


## Check for overdispersion of final model
## Obtain residual deviance and degrees of freedom
overdisp.glmer(m1)


## Overdispersion test (chi-square)

1-pchisq(602.673, df=529)   # model overdispersed


## Run overdispersion function

overdisp_fun(m1)     # not overdispersed


## Examine model fit:

## Check model residuals
chkres(m1)

## Plot the fitted data against the observed data
plot(contamination_eDNA$contamination ~ fitted(m1))

## Mineral oil significantly reduced human and cichlid DNA contamination between 
## sequencing runs.



## Human contamination only

m2 <- glmer(contamination_human ~ (1|plate) + run, 
            family=binomial,
            data = contamination_eDNA)

null2 <- glmer(contamination_human ~ (1|plate) + 1, 
               family=binomial,
               data = contamination_eDNA)

anova(null2,m2)

## LRT borderline significant; when compared to null model, AIC value of model 
## including sequencing run is lower, although difference <2 and p-value 
## not significant. This indicates that inclusion of sequencing run does not
## substantially explain any difference in human contamination.


## Final model:

summary(m2)
anova(m2)
drop1(m2, test="Chi")  # not significant


## Check for overdispersion of final model
## Obtain residual deviance and degrees of freedom
overdisp.glmer(m2)


## Overdispersion test (chi-square)

1-pchisq(691.029, df=529)   # model overdispersed


## Run overdispersion function

overdisp_fun(m2)     # not overdispersed


## Examine model fit:

## Check model residuals
chkres(m2)

## Plot the fitted data against the observed data
plot(contamination_eDNA$contamination ~ fitted(m2))

## Mineral oil did not significantly reduce human DNA contamination between 
## sequencing runs.



## Cichlid contamination only

m3 <- glmer(contamination_cichlid ~ (1|plate) + run, 
            family=binomial,
            data = contamination_eDNA)

null3 <- glmer(contamination_cichlid ~ (1|plate) + 1, 
               family=binomial,
               data = contamination_eDNA)

anova(null3,m3)

## LRT significant; when compared to null model, AIC value of model including
## sequencing run is lower and p-value significant indicating removal of
## sequencing run has reduced fit and simpler model cannot adequately 
## explain variation in data.

## Final model:

summary(m3)
anova(m3)
drop1(m3, test="Chi")  # significant


## Check for overdispersion of final model
## Obtain residual deviance and degrees of freedom
overdisp.glmer(m1)


## Overdispersion test (chi-square)

1-pchisq(602.673, df=529)   # model overdispersed


## Run overdispersion function

overdisp_fun(m1)     # not overdispersed


## Examine model fit:

## Check model residuals
chkres(m1)

## Plot the fitted data against the observed data
plot(contamination_eDNA$contamination ~ fitted(m1))

## Mineral oil significantly reduced cichlid DNA contamination between sequencing runs.



#'
#' Now we know that there is some GCN DNA and DNA from other species in the
#' negative controls but only one negative control contains more than 100 reads.
#' 
#' There is some contamination of negatives and eDNA samples with cichlid DNA 
#' and with human DNA. 
#' 
#' Mineral oil did not significantly reduce contamination in negative controls 
#' (all DNA) but did reduce cichlid and cichlid/human contamination of eDNA samples 
#' between sequencing runs, where contamination in run 2 was observedly lower. 
#' 




#' ---
#' 
#' ## 2) eDNA metabarcoding Vs qPCR for GCN detection
#' 
#' 
#' First, I must create a new dataframe containing the number of reads for
#' each taxonomic assignment and metadata on two-step PCR, qPCR and 
#' conventional GCN survey for each sample.
#' 
#' I need to merge the eDNA_only dataframe which has had the species-specific
#' false positive thresholds applied and controls removed with an imported
#' excel sheet containing metadata on each sample.
#' 
#' I will do this using the 'sample' column from the metadata dataframe and
#' the 'sample' column from the reads dataframe.
#' 


## Import metadata for samples
sample_metadata <- read.csv("../Data/GCN_metadata.csv", header=TRUE)


## Create two columns representing results of qPCR with different
## qPCR score thresholds i.e. number of positive qPCR replicates required 
## for a pond to be positive for GCN.
##
## Natural England follows the protocol set out in Biggs et al. (2015),
## where 1 or more positive qPCR replicates equates to a GCN positive pond
## However, 1/12 positive replicates is very low and reliability of this 
## result may be called into question. 
##
## Therefore, we will also analyse the data with a higher threshold of 
## 4 or more positive qPCR replicates to see how results of qPCR differ
## and how qPCR and metabarcoding detection differ.
##
## 4/12 qPCR replicates (1/3) was chosen as the threshold as this is more 
## comparable to metabarcoding where we perform three PCR1 replicates, and
## pool PCR products. Only one replicate may amplify in some cases.
## 


## Create qPCR_1rep column where samples with qPCR score >= 1 are coded as 
## P and samples where qPCR score = 0 are coded as N
sample_metadata$qPCR_thresh1 <- ifelse(sample_metadata$qPCR_score > 0, "P", "N")


## Create qPCR_4rep column where samples with qPCR score >= 4 are coded as 
## P and samples where qPCR score <= 3 are coded as N
sample_metadata$qPCR_thresh4 <- ifelse(sample_metadata$qPCR_score >= 4, "P", "N")


## Now, the sample columns from each data frame need to match type in order
## to perform the merge
## eDNA_only column is a character while metadata column is a factor
sample_metadata$sample <- as.character(sample_metadata$sample)

meta_reads <- merge(sample_metadata, eDNA_only, all=TRUE)

summary(meta_reads)
head(meta_reads)
names(meta_reads)
str(meta_reads)



#'
#' Next, I need to create several new columns to compare methods. Based on
#' the number of GCN reads after the false positive threshold was applied,
#' I need to know whether samples are positive or negative for GCN: mbTA_GCN
#' 
#' From these reads, I also need to know the proportion of GCN reads in a 
#' sample, expressed as a whole integer: GCNprop
#' 

## Create mbTA_GCN column where samples with GCN reads are coded as P and 
## samples without GCN reads are coded as N
meta_reads$mbTA_GCN <- ifelse(meta_reads$Triturus_cristatus == 0, "N", "P")

## Check this has worked
head(meta_reads[,70:80])


## The data frame at this point will be required for analysis in Chapter 2
## Export this data frame as csv file
write.csv(meta_reads, "../../Chapter 2 - GCN interactions/Data/MetabarcodingData.csv", 
          row.names=FALSE)


## Now create GCNprop
meta_reads$GCNprop <- meta_reads$Triturus_cristatus/meta_reads$reads_total * 100

## Check this has worked
head(meta_reads[,75:81])

## Round to nearest whole number
meta_reads$GCNprop <- round(meta_reads$GCNprop)
head(meta_reads[,75:81])

## Replace NAs with 0
meta_reads$GCNprop[is.na(meta_reads$GCNprop)] <- 0



#'
#'
#' ## Summary Statistics
#' 
#' Plot the difference in GCN detection between eDNA approaches.
#' 
#' Perform kappa coefficient, a measure of strength of agreement, and 
#' chi-square test, comparison of means, to compare qPCR and eDNA 
#' metabarcoding (with and without threshold applied) for detection of GCN.
#' 
#' First, create a dataframe containing the metabarcoding data, qPCR data
#' and conventional data.
#' 

## Create dataframe
GCN_methods <- data.frame(meta_reads$sample,
                          meta_reads$conv_P.N,
                          meta_reads$qPCR_thresh1,
                          meta_reads$qPCR_thresh4,
                          meta_reads$mbNT_GCN,
                          meta_reads$mbTA_GCN)

colnames(GCN_methods) <- c("sample","EggSearch","qPCR_NT","qPCR_TA","mbNT","mbTA")


## Calculate proportion of negative and positive detections by each method

trad_table <- table(GCN_methods$EggSearch)
prop <- prop.table(trad_table) 
trad_prop <- data.frame(prop)

qPCR_NT_table <- table(GCN_methods$qPCR_NT)
prop <- prop.table(qPCR_NT_table) 
qPCR_NT_prop <- data.frame(prop)

qPCR_TA_table <- table(GCN_methods$qPCR_TA)
prop <- prop.table(qPCR_TA_table) 
qPCR_TA_prop <- data.frame(prop)

mbNT_table <- table(GCN_methods$mbNT)
prop <- prop.table(mbNT_table) 
mbNT_prop <- data.frame(prop)

mbTA_table <- table(GCN_methods$mbTA)
prop <- prop.table(mbTA_table) 
mbTA_prop <- data.frame(prop)


## Create new dataframe of proportions
method_prop <- cbind(trad_prop, qPCR_NT_prop[,2], qPCR_TA_prop[,2], mbNT_prop[,2], mbTA_prop[,2])
colnames(method_prop) <- c("GCN","Egg search","qPCR NT","qPCR TA","Metabarcoding NT","Metabarcoding TA")


## Use melt in reshape package to show positive and negative proportion for 
## each method and enable plotting proportion by method and coloured by GCN
## detection
method_prop <- melt(method_prop, id.vars='GCN')
colnames(method_prop) <- c("GCN","Method","Proportion")

## Add frequencies of positive and negative ponds by each method in new 
## column
method_prop$Frequency <- c(448,58,267,265,358,174,350,182,383,149)


#'
#' Now, plot GCN detection by each method.
#' 

## Plot of proportion of ponds GCN present and absent in as determined by
## metabarcoding and qPCR
p9 <- ggplot(method_prop[order(method_prop$GCN,decreasing=T),],
             aes(x=Method, y=Proportion, fill=GCN))
p9 <- p9 + geom_bar(stat="identity", colour="black",  alpha=.7)
p9 <- p9 + geom_text(aes(label=Frequency), size=10, position=position_fill(), vjust=1.5)
p9 <- p9 + scale_y_continuous(labels = scales::percent, limits = c(0, 1.05),
                              expand = c(0, 0))
p9 <- p9 + labs(x="Method", y="Ponds surveyed")
p9 <- p9 + scale_fill_manual(name=expression(italic("T. cristatus")),
                             values=c("grey45","orange"),
                             labels=c("Negative",
                                      "Positive"))
p9 <- p9 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold"),
                 text = element_text(size=28),
                 legend.key.size = unit(2, 'lines'))
p9



## Now, examine which methods agree and disagree with one another for
## individual samples

## trad vs qPCR NT
GCN_methods$ADtradqpcrNT <- factor(ifelse(GCN_methods$Trad==GCN_methods$qPCR_NT,"Agree","Disagree"))
summary(GCN_methods$ADtradqpcrNT)

GCN_methods$ExptradqpcrNT <- factor(ifelse(GCN_methods$Trad=="N"&GCN_methods$qPCR_NT=="N","Negative",
                                  ifelse(GCN_methods$Trad=="P"&GCN_methods$qPCR_NT=="P","Positive",
                                  ifelse(GCN_methods$Trad=="P"&GCN_methods$qPCR_NT=="N","Traditional",
                                  ifelse(GCN_methods$Trad=="N"&GCN_methods$qPCR_NT=="P","qPCR NT",
                                         "Fail")))))
summary(GCN_methods$ExptradqpcrNT)


## trad vs qPCR TA
GCN_methods$ADtradqpcrTA <- factor(ifelse(GCN_methods$Trad==GCN_methods$qPCR_TA,"Agree","Disagree"))
summary(GCN_methods$ADtradqpcrTA)

GCN_methods$ExptradqpcrTA <- factor(ifelse(GCN_methods$Trad=="N"&GCN_methods$qPCR_TA=="N","Negative",
                                         ifelse(GCN_methods$Trad=="P"&GCN_methods$qPCR_TA=="P","Positive",
                                                ifelse(GCN_methods$Trad=="P"&GCN_methods$qPCR_TA=="N","Traditional",
                                                       ifelse(GCN_methods$Trad=="N"&GCN_methods$qPCR_TA=="P","qPCR TA",
                                                              "Fail")))))
summary(GCN_methods$ExptradqpcrTA)



## trad vs metabarcoding NT
GCN_methods$ADtradmbNT <- factor(ifelse(GCN_methods$Trad==GCN_methods$mbNT,"Agree","Disagree"))
summary(GCN_methods$ADtradmbNT)

GCN_methods$ExptradmbNT <- factor(ifelse(GCN_methods$Trad=="N"&GCN_methods$mbNT=="N","Negative",
                                  ifelse(GCN_methods$Trad=="P"&GCN_methods$mbNT=="P","Positive",
                                  ifelse(GCN_methods$Trad=="P"&GCN_methods$mbNT=="N","Traditional",
                                  ifelse(GCN_methods$Trad=="N"&GCN_methods$mbNT=="P","metabarcoding NT",
                                          "Fail")))))
summary(GCN_methods$ExptradmbNT)


## trad vs metabarcoding TA
GCN_methods$ADtradmbTA <- factor(ifelse(GCN_methods$Trad==GCN_methods$mbTA,"Agree","Disagree"))
summary(GCN_methods$ADtradmbTA)

GCN_methods$ExptradmbTA <- factor(ifelse(GCN_methods$Trad=="N"&GCN_methods$mbTA=="N","Negative",
                                         ifelse(GCN_methods$Trad=="P"&GCN_methods$mbTA=="P","Positive",
                                         ifelse(GCN_methods$Trad=="P"&GCN_methods$mbTA=="N","Traditional",
                                         ifelse(GCN_methods$Trad=="N"&GCN_methods$mbTA=="P","metabarcoding TA",
                                              "Fail")))))
summary(GCN_methods$ExptradmbTA)


## qPCR NT vs metabarcoding NT
GCN_methods$ADqPCRNTmbNT <- factor(ifelse(GCN_methods$qPCR_NT==GCN_methods$mbNT,"Agree","Disagree"))
summary(GCN_methods$ADqPCRNTmbNT)

GCN_methods$ExpqPCRNTmbNT <- factor(ifelse(GCN_methods$qPCR_NT=="N"&GCN_methods$mbNT=="N","Negative",
                                  ifelse(GCN_methods$qPCR_NT=="P"&GCN_methods$mbNT=="P","Positive",
                                  ifelse(GCN_methods$qPCR_NT=="P"&GCN_methods$mbNT=="N","qPCR NT",
                                  ifelse(GCN_methods$qPCR_NT=="N"&GCN_methods$mbNT=="P","metabarcoding NT",
                                          "Fail")))))
summary(GCN_methods$ExpqPCRNTmbNT)


## qPCR TA vs metabarcoding NT
GCN_methods$ADqPCRTAmbNT <- factor(ifelse(GCN_methods$qPCR_TA==GCN_methods$mbNT,"Agree","Disagree"))
summary(GCN_methods$ADqPCRTAmbNT)

GCN_methods$ExpqPCRTAmbNT <- factor(ifelse(GCN_methods$qPCR_TA=="N"&GCN_methods$mbNT=="N","Negative",
                                           ifelse(GCN_methods$qPCR_TA=="P"&GCN_methods$mbNT=="P","Positive",
                                                  ifelse(GCN_methods$qPCR_TA=="P"&GCN_methods$mbNT=="N","qPCR TA",
                                                         ifelse(GCN_methods$qPCR_TA=="N"&GCN_methods$mbNT=="P","metabarcoding NT",
                                                                "Fail")))))
summary(GCN_methods$ExpqPCRTAmbNT)


## qPCR NT vs metabarcoding TA
GCN_methods$ADqPCRNTmbTA <- factor(ifelse(GCN_methods$qPCR_NT==GCN_methods$mbTA,"Agree","Disagree"))
summary(GCN_methods$ADqPCRNTmbTA)

GCN_methods$ExpqPCRNTmbTA <- factor(ifelse(GCN_methods$qPCR_NT=="N"&GCN_methods$mbTA=="N","Negative",
                                  ifelse(GCN_methods$qPCR_NT=="P"&GCN_methods$mbTA=="P","Positive",
                                  ifelse(GCN_methods$qPCR_NT=="P"&GCN_methods$mbTA=="N","qPCR NT",
                                  ifelse(GCN_methods$qPCR_NT=="N"&GCN_methods$mbTA=="P","metabarcoding TA",
                                          "Fail")))))
summary(GCN_methods$ExpqPCRNTmbTA)


## qPCR TA vs metabarcoding TA
GCN_methods$ADqPCRTAmbTA <- factor(ifelse(GCN_methods$qPCR_TA==GCN_methods$mbTA,"Agree","Disagree"))
summary(GCN_methods$ADqPCRTAmbTA)

GCN_methods$ExpqPCRTAmbTA <- factor(ifelse(GCN_methods$qPCR_TA=="N"&GCN_methods$mbTA=="N","Negative",
                                           ifelse(GCN_methods$qPCR_TA=="P"&GCN_methods$mbTA=="P","Positive",
                                                  ifelse(GCN_methods$qPCR_TA=="P"&GCN_methods$mbTA=="N","qPCR TA",
                                                         ifelse(GCN_methods$qPCR_TA=="N"&GCN_methods$mbTA=="P","metabarcoding TA",
                                                                "Fail")))))
summary(GCN_methods$ExpqPCRTAmbTA)


## qPCR NT VS qPCR TA 
GCN_methods$ADqPCRNTqPCRTA <- factor(ifelse(GCN_methods$qPCR_NT==GCN_methods$qPCR_TA,"Agree","Disagree"))
summary(GCN_methods$ADqPCRNTqPCRTA)

GCN_methods$ExpqPCRNTqPCRTA <- factor(ifelse(GCN_methods$qPCR_NT=="N"&GCN_methods$qPCR_TA=="N","Negative",
                                         ifelse(GCN_methods$qPCR_NT=="P"&GCN_methods$qPCR_TA=="P","Positive",
                                                ifelse(GCN_methods$qPCR_NT=="P"&GCN_methods$qPCR_TA=="N","qPCR NT",
                                                       ifelse(GCN_methods$qPCR_NT=="N"&GCN_methods$qPCR_TA=="P","qPCR TA",
                                                              "Fail")))))
summary(GCN_methods$ExpqPCRNTqPCRTA)


## metabarcoding NT VS metabarcoding TA 
GCN_methods$ADmbNTmbTA <- factor(ifelse(GCN_methods$mbNT==GCN_methods$mbTA,"Agree","Disagree"))
summary(GCN_methods$ADmbNTmbTA)

GCN_methods$ExpmbNTmbTA <- factor(ifelse(GCN_methods$mbNT=="N"&GCN_methods$mbTA=="N","Negative",
                                 ifelse(GCN_methods$mbNT=="P"&GCN_methods$mbTA=="P","Positive",
                                 ifelse(GCN_methods$mbNT=="P"&GCN_methods$mbTA=="N","metabarcoding (NT)",
                                 ifelse(GCN_methods$mbNT=="N"&GCN_methods$mbTA=="P","metabarcoding (TA)",
                                        "Fail")))))
summary(GCN_methods$ExpmbNTmbTA)



#'
#' Create a venn diagram which examines overlap between different monitoring
#' approaches for great crested newt.
#' 

source("http://www.bioconductor.org/biocLite.R")
class(biocLite)
biocLite("limma")
library("limma")

Trad <- (GCN_methods$Trad == "P")
qPCR_NT <- (GCN_methods$qPCR_NT == "P")
qPCR_TA <- (GCN_methods$qPCR_TA == "P")
mb_NT <- (GCN_methods$mbNT == "P")
mb_TA <- (GCN_methods$mbTA == "P")

venn.df <- cbind(Trad, qPCR_NT, qPCR_TA, mb_NT, mb_TA)
a <- vennCounts(venn.df)
a

svg(filename="Venn.svg", width = 13, height = 10)
vennDiagram(a, names = c("Conventional", "qPCR NT", "qPCR TA", 
                         "Metabarcoding NT", "Metabarcoding TA"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("orange","purple","limegreen","deepskyblue","gold"))
dev.off()



#'
#' Next, examine strength of agreement between the eDNA approaches. Compare
#' only the metabarcoding detection with the GCN false positive threshold
#' applied to the qPCR detection.
#' 
#' Use kappa coefficient, a test which has been performed in several eDNA 
#' papers comparing eDNA to conventional survey methods.
#' 


## Metabarcoding NT and qPCR NT "agree" on detection 161 times
## Metabarcoding NT and qPCR NT "agree" on no detection 246 times
## Metabarcoding NT finds GCN 21 times that qPCR NT doesn't
## qPCR NT finds GCN 104 times that metabarcoding NT doesn't
## Neither method completely fails

## Set paramaters
a <- 161
b <- 21
c <- 104
d <- 246
N <- 532

## Run kappa coefficient
method.k <- kappa.function(method)


## Summarise output of kappa function

Output <- c("P(a)","P(e)","kappa", "P(pos)", "P(neg)", "P(index)", "B(index)", "PABAK")
Values <- method.k
kappa <- data.frame(Output, Values)
print(kappa)

#'
#' The kappa coefficient is 0.53, indicating moderate agreement 
#' between eDNA metabarcoding NT and qPCR NT for detection of GCN.
#' 


## Metabarcoding TA and qPCR NT "agree" on detection 138 times
## Metabarcoding TA and qPCR NT "agree" on no detection 256 times
## Metabarcoding TA finds GCN 11 times that qPCR NT doesn't
## qPCR NT finds GCN 127 times that metabarcoding TA doesn't
## Neither method completely fails

## Set paramaters
a <- 138
b <- 11
c <- 127
d <- 256
N <- 532


## Run kappa coefficient
method.k <- kappa.function(method)


## Summarise output of kappa function

Output <- c("P(a)","P(e)","kappa", "P(pos)", "P(neg)", "P(index)", "B(index)", "PABAK")
Values <- method.k
kappa <- data.frame(Output, Values)
print(kappa)

#'
#' The kappa coefficient is 0.48, indicating moderate agreement 
#' between eDNA metabarcoding TA and qPCR NT for detection of GCN.
#' 


## Metabarcoding NT and qPCR TA "agree" on detection 134 times
## Metabarcoding NT and qPCR TA "agree" on no detection 310 times
## Metabarcoding NT finds GCN 48 times that qPCR TA doesn't
## qPCR TA finds GCN 40 times that metabarcoding NT doesn't
## Neither method completely fails

## Set paramaters
a <- 134
b <- 48
c <- 40
d <- 310
N <- 532


## Run kappa coefficient
method.k <- kappa.function(method)


## Summarise output of kappa function

Output <- c("P(a)","P(e)","kappa", "P(pos)", "P(neg)", "P(index)", "B(index)", "PABAK")
Values <- method.k
kappa <- data.frame(Output, Values)
print(kappa)

#'
#' The kappa coefficient is 0.63, indicating good agreement 
#' between eDNA metabarcoding NT and qPCR TA for detection of GCN.
#' 


## Metabarcoding TA and qPCR TA "agree" on detection 123 times
## Metabarcoding TA and qPCR TA "agree" on no detection 332 times
## Metabarcoding TA finds GCN 26 times that qPCR TA doesn't
## qPCR TA finds GCN 51 times that metabarcoding TA doesn't
## Neither method completely fails

## Set paramaters
a <- 123
b <- 26
c <- 51
d <- 332
N <- 532


## Run kappa coefficient
method.k <- kappa.function(method)


## Summarise output of kappa function

Output <- c("P(a)","P(e)","kappa", "P(pos)", "P(neg)", "P(index)", "B(index)", "PABAK")
Values <- method.k
kappa <- data.frame(Output, Values)
print(kappa)


#'
#' The kappa coefficient is 0.66, indicating good agreement 
#' between eDNA metabarcoding TA and qPCR TA for detection of GCN.
#' 


#'
#' Is there a statistically significant difference in GCN detection between 
#' the qPCR (NT and TA) and metabarcoding (NT and TA)?
#' 
#' Null hypothesis: There is no difference in detection.
#' Alternate hypothesis: There is a difference in detection.
#' 

chi.data <- matrix(method_prop[c(3:4,7:8),4], 2, 2, byrow=TRUE,
                   dimnames=list(c("qPCR NT","Metabarcoding NT"), c("GCN Negative", "GCN Positive")))

chisq.test(chi.data)


chi.data <- matrix(method_prop[c(3:4,9:10),4], 2, 2, byrow=TRUE,
                   dimnames=list(c("qPCR NT","Metabarcoding TA"), c("GCN Negative", "GCN Positive")))

chisq.test(chi.data)


chi.data <- matrix(method_prop[c(5:6,7:8),4], 2, 2, byrow=TRUE,
                   dimnames=list(c("qPCR TA","Metabarcoding NT"), c("GCN Negative", "GCN Positive")))

chisq.test(chi.data)


chi.data <- matrix(method_prop[c(5:6,9:10),4], 2, 2, byrow=TRUE,
                   dimnames=list(c("qPCR TA","Metabarcoding TA"), c("GCN Negative", "GCN Positive")))

chisq.test(chi.data)


#'
#' Significant difference in detection which indicates that qPCR is more
#' sensitive for detection of GCN, except where threshold applied to qPCR
#' and metabarcoding with/without threshold.
#'
#' However, metabarcoding TA detected GCN in eDNA samples which were negative
#' by qPCR. Is there more disagreement between the two methods at lower
#' qPCR scores?
#'   


## Calculate proportion of disagreement for qPCR NT positive and negative 
## samples with metabarcoding NT

## Number of negative metabarcoding NT detections in qPCR NT positive ponds
qPCR_NT_pos_disagree <- 104/265

## Number of positive metabarcoding NT detections in qPCR NT negative ponds
qPCR_NT_neg_disagree <- 21/267


## Calculate proportion of disagreement for qPCR TA positive and negative 
## samples with metabarcoding TA

## Number of negative metabarcoding TA detections in qPCR TA positive ponds
qPCR_TA_pos_disagree <- 51/174

## Number of positive metabarcoding TA detections in qPCR TA negative ponds
qPCR_TA_neg_disagree <- 26/358




#' ---
#' 
#' ## 3) qPCR score Vs number of GCN reads
#' 
#' Visually and statistically examine whether there is a relationship between
#' qPCR score (the number of positive replicates for GCN) and the number of
#' GCN reads.
#' 
#' The average number of reads for a given qPCR score will be used instead
#' of actual number of reads as samples were only sequenced once by 
#' metabarcoding but 12 replicates were performed for each sample by qPCR.
#' Furthermore, the number of samples belonging to each qPCR score is uneven.
#' 
#' For a sample to be classed as positive by qPCR for GCN, 3/12 replicates
#' must be positive for the species.
#' 
#' First, calculate the average number of reads for each qPCR score and 
#' collate this data in a new dataframe.
#' 

GCN_mbNT <- vert_seq[96,]

## prepare the assignment names and transpose the dataframe
rownames(GCN_mbNT) <- GCN_mbNT[,1]
GCN_mbNT <- GCN_mbNT[,2:761]
GCN_mbNT <- data.frame(t(GCN_mbNT))


## move the sample names into the first column
GCN_mbNT$sample <- rownames(GCN_mbNT)
GCN_mbNT <- GCN_mbNT[,c(length(GCN_mbNT),1:(length(GCN_mbNT)-1))]
row.names(GCN_mbNT) <- NULL
colnames(GCN_mbNT) <- c("sample", "Triturus_cristatus_NT")
GCN_mbNT <- GCN_mbNT[-c(509:622,647:760),]


## Now, the sample columns from each data frame need to match type in order
## to perform the merge
## eDNA_only column is a character while metadata column is a factor
mb_qPCR_df <- merge(meta_reads, GCN_mbNT, all=TRUE)


## Create new dataframe
new.df <- data.frame(mb_qPCR_df$sample,
                     mb_qPCR_df$qPCR_score,
                     mb_qPCR_df$qPCR_thresh1,
                     mb_qPCR_df$qPCR_thresh4,
                     mb_qPCR_df$Triturus_cristatus_NT,
                     mb_qPCR_df$mbNT_GCN,
                     mb_qPCR_df$Triturus_cristatus,
                     mb_qPCR_df$mbTA_GCN)

colnames(new.df) <- c("sample","qPCR_score","qPCR_NT","qPCR_TA","GCN_reads_NT","GCN_mbNT","GCN_reads_TA","GCN_mbTA")

str(new.df)


## Convert qPCR score to a factor
new.df$qPCR_score <- as.factor(new.df$qPCR_score)
str(new.df)

## Calculate maximum great crested newt read count in qPCR negative samples
## and qPCR positive samples
qPCR_NT_neg <- subset(new.df, new.df$qPCR_NT == "N")
qPCR_NT_pos <- subset(new.df, new.df$qPCR_NT == "P")

max(qPCR_NT_neg$GCN_reads_NT)
max(qPCR_NT_pos$GCN_reads_NT)
max(qPCR_NT_neg$GCN_reads_TA)
max(qPCR_NT_pos$GCN_reads_TA)


qPCR_TA_neg <- subset(new.df, new.df$qPCR_TA == "N")
qPCR_TA_pos <- subset(new.df, new.df$qPCR_TA == "P")

max(qPCR_TA_neg$GCN_reads_NT)
max(qPCR_TA_pos$GCN_reads_NT)
max(qPCR_TA_neg$GCN_reads_TA)
max(qPCR_TA_pos$GCN_reads_TA)



## Calculate average GCN reads for each qPCR score
averages_NT <- tapply(new.df$GCN_reads_NT, new.df$qPCR_score, mean)
averages_NT

averages_TA <- tapply(new.df$GCN_reads_TA, new.df$qPCR_score, mean)
averages_TA


## Calculate standard deviation for each qPCR score
se_NT <- tapply(new.df$GCN_reads_NT, new.df$qPCR_score, se)
se_NT

se_TA <- tapply(new.df$GCN_reads_TA, new.df$qPCR_score, se)
se_TA


## Create new dataframe containing only averages, sd, qPCR score, 
## qPCR detection

average.df <- data.frame(qPCR_score = seq(from=0, to=12, by=1),
                         qPCR_NT = "",
                         qPCR_TA = "",
                         average_reads_NT = averages_NT,
                         average_reads_TA = averages_TA,
                         se_NT = se_NT,
                         se_TA = se_TA)


average.df$qPCR_NT <- factor(ifelse(average.df$qPCR_score == "0", "N", "P"))

average.df$qPCR_TA <- factor(ifelse(average.df$qPCR_score == "0", "N",
                             ifelse(average.df$qPCR_score == "1", "N", 
                             ifelse(average.df$qPCR_score == "2", "N",
                             ifelse(average.df$qPCR_score == "3", "N",
                                    "P")))))

head(average.df)


## Plot qPCR score compared to average number of GCN reads obtained for 
## each qPCR score

## Define errorbars
limits_NT <- aes(ymax = average_reads_NT + se_NT, ymin = average_reads_NT - se_NT)
limits_TA <- aes(ymax = average_reads_TA + se_TA, ymin = average_reads_TA - se_TA)

## Because the bars and errorbars have different widths, we need to specify 
## how wide the objects we are dodging are
dodge <- position_dodge(width=0.9)


## Calculate slope and line of best fit
coef(lm(average_reads_NT ~ qPCR_score, data = average.df))
coef(lm(average_reads_TA ~ qPCR_score, data = average.df))


## Plot relationship

p10 <- ggplot(average.df, aes(x=factor(qPCR_score), y=average_reads_TA))
p10 <- p10 + geom_point(size=2) 
p10 <- p10 + geom_errorbar(limits_NT, position=dodge, width=0.25)
p10 <- p10 + geom_abline(intercept = -751.268, slope = 374.158, size = 0.7, colour = "red")
p10 <- p10 + labs(y="Average number of crested newt reads per qPCR score", x = "qPCR score") 
p10 <- p10 + scale_y_continuous(expand = c(0,0), limits = c(0,12500), breaks = seq(0, 12000, by = 2000)) 
p10 <- p10 + theme(panel.background = element_rect(fill = 'white'),
                   plot.title = element_text(hjust = 0),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   legend.position = "none",
                   text = element_text(size=20))
p10




#'
#' Test whether the correlation between number of GCN reads and qPCR score 
#' is significant.
#' 


## spearman correlation (non-normal distributions):

## Without metabarcoding threshold, qPCR score 0-12
cor.test(average.df$qPCR_score, average.df$average_reads_NT, method="spearman")
r_squared <- cor(average.df$qPCR_score, average.df$average_reads_NT, method="spearman")^2 
r_squared


## With metabarcoding threshold, qPCR score 0-12
cor.test(average.df$qPCR_score, average.df$average_reads_TA, method="spearman")
r_squared <- cor(average.df$qPCR_score, average.df$average_reads_TA, method="spearman")^2 
r_squared



#'
#' The non-parametric test indicates a significant correlation
#' between average number of GCN reads and qPCR score (*P* < 0.05), 
#' regardless of whether metabarcoding threshold are applied.
#' 
#' However, PCR and sequencing are stochastic processes and variation in 
#' amplification/sequencing success can be introduced by multiple routes.
#' Therefore, run and PCR plate may absorb some of the variation in production
#' of GCN reads as samples were amplified on different PCR machines and 
#' sequenced at different times.
#' 
#' A GLMM would allow run and plate to be treated as random effects, and 
#' discount the variation caused by these variables in GCN read production.
#' 
#' Use a GLMM to explore the relationship between GCN read production and
#' qPCR score. Use percentage of GCN reads within an eDNA sample rather than 
#' actual number. Percentage puts the number of GCN reads within context of 
#' the total number of reads obtained for a sample.
#' 
#' In addition to qPCR score, sample inhibition and degradation, and post-PCR
#' eDNA concentration may affect GCN read production.
#' 

## GLMM
m1 <- glmer(GCNprop ~ (1|run) + (1|plate) 
            + qPCR_score + qPCR_degradation + qPCR_inhibition 
            + concentration, 
            family = poisson, 
            data = meta_reads) 


## Test model for overdispersion
## Obtain residual deviance and degrees of freedom
overdisp.glmer(m1)


## Overdispersion test (chi-square)

1-pchisq(3853.567, df=525)   # model overdispersed


## Must switch to quasi-Poisson or negative binomial GLMM to resolve 
## overdispersion

## Rob Thomas ' Data Analysis with R' suggests that for overdispersion
## statistics between 2-15, use a quasi-Poisson model


## Penalised quasi-likelihood model
## Similar to poisson distributed model but designed to handle
## overdispersed data

m1.qp <- glmmPQL(GCNprop ~ qPCR_score + qPCR_degradation 
                 + qPCR_inhibition + concentration,
                 random = list(~1|run, ~1|plate),
                 family = quasipoisson(link="log"), 
                 data = meta_reads)

summary(m1.qp)
Anova(m1.qp, type = "III")


## Check for overdispersion of quasi-likelihood model
## Run custom function to statisticall test for model overdispersion using
## Pearson residuals
overdisp_fun(m1.qp)

## Can't perform function


## Evaluate model fit:
chkres.PQL(m1.qp)   
plot(meta_reads$GCNprop ~ fitted(m1.qp))

## model still overdispersed based on residual plots and plot of 
## predicted data from model against the actual data



## Try a zero inflated model to allow for overdispersion in data due
## to occurrence of 0 reads and 0 ng/ul concentration in majority of 
## samples 

## Install packages required to fit zero inflated models

install.packages("R2admb")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")

library(glmmADMB)


## Convert random factors to factors in data frame for ZIP model
meta_reads$f.run <- as.factor(meta_reads$run)
meta_reads$f.plate <- as.factor(meta_reads$plate)


## Zero inflated poisson model

m1.zip <- glmmadmb(GCNprop ~ (1|f.run) + (1|f.plate) 
                   + qPCR_score + concentration 
                   + qPCR_degradation + qPCR_inhibition, 
                   data = meta_reads,
                   zeroInflation = TRUE,
                   family = "poisson")

summary(m1.zip)

## Check for overdispersion of zero-inflated Poisson model
## Run custom function to statisticall test for model overdispersion using
## Pearson residuals
overdisp_fun(m1.zip)

## model still overdispersed


## Evaluate model fit:
chkres.PQL(m1.zip)   
plot(meta_reads$GCNprop ~ fitted(m1.zip)) 


## Better fit than quasi-likelihood model but still overdispersed



## Try Zero inflated negative binomial model

m1.zinb <- glmmadmb(GCNprop ~ (1|f.run) + (1|f.plate) 
                   + qPCR_score + concentration 
                   + qPCR_degradation, 
                   data = meta_reads,
                   zeroInflation = TRUE,
                   family = "nbinom")

summary(m1.zinb)


## COULD NOT FIT MODEL


## Try a negative binomial distribution to resolve overdispersion

m1.nb <- glmer.nb(GCNprop ~ (1|run) + (1|plate) 
                  + qPCR_score + qPCR_degradation 
                  + qPCR_inhibition + concentration, 
                  data = meta_reads) 

summary(m1.nb)


## Check for overdispersion of negative binomial model
## Obtain residual deviance and degrees of freedom
overdisp.glmer(m1.nb)

## Overdispersion test (chi-square)

1-pchisq(281.126, df=525)   # model not overdispersed



## Use drop1() to determine which variable should be removed next based on
## AIC value. Values with lower AIC explain less variation in response 
## variable

drop1(m1.nb, test = "Chi") 


## Drop qPCR inhibition

m2.nb <- glmer.nb(GCNprop ~ (1|run) + (1|plate) 
                  + qPCR_score + qPCR_degradation
                  + concentration, 
                  data = meta_reads) 

anova(m1.nb, m2.nb, test = "Chi")


## p-value of the LRT not significant and AIC value of m2.nb lower
## so removal of qPCR inhibition has improved model fit.

drop1(m2.nb, test = "Chi")


## Drop qPCR degradation

m3.nb <- glmer.nb(GCNprop ~ (1|run) + (1|plate) 
                  + qPCR_score + concentration, 
                  data = meta_reads) 

anova(m2.nb, m3.nb, test = "Chi")


## p-value of the LRT not significant although AIC value of m3.nb not that
## different so removal of qPCR degradation has improved model fit.

drop1(m3.nb, test = "Chi")



## Remove concentration to see if final model has been achieved.

m4.nb <- glmer.nb(GCNprop ~ (1|run) + (1|plate) + qPCR_score, 
                  data = meta_reads) 

anova(m3.nb, m4.nb, test = "Chi")

## p-value of the LRT highly significant and AIC value of m4.nb greatly
## increased so removal of concentration has worsened model fit.



## Compare m3.nb to null model to ensure this is best possible fit to the 
## data

m5.nb <- glmer.nb(GCNprop ~ (1|run) + (1|plate), 
                  data = meta_reads)

anova(m3.nb, m5.nb, test = "Chi")

## same as above, final model has been reached



## Final model:

summary(m3.nb)
anova(m3.nb)
drop1(m3.nb, test = "Chi")



## Test final model for overdispersion
## Obtain residual deviance and degrees of freedom
overdisp.glmer(m3.nb)


## Overdispersion test (chi-square)

1-pchisq(280.646, df=527)   # model not overdispersed


## Run overdispersion function

overdisp_fun(m3.nb)    # model supposedly overdispersed


## Examine model fit:

## Check model residuals
chkres(m3.nb)

## Plot the fitted data against the observed data
plot(meta_reads$GCNprop ~ fitted(m3.nb))

## Best fit to data so far


## An additional test of model fit is the Hosmer and Lemeshow Goodness 
## of Fit Test, which is similar to a chi-square test and indicates the 
## extent to which the model provides better fit than a null model with 
## no predictors.
## If chi-square goodness of fit is not significant, then the model has 
## adequate fit. Similarly, if the test is significant, the model does not 
## adequately fit the data

hoslem.test(meta_reads$GCNprop, fitted(m3.nb))

# Model fits well as p-value is not significant 
# i.e. no significant difference between the model and the observed data




#' 
#' Many metabarcoding studies have been unable to find any relationship 
#' between eDNA concentration and the number of reads for individual species
#' within communities. 
#' 
#' Here, this does not hold true for GCN. We actually find a negative 
#' relationship. This is likely due to the presence of 0s in our dataset.
#' However, it may be suggestive of DNA competition in community samples,
#' where GCN DNA is being outcompeted for PCR amplification by DNA from 
#' other vertebrate species. 
#' 


## PLOT MODEL FIT:

## Obtain predicted values for the full data set: method.dat
## se.fit = TRUE will obtain the standard error for each of these predictions

fit <- data.frame(predictSE(m3.nb, newdata=meta_reads, 
                            se.fit = TRUE, print.matrix=T))


## Create new data set with fitted values and original data for plots 1 
## and 2
box.dat <- cbind(meta_reads, fit) 


## PLOT 1: qPCR score and proportion GCN reads 

p11 <- ggplot(box.dat) + ggtitle('(a)')
p11 <- p11 + geom_jitter(aes(x=factor(qPCR_score), y=GCNprop), shape=1, show.legend=FALSE, cex=1)
p11 <- p11 + geom_boxplot(aes(x=factor(qPCR_score), y=fit),alpha=0.7, show.legend=FALSE)
p11 <- p11 + coord_cartesian(ylim=c(0,100))
p11 <- p11 + labs(x = "qPCR score",
                  y = expression(paste("Proportion of ", italic("T. cristatus"), " reads in total reads per sample (%)")))
p11 <- p11 + theme(panel.background = element_rect(fill = 'white'),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=24))
p11


## PLOT 2: eDNA concentration and proportion GCN reads

## Create a range of concentration values which increase by 
## to predict GCN read proportion as concentration increases
range <- seq(from=min(meta_reads$concentration), 
             to=max(meta_reads$concentration), by=0.235)


# Create new data frame where only concentration changes
# qPCR score (factor) is set as one level

d1 <- data.frame(run=rep('Run 1', length(range)),
                 plate=rep('1', length(range)),
                 qPCR_score=rep(1, length(range)),
                 concentration=range)


# Get predictions for this new dataset
fit <- data.frame(predictSE(m3.nb, newdata=d1, 
                            se.fit = TRUE, print.matrix=T)) 


# Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

# Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

dat.conc <- cbind(d1, fit, ciu, cil) # This is now the data frame


# Plot showing relationship between eDNA concentration and proportion
# of GCN reads with predicted values from model

p12 <- ggplot() + ggtitle('(b)')
p12 <- p12 + geom_jitter(aes(x=meta_reads$concentration, y=meta_reads$GCNprop), shape=1, cex=0.95, show.legend=FALSE)
p12 <- p12 + geom_line(aes(x=dat.conc$concentration, y=dat.conc$fit), colour="red", size = 1, show.legend=FALSE)
p12 <- p12 + geom_ribbon(aes(x=dat.conc$concentration, ymin = dat.conc$cil, ymax = dat.conc$ciu), alpha = 0.25, show.legend=FALSE)
p12 <- p12 + scale_colour_gradient(limits=c(0, 100), low="gray30", high="gray80")
p12 <- p12 + coord_cartesian(ylim=c(0,100))
p12 <- p12 + scale_x_continuous(breaks=seq(0,125,25))
p12 <- p12 + labs(x = expression(paste("Post-PCR eDNA concentration (ng/",mu,"l)")), y = "")
p12 <- p12 + theme(panel.background = element_rect(fill = 'white'), 
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   plot.title =element_text(face="bold", hjust=0),
                   text = element_text(size=24),
                   legend.position="none")
p12


# SHOW ALL PLOTS
# grid.arrange() is a function within the gridExtra package that allows to
# arrange how ggplots are viewed, similar to par(mfrow=c())
# Within this function, I can use arrangeGrob() to specify how many plots I
# want on each row

grid.arrange(arrangeGrob(p11,p12, nrow=1))

