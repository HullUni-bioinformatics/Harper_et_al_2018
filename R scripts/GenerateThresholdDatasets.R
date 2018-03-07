#' ---
#' Title: "Generate metabarcoding datasets with different false positive thresholds"
#' Author: "Lynsey Rebecca Harper"
#' Date: "5th December 2016"
#' ---
#' 
#' 
#' A total of 532 samples previously analysed by qPCR for great crested newt 
#' (GCN) were re-analysed by eDNA metabarcoding to compare method sensitivity 
#' and evaluate this method for detection of vertebrate communities in 
#' freshwater ponds.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. 
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

origBLAST <- read.csv("../Data/metabarcoding_PrepOutput.csv",
                      header=TRUE)

summary(origBLAST[,1:6])
head(origBLAST)
names(origBLAST)
str(origBLAST)


## Import the prepared unassigned file:

unassBLAST <- read.csv("../Data/Unassigned_PrepOutput.csv",
                       header=TRUE)

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

## Remove '_.nc.blast' after sample ID in both data frames
colnames(unassBLAST) <- gsub("_.nc.blast.blast", "", colnames(unassBLAST))
colnames(origBLAST) <- gsub("_.nc.blast", "", colnames(origBLAST))

## Bind data frames
merged_df <- rbind(origBLAST, unassBLAST)

## Merge read counts by taxonomic assignment
## Anything with the same name in each data frame will be merged and
## anything new in the unassigned BLAST added to the original BLAST
## data frame
mergedBLAST <- ddply(merged_df, .(Tax_assignment), numcolwise(sum))

## Export as csv file
write.csv(mergedBLAST, "../Data/MergedBLAST.csv", row.names=FALSE)


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
#' - *Ceratocystis* (Fungi)
#' - *Cloen dipterum* (invert)
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

## Reset row names
row.names(true_assign) <- NULL

## Make taxonomic assignment a character column
true_assign$Tax_assignment <- as.character(true_assign$Tax_assignment)

## Rename the genus/family assignments
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
#' ## Investigate effect of False Positive Sequence Thresholds
#' 
#' Now the data is in a form that can be manipulated easily, I must tidy the
#' metabarcoding dataset using false positive sequence thresholds. These 
#' thresholds remove assignments if they don't exceed a designated read count.
#' Thresholds are determined by examining PCR positive controls. Any DNA
#' detected in these controls which is not cichlid DNA, can be classed as
#' a false positive.
#' 
#' Arguably, this is a more effective way to deal with contamination and 
#' false positives in metabarcoding as negative controls have no template 
#' DNA for contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially and paints a false picture of contamination.
#' 
#' 
#' There are two main options:
#' 
#' 1. Apply a blanket false positive sequence threshold to remove a desirable
#'    proportion of contamination from the dataset. Any assignment which
#'    exceeds this threshold will be removed.
#'    
#' 2. Apply species-specific false positive thresholds i.e. the frequency 
#'    required to remove a given species from the PCR positive control. Only
#'    that species will be removed. 
#' 
#' A blanket threshold will not remove major contaminants, such as human. 
#' Conversely, a species-specific threshold for human would ensure this 
#' contamination is also removed from the positive PCR control and presumably 
#' in my eDNA samples. 
#' 
#' However, even with application of false positive thresholds, some 
#' contamination will still be observed in PCR negative controls. This may 
#' come from the environment, other lab work or consumables.
#' 
#' We will compare both approaches. I will test a range of blanket false
#' positive sequence thresholds (0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.3)
#' previously examined using an excel spreadsheet. A false positive threshold 
#' of 0.005 removes 97.5% contamination from my dataset but a threshold of 
#' 0.3 was required to entirely remove contamination.
#' 
#' I will examine the effect increasing blanket sequence threshold has on
#' my downstream analyses where a species-specific threshold approach has 
#' been used.
#'  
#' Calculate the total read counts, frequency of each taxonomic in an 
#' assignment, and then apply a condition to the dataframe to set any 
#' assignments to a species below the false positive threshold to 0.
#' 

## Create new dataframe to contain only samples and read counts
reads <- subset(vert_seq[,2:761])

## Calculate total number of reads for each sample
counts <- rbind(reads, colSums(reads))

## Create new dataframe containing the frequency of reads in each sample
freq <- counts/c(counts[100,])

## Combine read frequencies with species names
species <- data.frame(vert_seq$Tax_assignment)
colnames(species) <- "Tax_assignment"

## Create an empty row which will contain total read frequency for each
## sample
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"

## Bind data frames
species <- rbind(species, total)
freq <- cbind(species, freq)


## Create a series of dataframes with different sequence thresholds
## Replace all assignments that are less than or equal to the false 
## positive threshold for a given species with zero.
## Don't include negative controls as any reads in these constitute 
## true contamination

## 0.0005
threshold_0.0005 <- freq
threshold_0.0005[,c(2:509,624:761)][threshold_0.0005[,c(2:509,624:761)] <= 0.0005] <- 0
  
## 0.001
threshold_0.001 <- freq
threshold_0.001[,c(2:509,624:761)][threshold_0.001[,c(2:509,624:761)] <= 0.001] <- 0

## 0.005
threshold_0.005 <- freq
threshold_0.005[,c(2:509,624:761)][threshold_0.005[,c(2:509,624:761)] <= 0.005] <- 0
      
## 0.01
threshold_0.01 <- freq
threshold_0.01[,c(2:509,624:761)][threshold_0.01[,c(2:509,624:761)] <= 0.01] <- 0

## 0.05
threshold_0.05 <- freq
threshold_0.05[,c(2:509,624:761)][threshold_0.05[,c(2:509,624:761)] <= 0.05] <- 0

## 0.1
threshold_0.1 <- freq
threshold_0.1[,c(2:509,624:761)][threshold_0.1[,c(2:509,624:761)] <= 0.1] <- 0

## 0.3
threshold_0.3 <- freq
threshold_0.3[,c(2:509,624:761)][threshold_0.3[,c(2:509,624:761)] <= 0.3] <- 0


#'
#' Now, I need to get each dataframe back into read count format with 
#' taxonomic assignment for further analyses.
#' 


##########
# 0.0005 #
##########

## Remove last row containing total frequencies and add total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
sample_threshold_0.0005 <- threshold_0.0005[1:99,1:761]
total_counts <- counts[100,]

## Create empty row which will contain total read counts
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

## Bind data frames using smartbind (recognises same column names
## across data frames) into temporary data frame
conversion <- smartbind(sample_threshold_0.0005, total_counts)

## Convert frequences back to read counts
threshold_0.0005 <- conversion*c(conversion[100,])

## Remove 'Total' row and empty Tax_assignment row and add metadata 
## to final dataframe
threshold_0.0005 <- threshold_0.0005[-c(100), -c(1)]
threshold_0.0005 <- cbind(vert_seq[,1], threshold_0.0005)

## Replace any NA values in dataframe with 0
threshold_0.0005[is.na(threshold_0.0005)] <- 0

## Export as csv file to data working directory for visual inspection
## and comparison with initial input file
write.csv(threshold_0.0005, "../Data/0.0005ThresholdApplied.csv",
          row.names=FALSE)


#########
# 0.001 #
#########

## Remove last row containing total frequencies and add total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
sample_threshold_0.001 <- threshold_0.001[1:99,1:761]
total_counts <- counts[100,]

## Create empty row which will contain total read counts
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

## Bind data frames using smartbind (recognises same column names
## across data frames) into temporary data frame
conversion <- smartbind(sample_threshold_0.001, total_counts)

## Convert frequences back to read counts
threshold_0.001 <- conversion*c(conversion[100,])

## Remove 'Total' row and empty Tax_assignment row and add metadata 
## to final dataframe
threshold_0.001 <- threshold_0.001[-c(100), -c(1)]
threshold_0.001 <- cbind(vert_seq[,1], threshold_0.001)

## Replace any NA values in dataframe with 0
threshold_0.001[is.na(threshold_0.001)] <- 0

## Export as csv file to data working directory for visual inspection
## and comparison with initial input file
write.csv(threshold_0.001, "../Data/0.001ThresholdApplied.csv",
          row.names=FALSE)


#########
# 0.005 #
#########

## Remove last row containing total frequencies and add total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
sample_threshold_0.005 <- threshold_0.005[1:99,1:761]
total_counts <- counts[100,]

## Create empty row which will contain total read counts
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

## Bind data frames using smartbind (recognises same column names
## across data frames) into temporary data frame
conversion <- smartbind(sample_threshold_0.005, total_counts)

## Convert frequences back to read counts
threshold_0.005 <- conversion*c(conversion[100,])

## Remove 'Total' row and empty Tax_assignment row and add metadata 
## to final dataframe
threshold_0.005 <- threshold_0.005[-c(100), -c(1)]
threshold_0.005 <- cbind(vert_seq[,1], threshold_0.005)

## Replace any NA values in dataframe with 0
threshold_0.005[is.na(threshold_0.005)] <- 0

## Export as csv file to data working directory for visual inspection
## and comparison with initial input file
write.csv(threshold_0.005, "../Data/0.005ThresholdApplied.csv",
          row.names=FALSE)


#########
# 0.01  #
#########

## Remove last row containing total frequencies and add total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
sample_threshold_0.01 <- threshold_0.01[1:99,1:761]
total_counts <- counts[100,]

## Create empty row which will contain total read counts
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

## Bind data frames using smartbind (recognises same column names
## across data frames) into temporary data frame
conversion <- smartbind(sample_threshold_0.01, total_counts)

## Convert frequences back to read counts
threshold_0.01 <- conversion*c(conversion[100,])

## Remove 'Total' row and empty Tax_assignment row and add metadata 
## to final dataframe
threshold_0.01 <- threshold_0.01[-c(100), -c(1)]
threshold_0.01 <- cbind(vert_seq[,1], threshold_0.01)

## Replace any NA values in dataframe with 0
threshold_0.01[is.na(threshold_0.01)] <- 0

## Export as csv file to data working directory for visual inspection
## and comparison with initial input file
write.csv(threshold_0.01, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 1/Data/0.01ThresholdApplied.csv",
          row.names=FALSE)


########
# 0.05 #
########

## Remove last row containing total frequencies and add total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
sample_threshold_0.05 <- threshold_0.05[1:99,1:761]
total_counts <- counts[100,]

## Create empty row which will contain total read counts
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

## Bind data frames using smartbind (recognises same column names
## across data frames) into temporary data frame
conversion <- smartbind(sample_threshold_0.05, total_counts)

## Convert frequences back to read counts
threshold_0.05 <- conversion*c(conversion[100,])

## Remove 'Total' row and empty Tax_assignment row and add metadata 
## to final dataframe
threshold_0.05 <- threshold_0.05[-c(100), -c(1)]
threshold_0.05 <- cbind(vert_seq[,1], threshold_0.05)

## Replace any NA values in dataframe with 0
threshold_0.05[is.na(threshold_0.05)] <- 0

## Export as csv file to data working directory for visual inspection
## and comparison with initial input file
write.csv(threshold_0.05, "../Data/0.05ThresholdApplied.csv",
          row.names=FALSE)


#######
# 0.1 #
#######

## Remove last row containing total frequencies and add total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
sample_threshold_0.1 <- threshold_0.1[1:99,1:761]
total_counts <- counts[100,]

## Create empty row which will contain total read counts
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

## Bind data frames using smartbind (recognises same column names
## across data frames) into temporary data frame
conversion <- smartbind(sample_threshold_0.1, total_counts)

## Convert frequences back to read counts
threshold_0.1 <- conversion*c(conversion[100,])

## Remove 'Total' row and empty Tax_assignment row and add metadata 
## to final dataframe
threshold_0.1 <- threshold_0.1[-c(100), -c(1)]
threshold_0.1 <- cbind(vert_seq[,1], threshold_0.1)

## Replace any NA values in dataframe with 0
threshold_0.1[is.na(threshold_0.1)] <- 0

## Export as csv file to data working directory for visual inspection
## and comparison with initial input file
write.csv(threshold_0.1, "../Data/0.1ThresholdApplied.csv",
          row.names=FALSE)


#######
# 0.3 #
#######

## Remove last row containing total frequencies and add total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
sample_threshold_0.3 <- threshold_0.3[1:99,1:761]
total_counts <- counts[100,]

## Create empty row which will contain total read counts
total <- data.frame("Total")
colnames(total) <- "Tax_assignment"
total_counts <- cbind(total, total_counts)

## Bind data frames using smartbind (recognises same column names
## across data frames) into temporary data frame
conversion <- smartbind(sample_threshold_0.3, total_counts)

## Convert frequences back to read counts
threshold_0.3 <- conversion*c(conversion[100,])

## Remove 'Total' row and empty Tax_assignment row and add metadata 
## to final dataframe
threshold_0.3 <- threshold_0.3[-c(100), -c(1)]
threshold_0.3 <- cbind(vert_seq[,1], threshold_0.3)

## Replace any NA values in dataframe with 0
threshold_0.3[is.na(threshold_0.3)] <- 0

## Export as csv file to data working directory for visual inspection
## and comparison with initial input file
write.csv(threshold_0.3, "../Data/0.3ThresholdApplied.csv",
          row.names=FALSE)


#'
#' Now add some metadata on assignment levels and keep only species 
#' assignments.
#' 
#' Create new data frames for presence/absence detection of species in
#' eDNA samples (1 = presence). 
#' 
#' Exclude cichlid assignment as this is the positive control. 
#' Exclude positive and negative PCR controls.
#' 

## Import prepared output file
prep_df <- read.csv("../Data/ThresholdAnalysis_PrepInfo.csv",
                    header=TRUE)


##########
# 0.0005 #
##########

## Rename first column in threshold dataframes
colnames(threshold_0.0005)[1] <- "Tax_assignment"

## Merge dataframes
merge_0.0005 <- merge(prep_df,threshold_0.0005, all=TRUE)
head(merge_0.0005[,1:6])

## Remove cichlid and PCR controls
binary_0.0005 <- data.frame(merge_0.0005[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0.0005[,3:534][binary_0.0005[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0.0005, binary_0.0005$AssignmentLevel=="Species")
new_0.0005 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0.0005) <- new_0.0005[,1]
new_0.0005 <- new_0.0005[,-1]
pre_0.0005 <- data.frame(t(new_0.0005))

## Move the sample names into the first column
pre_0.0005 <- cbind(sample = rownames(pre_0.0005), pre_0.0005, row.names = NULL)

## Import metadata for samples
sample_metadata <- read.csv("../Data/GCN_metadata.csv", header=TRUE)

## Make sure sample columns in each data frame match type
str(pre_0.0005)
str(sample_metadata)

## Merge data frames
final_0.0005 <- merge(sample_metadata, pre_0.0005, all=TRUE)

## Check new data frame
summary(final_0.0005)
head(final_0.0005)
names(final_0.0005)
str(final_0.0005)

## Export as csv file
write.csv(final_0.0005, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_threshold_0.0005.csv", 
          row.names=FALSE)

## REPEAT FOR DIFFERENT THRESHOLDS


#########
# 0.001 #
#########

## Rename first column in threshold dataframes
colnames(threshold_0.001)[1] <- "Tax_assignment"

## Merge dataframes
merge_0.001 <- merge(prep_df,threshold_0.001, all=TRUE)
head(merge_0.001[,1:6])

## Remove cichlid and PCR controls
binary_0.001 <- data.frame(merge_0.001[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0.001[,3:534][binary_0.001[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0.001, binary_0.001$AssignmentLevel=="Species")
new_0.001 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0.001) <- new_0.001[,1]
new_0.001 <- new_0.001[,-1]
pre_0.001 <- data.frame(t(new_0.001))

## Move the sample names into the first column
pre_0.001 <- cbind(sample = rownames(pre_0.001), pre_0.001, row.names = NULL)

## Make sure sample columns in each data frame match type
str(pre_0.001)
str(sample_metadata)

## Merge data frames
final_0.001 <- merge(sample_metadata, pre_0.001, all=TRUE)

## Check new data frame
summary(final_0.001)
head(final_0.001)
names(final_0.001)
str(final_0.001)

## Export as csv file
write.csv(final_0.001, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_threshold_0.001.csv", 
          row.names=FALSE)


#########
# 0.005 #
#########

## Rename first column in threshold dataframes
colnames(threshold_0.005)[1] <- "Tax_assignment"

## Merge dataframes
merge_0.005 <- merge(prep_df,threshold_0.005, all=TRUE)
head(merge_0.005[,1:6])

## Remove cichlid and PCR controls
binary_0.005 <- data.frame(merge_0.005[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0.005[,3:534][binary_0.005[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0.005, binary_0.005$AssignmentLevel=="Species")
new_0.005 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0.005) <- new_0.005[,1]
new_0.005 <- new_0.005[,-1]
pre_0.005 <- data.frame(t(new_0.005))

## Move the sample names into the first column
pre_0.005 <- cbind(sample = rownames(pre_0.005), pre_0.005, row.names = NULL)

## Make sure sample columns in each data frame match type
str(pre_0.005)
str(sample_metadata)

## Merge data frames
final_0.005 <- merge(sample_metadata, pre_0.005, all=TRUE)

## Check new data frame
summary(final_0.005)
head(final_0.005)
names(final_0.005)
str(final_0.005)

## Export as csv file
write.csv(final_0.005, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_threshold_0.005.csv", 
          row.names=FALSE)


########
# 0.01 #
########

## Rename first column in threshold dataframes
colnames(threshold_0.01)[1] <- "Tax_assignment"

## Merge dataframes
merge_0.01 <- merge(prep_df,threshold_0.01, all=TRUE)
head(merge_0.01[,1:6])

## Remove cichlid and PCR controls
binary_0.01 <- data.frame(merge_0.01[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0.01[,3:534][binary_0.01[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0.01, binary_0.01$AssignmentLevel=="Species")
new_0.01 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0.01) <- new_0.01[,1]
new_0.01 <- new_0.01[,-1]
pre_0.01 <- data.frame(t(new_0.01))

## Move the sample names into the first column
pre_0.01 <- cbind(sample = rownames(pre_0.01), pre_0.01, row.names = NULL)

## Make sure sample columns in each data frame match type
str(pre_0.01)
str(sample_metadata)

## Merge data frames
final_0.01 <- merge(sample_metadata, pre_0.01, all=TRUE)

## Check new data frame
summary(final_0.01)
head(final_0.01)
names(final_0.01)
str(final_0.01)

## Export as csv file
write.csv(final_0.01, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_threshold_0.01.csv", 
          row.names=FALSE)


########
# 0.05 #
########

## Rename first column in threshold dataframes
colnames(threshold_0.05)[1] <- "Tax_assignment"

## Merge dataframes
merge_0.05 <- merge(prep_df,threshold_0.05, all=TRUE)
head(merge_0.05[,1:6])

## Remove cichlid and PCR controls
binary_0.05 <- data.frame(merge_0.05[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0.05[,3:534][binary_0.05[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0.05, binary_0.05$AssignmentLevel=="Species")
new_0.05 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0.05) <- new_0.05[,1]
new_0.05 <- new_0.05[,-1]
pre_0.05 <- data.frame(t(new_0.05))

## Move the sample names into the first column
pre_0.05 <- cbind(sample = rownames(pre_0.05), pre_0.05, row.names = NULL)

## Make sure sample columns in each data frame match type
str(pre_0.05)
str(sample_metadata)

## Merge data frames
final_0.05 <- merge(sample_metadata, pre_0.05, all=TRUE)

## Check new data frame
summary(final_0.05)
head(final_0.05)
names(final_0.05)
str(final_0.05)

## Export as csv file
write.csv(final_0.05, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_threshold_0.05.csv", 
          row.names=FALSE)


#######
# 0.1 #
#######

## Rename first column in threshold dataframes
colnames(threshold_0.1)[1] <- "Tax_assignment"

## Merge dataframes
merge_0.1 <- merge(prep_df,threshold_0.1, all=TRUE)
head(merge_0.1[,1:6])

## Remove cichlid and PCR controls
binary_0.1 <- data.frame(merge_0.1[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0.1[,3:534][binary_0.1[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0.1, binary_0.1$AssignmentLevel=="Species")
new_0.1 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0.1) <- new_0.1[,1]
new_0.1 <- new_0.1[,-1]
pre_0.1 <- data.frame(t(new_0.1))

## Move the sample names into the first column
pre_0.1 <- cbind(sample = rownames(pre_0.1), pre_0.1, row.names = NULL)

## Make sure sample columns in each data frame match type
str(pre_0.1)
str(sample_metadata)

## Merge data frames
final_0.1 <- merge(sample_metadata, pre_0.1, all=TRUE)

## Check new data frame
summary(final_0.1)
head(final_0.1)
names(final_0.1)
str(final_0.1)

## Export as csv file
write.csv(final_0.1, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_threshold_0.1.csv", 
          row.names=FALSE)


#######
# 0.3 #
#######

## Rename first column in threshold dataframes
colnames(threshold_0.3)[1] <- "Tax_assignment"

## Merge dataframes
merge_0.3 <- merge(prep_df,threshold_0.3, all=TRUE)
head(merge_0.3[,1:6])

## Remove cichlid and PCR controls
binary_0.3 <- data.frame(merge_0.3[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0.3[,3:534][binary_0.3[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0.3, binary_0.3$AssignmentLevel=="Species")
new_0.3 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0.3) <- new_0.3[,1]
new_0.3 <- new_0.3[,-1]
pre_0.3 <- data.frame(t(new_0.3))

## Move the sample names into the first column
pre_0.3 <- cbind(sample = rownames(pre_0.3), pre_0.3, row.names = NULL)

## Make sure sample columns in each data frame match type
str(pre_0.3)
str(sample_metadata)

## Merge data frames
final_0.3 <- merge(sample_metadata, pre_0.3, all=TRUE)

## Check new data frame
summary(final_0.3)
head(final_0.3)
names(final_0.3)
str(final_0.3)

## Export as csv file
write.csv(final_0.3, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_threshold_0.3.csv", 
          row.names=FALSE)


################
# No threshold #
################

#' 
#' NB: vert_seq dataframe is dataset with no threshold applied.
#'

## Rename vert_seq
threshold_0 <- vert_seq

## Merge dataframes
merge_0 <- merge(prep_df,threshold_0, all=TRUE)
head(merge_0[,1:6])

## Remove cichlid and PCR controls
binary_0 <- data.frame(merge_0[-c(83), c(1:510,625:648)])

## Convert read counts to presence or absence
binary_0[,3:534][binary_0[,3:534] > 0] <- 1

## Remove assignment levels higher than species. These don't provide 
## detailed information on biodiversity unless a detected genus or
## family only contains one species in the UK. These instances have 
## already been identified and changed to species.
detections <- subset(binary_0, binary_0$AssignmentLevel=="Species")
new_0 <- detections[,c(1,3:length(detections))]

## Make assignment names row names and transpose the dataframe
rownames(new_0) <- new_0[,1]
new_0 <- new_0[,-1]
pre_0 <- data.frame(t(new_0))

## Move the sample names into the first column
pre_0 <- cbind(sample = rownames(pre_0), pre_0, row.names = NULL)

## Make sure sample columns in each data frame match type
str(pre_0)
str(sample_metadata)

## Merge data frames
final_0 <- merge(sample_metadata, pre_0, all=TRUE)

## Check new data frame
summary(final_0)
head(final_0)
names(final_0)
str(final_0)

## Export as csv file
write.csv(final_0, "C:/Users/518740/Dropbox/PhD Write-up/Chapter 2 - GCN interactions/Data/data_no_threshold.csv", 
          row.names=FALSE)
