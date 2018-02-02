#' ---
#' title: "Kappa Coefficient"
#' author: "Lynsey Rebecca Harper, 0901795"
#' date: "29th July 2015"
#' ---

# Remove anything currently in R's brain
rm(list=ls())

## set working directory
setwd("/Users/Lynsey/Documents/Masters Project/Data Analysis/Torchlight Survey")


#'
#' Agreement between torchlight survey and eDNA analysis would have been 
#' measured using Cohen’s kappa coefficient (Cohen 1960):
#' 
#' $κ = \frac {P_a - P_e}{1 - P_e}$
#' 
#' where $P_a$ is the relative agreement among rates and $P_e$ is the 
#' hypothetical probability of chance agreement, using the observed data 
#' to calculate the probabilities of each method randomly giving a positive 
#' detection. If the methods are in complete agreement, then $κ$ = 1. 
#' If there is no agreement other than that which would be expected by chance, 
#' $κ$ = 0. 
#' 

#' ---
#' 
#' ### Thomson PCR and Torchlight Survey
#' 

# Use matrix function to enter 2x2 contingency table into R

## Visit 1
visit1 <- matrix(c(22, 2, 0, 0), 2, 2, byrow = TRUE)


## Visit 2
visit2 <- matrix(c(24, 0, 0, 0), 2, 2, byrow = TRUE)


## Visit 3
visit3 <- matrix(c(21, 3, 0, 0), 2, 2, byrow = TRUE)


## Visit 4
visit4 <- matrix(c(23, 1, 0, 0), 2, 2, byrow = TRUE)


## Visit 5
visit5 <- matrix(c(20, 4, 0, 0), 2, 2, byrow = TRUE)

## Load psych package which contains cohen.kappa() function to calculate
## kappa coefficient
library(psych)
```{r, warning=F}

cohen.kappa(visit1, alpha = 0.05)
cohen.kappa(visit2, alpha = 0.05)
cohen.kappa(visit3, alpha = 0.05)
cohen.kappa(visit4, alpha = 0.05)
cohen.kappa(visit5, alpha = 0.05)


#'
#' Estimates of 0 were produced for each visit due to a prevalence effect.
#' A prevalence effect occurs where there is lack of disagreement between 
#' two methods, resulting in $P_a$ = $P_e$. This is also known as the kappa 
#' paradox and can lead to misleading results. I must calculate additional 
#' values to report alongside $κ$ due to the prevalence effect. $κ$ can be 
#' adjusted for these imbalances using the Prevalence Adjusted Bias Adjusted 
#' Kappa (PABAK) (Byrt, Bishop & Carlin 1993). This statistic depends solely 
#' on the observed proportion of agreement between the two methods:
#' 
#' PABAK = 2$P_a$ – 1
#' 
#' 
#' I have made my own function in a separate R file to calculate $κ$ and other values. 
#' The other values are:
#' 
#' + Proportion of agreement ($P_a$)
#' 
#' + Expected proportion of agreement ($P_e$)
#' 
#' + Proportion of positive agreement ($P_{pos}$) and negative agreement ($P_{neg}$)
#' 
#' + Prevalence index ($P_{index}$)
#' 
#' + Bias index ($B_{index}$)
#' 
#' + Prevalence Adjusted Bias Adjusted Kappa (PABAK)
#' 

## Source file containing kappa coefficient function
source("kappa.coefficient.function.R")


## VISIT 1

## Set paramaters
a <- 22
b <- 2
c <- 0
d <- 0
N <- 24

visit1.k <- kappa.function(visit1)


## VISIT 2

## Set paramaters
a <- 24
b <- 0
c <- 0
d <- 0
N <- 24

visit2.k <- kappa.function(visit2)


## VISIT 3

## Set paramaters
a <- 21
b <- 3
c <- 0
d <- 0
N <- 24

visit3.k <- kappa.function(visit3)


## VISIT 4

## Set paramaters
a <- 23
b <- 1
c <- 0
d <- 0
N <- 24

visit4.k <- kappa.function(visit4)


## VISIT 5

## Set paramaters
a <- 20
b <- 4
c <- 0
d <- 0
N <- 24

visit5.k <- kappa.function(visit5)


#' 
#' I will now summarise the values for each torchlight survey versus Thomson
#' PCR eDNA analysis in a table.
#' 

kappa <- matrix(c(visit1.k, visit2.k, visit3.k, visit4.k, visit5.k), ncol=8, byrow=TRUE)
colnames(kappa) <- c("P(a)","P(e)","kappa", "P(pos)", "P(neg)", "P(index)", "B(index)", "PABAK")
rownames(kappa) <- c("Visit 1","Visit 2","Visit 3","Visit 4", "Visit 5")
kappa <- as.table(kappa)
kappa



#' ---
#' 
#' ### Warwick PCR and Torchlight Survey
#'

## VISIT 1

a <- 0
b <- 0
c <- 12
d <- 2
N <- 14

visit1.k <- kappa.function(visit1)


## VISIT 2

a <- 0
b <- 0
c <- 14
d <- 0
N <- 14

visit2.k <- kappa.function(visit2)


## VISIT 3

a <- 0
b <- 0
c <- 12
d <- 2
N <- 14

visit3.k <- kappa.function(visit3)


## VISIT 4

a <- 0
b <- 0
c <- 13
d <- 1
N <- 14

visit4.k <- kappa.function(visit4)


## VISIT 5

a <- 0
b <- 0
c <- 11
d <- 3
N <- 14

visit5.k <- kappa.function(visit5)


#' 
#' I will now summarise the values for each torchlight survey versus Warwick
#' PCR eDNA analysis in a table.
#' 

kappa <- matrix(c(visit1.k, visit2.k, visit3.k, visit4.k, visit5.k), ncol=8, byrow=TRUE)
colnames(kappa) <- c("P(a)","P(e)","kappa", "P(pos)", "P(neg)", "P(index)", "B(index)", "PABAK")
rownames(kappa) <- c("Visit 1","Visit 2","Visit 3","Visit 4", "Visit 5")
kappa <- as.table(kappa)
kappa


