#' ---
#' title: "Kappa Coefficient Function"
#' author: "Lynsey Rebecca Harper, 0901795"
#' date: "29th July 2015"
#' ---
#' 
#' 
#' Run one step of kappa coefficient calculation.
#' 
#' Arguments:
#' 
#' - $P_a$ -- the observed agreement
#' - $P_e$ -- the expected agreement
#' 
#' Returns:
#' 
#' - the kappa coefficient ($κ$) i.e. level of agreement between two methods
#'
#' The function() command has 1 line of code that defines what it does.
#' For multiple commands, curly arrow syntax must be used.


## Set up function doing work in simulation with added improvements to handle 
## data frame
kappa.function <- function(Pa, Pe)
{
  ## Calculation for Pa
  Pa <- (a+d)/N
  
  ## Calculation for Pe
  Pe <- ((a+c)*(a+b) + (b+d)*(c+d))/N^2
  
  ## Calculation for kappa coefficient
  kappa <- (Pa-Pe)/(1-Pe)
  
  ## Check for kappa = 0 and specify equations to be used if so
  if(kappa == "NaN")      ## If NaN value returned for kappa...
  {
    ## Replace NaN value with 0
    kappa <- replace(kappa, is.na(kappa), 0)
  }
    
  ## Proportion of positive agreement (ppos)
  ppos <- (2*a)/(N+a-d)
    
  ## Proportion of negative agreement (pneg)
  pneg <- (2*d)/(N-a+d)
    
  ## If pneg returns NaN value, replace with 0
  if (pneg == "NaN")
  {
    ## Replace NaN value with 0
    pneg <- replace(pneg, is.na(pneg), 0)
  }
    
  ## Prevalence index (pi)
  pi <- (a-d)/N
    
  # Bias index (bi)
  bi <- (b-c)/N
    
  ## Prevalence Adjusted Bias Adjusted Kappa (PABAK)
  pabak <- 2 * Pa - 1
  
  ## store values returned by kappa.function
  summary <- c(Pa, Pe, kappa, ppos, pneg, pi, bi, pabak)
  
  ## return values
  summary
}


## Make sure function is only using information given to it in arguments
## Use findGlobals function to check that it doesn't use any external information
library(codetools)
findGlobals(kappa.function, merge=FALSE)

