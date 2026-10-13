# R-code to generate data for F79MA Assignment 1. This code MUST be 
# placed at the start of your own R-script. You must edit
# the argument to the set.seed( ) function to fit your own
# registration number
#
# Please also ensure you have created all required variables
# (Quantity1, Quantity2, Quantity3, Quantity4, Quantity5, Quantity6)
# and include code for all of your plots
#

set.seed(7830)

RNGkind(sample.kind = "Rejection")
lambda0 <- 2
lambda1 <- 3
lambda2 <- 4
#
simulate_YA_YB <- function(n){ # n is the sample size
  set.seed(1234) ## Change this number
  X0 <- rpois(n,lambda0)
  X1 <- rpois(n,lambda1)
  X2 <- rpois(n,lambda2)
  YA <- X0 + X1 # sample of tremor numbers in zone A
  YB <- X0 + X2 # sample of tremor numbers in zone B
  return(list(YA = YA,YB = YB))
}
#####################################################################
# Please insert your R code after this line
#####################################################################