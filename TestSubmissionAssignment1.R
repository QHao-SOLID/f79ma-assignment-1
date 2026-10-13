###################################################
#   Script to test submissions for assignment    ##
###################################################
#   DO NOT INCLUDE THIS CODE INTO YOUR R-FILE   ###
###################################################

# remove all R objects in the current workspace
# never include this line in the R file for your
# submission
rm(list = ls())

# close all graphical displays 
graphics.off()

# Load R file 
# The following line needs to be changed if your
# R file is not called "F79MA_Assignment1.R"
source("F79MA_Assignment1.R")

# Check that all answers for this project are
# correctly defined
print(c(Quantity1, Quantity2, Quantity3, Quantity4,
        Quantity5, Quantity6))