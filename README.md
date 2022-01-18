# A <em> de novo <em> genome assembly and annotation of the southern flying squirrel (Glaucomys volans) (Wolf et al. 2021)
### This repository the raw data and R scripts used to perform analyses and generate figures for my manuscript published in Ecology and Evolution. Detailed below are the steps this script walks you through, starting with raw data, data manipulation/filtering, running statistical models, and generating figures.

If you are interested in the manuscript, see https://doi.org/10.1002/ece3.8106

1. Read in raw data in .csv format
1. Reducing dataset to relevant variables
1. Ensuring all relevant data columns are recognized as factors 
1. Create a series of models, building up from null model to full model by adding one term at a time
1. Run final Generalized Linear Model 
1. Generate figures and combine them using the patchwork package in R
