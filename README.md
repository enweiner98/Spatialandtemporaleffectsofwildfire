REPOSITORY ORIENTATION

This repository contains 2 types of files:
   > Data files (.csv)
   > Coding files (.R)

Open **'Thesis.Rproj'** in R Studio

SETTING UP COVARIATE DATA 
From working directory, open **'autoOcc_setup.R'** file, which contains code to:
  > Import covariate data +
  > Center and scale covariates for spatial and temporal models +
  > Create "covariate frames" ("covFrames") for use in occupancy models (later referenced in species-specific .R files)
  > Perform a correlation analysis
**Run lines 4-210** to: Import covariate data, scale and center covariates, and build covariate frames
**Run lines 214-278** to: Build and view correlation matrices and plots for continuous temporal and spatial covariates

MODELING OCCUPANCY FOR EACH TARGET SPECIES
From working directory, open **'coyote_final.R'** file, which contains code to:
  > Import spatial and temporal detection data + 
  > Build a species-specific dataframe for structured for use with 'autoOcc' +  
  > Run spatial and temporal occupancy models and perform AIC comparisons + 
  > Produce plots and summary statistics
There is a separate .R file for analysis of each target species (10 files total), but all species-specific .R files are structured in the same order: 
  1. Import spatial and temporal detection data from working directory  + 
  2. Produce species-specific summary statistics +
  3. Create a species-specific autoOcc dataframe + 
  4. Run spatial models, perform AIC comparison, and produce covariate plots +  
  5. Run temporal models, perform AIC comparison, and produce covariate plots
