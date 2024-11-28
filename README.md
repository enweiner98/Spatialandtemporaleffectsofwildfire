REPOSITORY ORIENTATION

This repository contains 3 types of files:
   1. R project file (.Rproj)
   2. Data files (.csv)
   3. Coding files (.R)

> Open **'Thesis.Rproj'** in R Studio

SETTING UP COVARIATE DATA 
> From working directory, open **'autoOcc_setup.R'** file, which contains code to:
  1. Import covariate data +
  2. Center and scale covariates for spatial and temporal models +
  3. Create "covariate frames" ("covFrames") for use in occupancy models (later referenced in species-specific .R files)
  4. Perform a correlation analysis

> To import covariate data, scale and center covariates, and build covariate frames, **run lines 4-210** 

> To build and view correlation matrices and plots for continuous temporal and spatial covariates, **run lines 214-278** 

MODELING OCCUPANCY FOR EACH TARGET SPECIES
> From working directory, open **'coyote_final.R'** file, which contains code to:
  1. Import spatial and temporal detection data + 
  2. Build a species-specific dataframe for structured for use with 'autoOcc' +  
  3. Run spatial and temporal occupancy models and perform AIC comparisons + 
  4. Produce plots and summary statistics
> There is a separate .R file for analysis of each target species (10 files total), but all species-specific .R files are structured in the same order: 
  1. Import spatial and temporal detection data from working directory  + 
  2. Produce species-specific summary statistics +
  3. Create a species-specific autoOcc dataframe + 
  4. Run spatial models, perform AIC comparison, and produce covariate plots +  
  5. Run temporal models, perform AIC comparison, and produce covariate plots

All required data files should be present in the repository. Make sure to run the code in autoOcc_setup.R before attempting to run code in species-specific files.
