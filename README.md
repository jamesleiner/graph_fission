# Graph Fission and Cross-Validation

This repository contains code to replicate all figures and experiments found in our AISTATS submission Graph Fission and Cross-Validation

## Pre-requisites
Prior to running the code, the following packages are required:
* MASS
* matrixcalc
* pracma
* dplyr
* tidyr
* plotly
* genlasso

Code built in R version 4.0.3. Note the code uses the function "mclapply", which will not work in a Windows environment. This can be changed to lapply, but then the code will not be parallelized and will take signficiantly longer to run. 

## Instructions
To generate experiments and figures used in the paper
* Run the script run_experiments.R
  * Note that the code was built in a multi-core environment with 128 cores available. Change the option "mc.cores=128" in the mclapply calls to whatever number of cores is available in your computing environment.
* After all results have been saved, run the script create_graphs.R to create the figures used in the paper


* The experiments for the interactive multiple testing application (section 3) need to be run on a cluster. To do this, please
  * Create a results folder in your working directory to store output
