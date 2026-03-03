#!/bin/bash

##### GIVE SBATCH DIRECTIVES #####
#SBATCH --array 0-999
#SBATCH --job-name=BASH_script_S3
#SBATCH --mem-per-cpu 4g 
#SBATCH --partition=day
#SBATCH --time=12:00:00 
#SBATCH --mail-type ALL

##### SET UP JOB ENVIRONMENT #####

# Load standard environment:
module load StdEnv

# Specify environment variables to export:
export SLURM_EXPORT_ENV=ALL

# Load dSQ & Rstudio:
module purge # clears environment of any previously loaded modules
module load dSQ ### IS THIS IN THE RIGHT PLACE?
module load R/4.3.0-foss-2020b
module list # checks to see that the correct modules are loaded

##### EXECUTE THE BATCH JOB #####

# Set working directory (in which job list and R script are located):
cd /R_scripts

# Run the correct job array: (NOTE that the file paths below may need to change for the script to run on your cluster!)
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /R_scripts/script_S3_joblist.txt --status-dir /R_scripts


