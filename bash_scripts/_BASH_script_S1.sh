#!/bin/bash

##### GIVE SBATCH DIRECTIVES #####

#SBATCH --job-name=BASH_script_S1
#SBATCH --partition=week
#SBATCH --time=2-
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=ALL

##### SET UP JOB ENVIRONMENT #####

# Load standard environment:
module load StdEnv

# Specify environment variables to export:
export SLURM_EXPORT_ENV=ALL

# Load Rstudio:
module purge # clears environment of any previously loaded modules
module load R/4.2.0-foss-2020b
module list # checks to see that the correct modules are loaded

##### EXECUTE THE BATCH JOB #####

# Set working directory (in which R script is located):
cd /R_scripts

# Run the correct R script:
Rscript script_S1.R
