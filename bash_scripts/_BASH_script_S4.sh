#!/bin/bash
##### GIVE SBATCH DIRECTIVES #####
#SBATCH --output /R_scripts/script_S4_array_%A/%N/slurm-%A_%a.out
#SBATCH --array 0-9215
#SBATCH --job-name dsq-script_S4_joblist
#SBATCH --mem-per-cpu 4g -p day -t 1- --mail-type ALL

# Set working directory (in which job list and R script are located):
cd /R_scripts

# Run the correct job array: (NOTE that the file paths below may need to change for the script to run on your cluster!)
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /R_scripts/script_S4_joblist.txt --status-dir /R_scripts
