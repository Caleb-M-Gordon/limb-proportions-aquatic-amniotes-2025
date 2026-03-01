#!/bin/bash
#SBATCH --output /home/cmg89/palmer_scratch/output_files/script_S4_array_%A/%N/slurm-%A_%a.out
#SBATCH --array 0-9215
#SBATCH --job-name dsq-script_S4_joblist
#SBATCH --mem-per-cpu 4g -p day -t 1- --mail-type ALL

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/project/bhullar/cmg89/Flipper_Project/script_S4_joblist.txt --status-dir /gpfs/gibbs/project/bhullar/cmg89/Flipper_Project

