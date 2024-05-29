#!/bin/bash
#SBATCH --job-name dtr_sims      # Set a name for your job. This is especially useful if you
#SBATCH --partition quick     # Slurm partition to use: quick, norm, 
#SBATCH -c 1        # Number of tasks to run 3
#SBATCH --ntasks-per-core=1	# Do not use hyperthreading (this flag typically used for parallel jobs)
#SBATCH --time 0-00:25       # Wall time limit in D-HH:MM
#SBATCH --mem 5000     # Memory limit for each tasks (in MB) # 1500
#SBATCH -o /home/error/dtr_sims.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e /home/error/dtr_sims.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --array 1-1000
#module load gcc/7.1.0-fasrc01 R/3.5.1-fasrc02
#export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER
#you can find bash in: /home/bash
module load R/4.2.2

Rscript '/home/dtr/code/dtr_msm_sims_dt.R' $1 $2
