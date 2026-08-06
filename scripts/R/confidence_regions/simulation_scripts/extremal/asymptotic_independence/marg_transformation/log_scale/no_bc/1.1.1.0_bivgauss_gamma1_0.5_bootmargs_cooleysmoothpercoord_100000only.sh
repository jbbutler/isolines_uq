#!/bin/bash
#SBATCH --job-name bivgauss_sims_100000
#SBATCH --nodes 1
#SBATCH --output bivgauss_sims_100000.out
#SBATCH --error bivgauss_sims_100000.err
#SBATCH -p high
#SBATCH --cpus-per-task 24
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user butlerj@berkeley.edu

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate isolinesR
export OMP_NUM_THREADS=1

DF_SAVE_PATH="~/isolines_uq/outputs/simulations/extremal/"

Rscript 1.1.1.0_bivgauss_gamma1_0.5_bootmargs_cooleysmoothpercoord_100000only.R --save_df_path "$DF_SAVE_PATH" --n_cores 24