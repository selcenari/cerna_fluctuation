#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --time=0-20:00:00
#SBATCH -J star_test
#SBATCH --partition=
# #SBATCH --constraint=
#SBATCH --account=
#SBATCH --ntasks=1


conda activate ceRNAnet
conda install -c bioconda bioconductor-sponge

# the --exclusive to srun makes srun use distinct CPUs for each job step
# -N1 -n1 allocates a single core to each task
#srun="srun --exclusive -N1 -n1"
srun="srun --exclusive -N1 -n1"

# --delay .2 prevents overloading the controlling node
# -j is the number of tasks parallel runs so we set it to $SLURM_NTASKS
# --joblog makes parallel create a log of tasks that it has already run
# --resume makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off
#today=$(date +%y%m%d-%H%M)

srun Rscript sponge_starbase.R