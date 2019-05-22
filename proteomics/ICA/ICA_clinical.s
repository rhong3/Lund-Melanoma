#!/bin/bash
#
#SBATCH --job-name=ICA_clinical
#SBATCH --partition=c01_17 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=07-00:00
#SBATCH --mem=100GB 
#SBATCH --mail-type=END                                                        
#SBATCH --mail-user=Runyu.Hong@nyumc.org
#SBATCH --mail-user=rh2740@nyu.edu 
#SBATCH --output=ICA_clinical_%j.out
 
module load r/intel/3.3.2

SRCDIR=/beegfs/rh2740/Fenyo/
RUNDIR=$SCRATCH/run-${SLURM_JOB_ID/.*}
mkdir -p $RUNDIR

cd $SLURM_SUBMIT_DIR

## srun R CMD BATCH example.R example.out
Rscript clinical_association_template.R clean_proteomics ../wna_clinical.tsv Clean_proteomics_ICA.Rdata . 0.002 
exit
