#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=tn2220@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load anaconda3/2020.07
module load nextflow/21.10.6


DATADIR=$SCRATCH/Mitonuclear/TEST1

cd $DATADIR
nextflow run main.nf

