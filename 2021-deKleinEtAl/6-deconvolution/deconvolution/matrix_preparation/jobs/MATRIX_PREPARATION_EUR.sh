#!/bin/bash
#SBATCH --job-name=MATRIX_PREPARATION_EUR
#SBATCH --output=/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/jobs/MATRIX_PREPARATION_EUR.out
#SBATCH --error=/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/jobs/MATRIX_PREPARATION_EUR.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.7.4-GCCcore-7.3.0-bare
source $HOME/env/bin/activate

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation.py -s CortexEUR-cis -n CortexEUR-cis

deactivate