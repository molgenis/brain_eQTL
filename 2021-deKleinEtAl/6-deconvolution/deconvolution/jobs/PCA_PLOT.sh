#!/bin/bash
#SBATCH --job-name=PCA_PLOT
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/PCA_PLOT.out
#SBATCH --error=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/PCA_PLOT.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.6.3-foss-2015b
source $HOME/venv/bin/activate

python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/covariate_cluster_plot.py

deactivate