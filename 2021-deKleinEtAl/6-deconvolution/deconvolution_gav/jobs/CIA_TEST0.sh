#!/bin/bash
#SBATCH --job-name=CIA_TEST0
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-decon-optimzer/jobs/output/CIA_TEST0.out
#SBATCH --error=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-decon-optimzer/jobs/output/CIA_TEST0.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module load Python/3.6.3-foss-2015b
source $HOME/venv/bin/activate

python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-decon-optimzer/custom_interaction_analyser.py -n cotex_eur_cis -s standard_settings -ne 50 -ns 2970

deactivate