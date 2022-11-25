#!/bin/bash
#SBATCH --job-name=DIM_REDUCTION
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/DIM_REDUCTION.out
#SBATCH --error=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/DIM_REDUCTION.out
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.6.3-foss-2015b
source $HOME/venv/bin/activate

python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/dim_reduction_plot.py -ex ../data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.pkl -m /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -mid rnaseq_id

deactivate