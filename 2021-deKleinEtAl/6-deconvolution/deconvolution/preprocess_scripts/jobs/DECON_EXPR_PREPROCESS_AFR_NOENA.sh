#!/bin/bash
#SBATCH --job-name=DECON_EXPR_PREPROCESS_AFR_NOENA
#SBATCH --output=/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/jobs/DECON_EXPR_PREPROCESS_AFR_NOENA.out
#SBATCH --error=/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/jobs/DECON_EXPR_PREPROCESS_AFR_NOENA.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.7.4-GCCcore-7.3.0-bare
source $HOME/env/bin/activate

/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.txt.gz -t /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-05-step6-covariate-removal/2020-05-26-step5-remove-covariates-per-dataset/2020-05-25-covariatefiles/2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz -gte /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR -p GTE-AFR- -e ENA -od /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts -of CortexAFR-cis-NoENA



deactivate