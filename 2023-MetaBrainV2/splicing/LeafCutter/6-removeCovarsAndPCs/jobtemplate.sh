#!/bin/bash
#SBATCH --job-name=TISSUE
#SBATCH --output=/groups/umcg-biogen/tmp02/output/2021-FreezeThree/2021-02-18-splicing/2023-02-11-all-samples-leafcutter-fractionfix/6-removeCovarsAndPCs/outputLowQual-run2/TISSUE.log
#SBATCH --error=/groups/umcg-biogen/tmp02/output/2021-FreezeThree/2021-02-18-splicing/2023-02-11-all-samples-leafcutter-fractionfix/6-removeCovarsAndPCs/outputLowQual-run2/TISSUE.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Java/11.0.16
ml Python/3.10.4-GCCcore-11.3.0

BASE_DIR="/groups/umcg-biogen/tmp02/output/2021-FreezeThree/2021-02-18-splicing/2023-02-11-all-samples-leafcutter-fractionfix/"
tissue=TISSUE
dataset=DATASET

# regress covariates
java -Xmx18g -jar /groups/umcg-biogen/tmp02/output/2021-FreezeThree/2021-02-18-splicing/Regressor.jar \
	${BASE_DIR}/5-removeOutlierSamples/outputLowQual-run2/${tissue}-PSI-filtered.txt.gz \
	${BASE_DIR}../2021-08-27-covariateFiles/2022-11-28-${dataset}-covars-MDS.txt \
	${BASE_DIR}/6-removeCovarsAndPCs/outputLowQual-run2/${tissue}-PSI-filtered-logit-outliersRemoved-covars

# PCA
python ${BASE_DIR}/6-removeCovarsAndPCs/pca.py \
	${BASE_DIR}/6-removeCovarsAndPCs/outputLowQual-run2/${tissue}-PSI-filtered-logit-outliersRemoved-covars.CovariatesRemovedOLS.txt.gz \
	${BASE_DIR}/6-removeCovarsAndPCs/outputLowQual-run2-pcs/${tissue}-PSI-filtered-logit-outliersRemoved-covars

# merge covariates
python ${BASE_DIR}/6-removeCovarsAndPCs/mergecovars.py \
	${BASE_DIR}../2021-08-27-covariateFiles/2022-11-28-${dataset}-covars-MDS.txt \
	${BASE_DIR}/6-removeCovarsAndPCs/outputLowQual-run2-pcs/${tissue}-PSI-filtered-logit-outliersRemoved-covars_PCs.txt \
	5 \
	${BASE_DIR}/6-removeCovarsAndPCs/outputLowQual-run2-pcs/${tissue}-covarsAnd5PCs.txt.gz

# regress PCA + covariates
java -Xmx18g -jar /groups/umcg-biogen/tmp02/output/2021-FreezeThree/2021-02-18-splicing/Regressor.jar \
	${BASE_DIR}/5-removeOutlierSamples/outputLowQual-run2/${tissue}-PSI-filtered.txt.gz \
	${BASE_DIR}/6-removeCovarsAndPCs/outputLowQual-run2-pcs/${tissue}-covarsAnd5PCs.txt.gz \
	${BASE_DIR}/6-removeCovarsAndPCs/outputLowQual-run2/${tissue}-PSI-filtered-logit-outliersRemoved-covarsAnd5PCs
