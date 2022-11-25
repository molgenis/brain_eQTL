#!/bin/bash
#SBATCH --job-name=DECON_TEST
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/DECON_TEST.out
#SBATCH --error=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/DECON_TEST.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.6.3-foss-2015b
source $HOME/venv/bin/activate

TMM="TMM:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-01-31-step2-filter-samples-from-alignmentQC/ROSMAP.geneCounts.2020-01-06.TMM.freeze2dot1.columnNamesFixed_samples_removed.txt.gz"
TMM_LOG2="TMM_LOG2:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz"
TMM_LOG2_GENECENTERED="TMM_LOG2_GENECENTERED:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.txt.gz"
TMM_LOG2_GENECENTERED_SAMPLEZTRANSFORMED="TMM_LOG2_GENECENTERED_SAMPLEZTRANSFORMED:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"
DATA_FILES=($TMM $TMM_LOG2 $TMM_LOG2_GENECENTERED $TMM_LOG2_GENECENTERED_SAMPLEZTRANSFORMED)

SINGLE_CELL="SINGLE_CELL:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/AMP-AD/single_cell_counts.txt.gz"
IHC="IHC:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/AMP-AD/IHC_counts.txt.gz"
CELLMAP_NOLOG2="CELLMAP_NOLOG2:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-14-CellMapPredictions/nolog2/CellMap_noLog2_counts.txt.gz"
CELLMAP_LOG2="CELLMAP_LOG2:/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-14-CellMapPredictions/log2/CellMap_log2_counts.txt.gz"
GT_FILES=($SINGLE_CELL $IHC $CELLMAP_LOG2 $CELLMAP_NOLOG2)

for DATA_FILE in ${DATA_FILES[@]}; do
    DATA_TYPE="${DATA_FILE%%:*}"
    DATA_FP="${DATA_FILE##*:}"
    echo $DATA_TYPE
    echo $DATA_FP
    for GT in ${GT_FILES[@]}; do
      GT_TYPE="${GT%%:*}"
      GT_FP="${GT##*:}"

      echo /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/partial_deconvolution.py -d $DATA_FP -si /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/all/links2-ExpressionSamplesLinked.txt -c AMP-AD -g $GT_FP -zscore -sum_to_one -o AMPAD_${DATA_TYPE}_${GT_TYPE}_zscore -visualise -e pdf
      python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/partial_deconvolution.py -d $DATA_FP -si /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/all/links2-ExpressionSamplesLinked.txt -c AMP-AD -g $GT_FP -zscore -sum_to_one -o AMPAD_${DATA_TYPE}_${GT_TYPE}_zscore -visualise -e pdf
      echo ""

      echo /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/partial_deconvolution.py -d $DATA_FP -si /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/all/links2-ExpressionSamplesLinked.txt -c AMP-AD -g $GT_FP -log2 -m 5 -sum_to_one -o AMPAD_${DATA_TYPE}_${GT_TYPE}_log2 -visualise -e pdf
      python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/partial_deconvolution.py -d $DATA_FP -si /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/all/links2-ExpressionSamplesLinked.txt -c AMP-AD -g $GT_FP -log2 -m 5 -sum_to_one -o AMPAD_${DATA_TYPE}_${GT_TYPE}_log2 -visualise -e pdf
      echo ""

      echo /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/partial_deconvolution.py -d $DATA_FP -si /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/all/links2-ExpressionSamplesLinked.txt -c AMP-AD -g $GT_FP -log2 -m 5 -zscore -sum_to_one -o AMPAD_${DATA_TYPE}_${GT_TYPE}_zscoreAndLog2 -visualise -e pdf
      python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/partial_deconvolution.py -d $DATA_FP -si /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-25-samplelinks/all/links2-ExpressionSamplesLinked.txt -c AMP-AD -g $GT_FP -log2 -m 5 -zscore -sum_to_one -o AMPAD_${DATA_TYPE}_${GT_TYPE}_zscoreAndLog2 -visualise -e pdf
      echo ""
    done
done

deactivate