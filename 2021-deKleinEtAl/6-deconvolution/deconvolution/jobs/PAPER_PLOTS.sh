#!/bin/bash
#SBATCH --job-name=PAPER_PLOTS
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/PAPER_PLOTS.out
#SBATCH --error=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/PAPER_PLOTS.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module load Python/3.6.3-foss-2015b
source $HOME/venv/bin/activate

python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n cis_new_output -s new_settings -a 0.00036 -p deconvolution_covariate_comparison -e png
python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n cis_new_output -s new_settings -a 0.00036 -p simple_eqtl_effect inter_eqtl_zscore_bars inter_eqtl_effect inter_eqtl_celltype_details -i 12 158 -e png
python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n ms_cis_new_output -s new_settings -a 0.000435 -p simple_eqtl_effect inter_eqtl_zscore_bars inter_eqtl_effect inter_eqtl_celltype_details -i 35 -e png

python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n cis_new_output -s new_settings -a 0.00036 -p deconvolution_covariate_comparison -e pdf
python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n cis_new_output -s new_settings -a 0.00036 -p simple_eqtl_effect inter_eqtl_zscore_bars inter_eqtl_effect inter_eqtl_celltype_details -i 12 158 -e pdf
python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/visualiser.py -n ms_cis_new_output -s new_settings -a 0.000435 -p simple_eqtl_effect inter_eqtl_zscore_bars inter_eqtl_effect inter_eqtl_celltype_details -i 35 -e pdf

deactivate