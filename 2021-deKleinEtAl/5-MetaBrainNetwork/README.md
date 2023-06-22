
Creating co-regulation network
=======

See also https://github.com/npklein/Genenetwork-cookbook

These are scripts to help create metabrainnetwork.nl. Run the three steps in order:

bash step1-PCA-outlier-detection.sh
bash step2-correlation-and-evd.sh
bash step3-backend-matrix-creation.sh 

All these scripts use lmod module system to load required software, if you run this on a system without lmod you will have to delete 

# Prerequisites:

Download RunV13.jar
wget https://github.com/npklein/metabrain/releases/download/2020-03-02/RunV13.jar


Download eQTL mapping pipeline
wget https://molgenis50.gcc.rug.nl/jenkins/job/systemsgenetics/473/nl.systemsgenetics$GeneNetworkBackend/artifact/nl.systemsgenetics/GeneNetworkBackend/1.0.7-SNAPSHOT/GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar

# step1-PCA-outlier-detection.sh

This step you want to possibly run multipe times. Among other things, it will output a PCA plot. Check this PC plot for outliers. 
If there are any, remove these from the file you use for the -s samples_to_include option, and rerun this step (with the raw gene count file, and change the -q option). 
If there are no outliers, go to step2-correlation-and-evd.sh.

## options
Mandatory options: 

```
  -t      TMPDIR where files will be written during runtime
  -e      Raw gene count file with columns = samples, rows = genes. For metabrain.nl kallisto was used, but other counts should work also
  -s      File with list of samples to include (exclude all others)
  -p      Base of the project_dir where config files will be written
  -o      Output directory where results will be written
  -j      Location of RunV13.jar
  -g      Location of https://github.com/npklein/metabrain/tree/master/2021-deKleinEtAl/5-MetaBrainNetwork directory on your machine
  -a      GTF file
  -m      Memory to give to the eqtlgen jar file when running, will be appendended to -Xmx, so 8g will be java -jar -Xmx8g. For 8000 samples I used 100g
  ```
  
  Optional options:
  
```
  -q      Current step, e.g. 1a or 1b. This way can distinguish between multiple PCA rounds
  -h      display help
```

## output

All output is written to your -o output directory.

0_remove_genes/
  Expression file with genes removed that have exact same counts as other genes (likely duplicates), are only found on scaffolds, or have duplicate names
1_selectSamples/
  Expression file from step 0, but with only samples included in the -s sample list and genes with no variance are removed
2_removeDuplicates
  Expression file from step 1, but with duplicate samples (with same counts) removed (used because we have samples from ENA, which might contain duplicates
3_quantileNormalized_{step}
  Expression file from step 2, but quantile normalized. The outdir {step} is what is given with the -q step option. 
  Also includes a PC plot. Check PC plot for outliers before going to next step
  
  
# step2-correlation-and-evd.sh

As with previous script, this script uses lmod module system. If your environment doesn't use this, remove lines starting with ml and install have these programs in your PATH.

This script submits jobs to a slurm cluster, and to do this uses commands that are specific to slurm. If you run this in a different cluster environment, adjust the following files:

  - scripts_per_step/6_CorrelationMatrix.sh
  - scripts_per_step/7_evd_on_correlation.sh

Also, in this script, change the lines

    jobid=$(sbatch --parsable 7_evd_on_correlation.sh)
    
    and
    
    while [ ! $(squeue -j $jobid 2> /dev/null | wc -l) -le 1 ];

To what your cluster uses to submit jobs (e.g. bsub 7_evd_on_correlation.sh). First line has to output a job id, second line has to output how many lines are still in the scheduler with jobs with that job id.


## options


Mandatory options:

```
  -t      TMPDIR where files will be written during runtime
  -e      Raw gene count file with columns = samples, rows = genes. Same as used for input step 1
  -p      Base of the project_dir where config files will be written
  -o      Output directory where results will be written
  -j      Location of eqtl-mapping-pipeline.jar
  -s      File with samples to include (all others will be excluded)
  -g      Location of https://github.com/npklein/metabrain/tree/master/2021-deKleinEtAl/5-MetaBrainNetwork directory on your machine
  -z      Covariate table. Columns are samples, rows are covariates. 
  -a      GTF file
  -v      Number of threads to use for correlation step and PCA step
  -m      Memory to use for some steps
  -q      qos to run sbatch jobs in
  -n      Name that will be used in output file
  -y      the previous step which has the correctly selected samples. e.g. step1b, will use this to find right input directory
```

Optional options:

```
-h      display help
```

## output
All output is written to your -o output directory.

4_deseqNormalized/
  Expression table normalized by median of ratio's as used and described by DESeq
5_covariatesRemoved/
  Expression table from step 4, but covariates regressed out
6_correlation_matrix/
  gene-gene correlation matrix
7_evd_on_correlation_matrix/
  evd on the correlation matrix (similar to doing PCA directly on matrix from 5_covariatesRemoved/


## Covariates
The covariates we used for metabrain.nl (from FastQC, STAR aligner, and GATK MultipleMetics). 
We identified these by correlating covariates to our expression data and taking the 20.

```
PCT_CODING_BASES
PCT_MRNA_BASES
PCT_INTRONIC_BASES
MEDIAN_3PRIME_BIAS
PCT_USABLE_BASES
PCT_INTERGENIC_BASES
PCT_UTR_BASES
PF_HQ_ALIGNED_READS
PCT_READS_ALIGNED_IN_PAIRS
PCT_CHIMERAS
PF_READS_IMPROPER_PAIRS
PF_HQ_ALIGNED_Q20_BASES
PF_HQ_ALIGNED_BASES
PCT_PF_READS_IMPROPER_PAIRS
PF_READS_ALIGNED
avg_mapped_read_length
avg_input_read_length
uniquely_mapped
total_reads
Total.Sequences_R1
```

# step3-backend-matrix-creation.sh
As with previous script, this script uses lmod module system. If your environment doesn't use this, remove lines starting with ml and install have these programs in your PATH.

This script submits jobs to a slurm cluster, and to do this uses commands that are specific to slurm. If you run this in a different cluster environment, adjust the following files:

  - scripts_per_step/9_CorrelateEigenvectors.sh
  - scripts_per_step/7_evd_on_correlation.sh

Also, in this script, change the lines

    job_ids+=("$(sbatch --parsable ${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh)")    
    
    and
    
    while [ ! $(squeue -j $jobid 2>/dev/null | wc -l ) -le 1 ];

To what your cluster uses to submit jobs (e.g. bsub 7_evd_on_correlation.sh). First line has to output a job id, second line has to output how many lines are still in the scheduler with jobs with that job id.



## options

Mandatory options
```
  -t      TMPDIR where files will be written during runtime
  -p      Base of the project_dir where config files will be written
  -o      Output directory where results will be written
  -j      Location of eqtl-mapping-pipeline.jar
  -g      Location of https://github.com/npklein/metabrain/tree/master/2021-deKleinEtAl/5-MetaBrainNetwork directory on your machine
  -a      GTF file
  -v      Number of threads to use for correlation step and PCA step
  -d      GeneNetwork directory (with backend data for predictions
  -m      Memory to use for some steps (that are submitted to the cluster)
  -n      Name that will be used in output file. !!!Needs to be same as for the step2-correlation-and-evd.sh script!!!
  -z      Number of eigenvectors to use
  -q      qos to use when submitting to cluster
  -e      file with 1st column ensg ID, 3rd column HNGC symbol
  -h      display help
```

Optional options

```
  -h      display help
```

## output
8_CenterScaledColumnsRows/
  Center scaled columns and rows of the evd matrix from 7_evd_on_correlation_matrix
9_eigenvector_correlation_matrix/
  Correlation matrix of the 8_CenterScaledColumnsRows/ matrix
10_GeneNetwork_predictions/
  GeneNetwork predictions
11_ImportToWebsite
  Files necesarry to create metabrain.nl
  
# Creating the website
