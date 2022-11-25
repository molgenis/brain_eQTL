# Partial Deconvolution Manual
This file describes how to perform partial deconvolution of bulk RNA-seq expression data using the code in this directory.  


## Introduction

This program performs partial deconvolution of a bulk RNA-seq expression matrix using a single cell reference profile. In order to speed up to process of multiple iterations, the program looks in the **-o** / **--outdir** folder if the data input options have been applied earlier. It compares the settings in the settings.json file and skips the [data_preprocessing](../partial_deconvolution/src/data_preprocessor.py) step if they are identical. This allows for fast iteration of different normalization inputs without the need of re-filtering the input files.

## Prerequisites  

This program is developed in Pycharm 2019.3 (Professional Edition), performance on other system is not guaranteed.

The program requires the following packages to be installed:  

 * pandas ([v1.0.1](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))  
 * numpy ([v1.18.3](https://pypi.org/project/numpy/#history); [BSD License](https://www.numpy.org/license.html))  
 * scipy ([v1.4.1](https://docs.scipy.org/doc/scipy/reference/release.html); [BSD License](https://www.scipy.org/scipylib/license.html)) 
 * matplotlib ([v3.2.1](https://github.com/matplotlib/matplotlib/releases); [BSD License](https://matplotlib.org/3.1.3/users/license.html))  
 * seaborn ([v0.10.0](https://github.com/mwaskom/seaborn); [BSD 3-Clause License](https://github.com/mwaskom/seaborn/blob/master/LICENSE))  

See 'Installing' on how to install these packages.

## Installing  

Installing and running the program can be done by executing the following steps:

**Step 1: acquire the source files**      
Either [clone](https://github.com/mvochteloo/brain_eQTL.git) or [download](https://github.com/mvochteloo/brain_eQTL/archive/master.zip) the source files.

**Step 2: creating the virtual environment [OPTIONAL]**    
1) Open JetBrains PyCharm  
2) Click on 'File' -> 'Settings' (Ctrl + Alt + S) -> 'Project: progr1' -> 'Project Interpreter'  
3) Add a new interperter by pressing the gear icon -> 'Add'  
4) Make sure 'New environment' is selected, click on 'OK' and again click on 'OK'  
5) Activate the virtual environment by either restarting PyCharm or executing the following command in the terminal:  
```console
source venv/bin/activate
```  

**Step 3: installing the required packages**   
  
1) Make sure the virtual environment is active and run the command  
```console  
pip install -r requirements.txt
```  
2) Wait for the command to finish and Pycharm has updated the indices  


## Command line arguments  
   
Syntax:
```console  
python3 ./partial_deconvolution.py
```  
Data input options:  
   
 * **-d** / **--data**: The bulk expression matrix. Rows are genes and columns are samples.
 * **-si** / **--signature**: The signature matrix. Rows are genes and columns are cell types.
 * **-t** / **--translate**: The gene to ensembl ID translate matrix. Must contain the columns 'ArrayAddress' and 'Symbol'.
 * **-sa** / **--sample**: The sample to cohort translate matrix. Must contain the columns 'RnaID' and 'GenotypeID'.
 * **-c** / **--cohort**: The cohort. Default: 'All'. You can use this setting to only analyse, for example, AMP-AD samples.
 * **-g** / **--ground_truth**: A matrix with the ground truth values. Must be in <type>_counts.<extension> format. Default: None'. Rows are samples and columns are cell types.
 
Data output options:
  
 * **-o** / **--outdir**: The name of the output directory. Default: 'output'.       

Signature matrix transformation options: 
   
 * **-m** / **--min_expr**: The minimal expression value per gene. Default: 0. All genes that have a maximum expression value lower than this value will be discarded. 
 * **-n** / **--normalize**: Divide each value by row row / column sum. Default: False.
 * **-zscore**: Z-score transform the profile values. Default: False.
 * **-log2**: Log2 transform the profile values. Default: False.
 * **-sum_to_one**: Sum-to-one the deconvolution weights. Default: False.

Deconvolution options: 
   
 * **-dm** / **--decon_method**: The deconvolution method to use. Default: 'NNLS'. Currently this is the only option.
 
Visualisation options:  
  
 * **-visualise**: Whether or not to visualise the data. Default: False.
 * **-e** / **--extension**: The figure file extension. Default: 'png'.     

## Author  

Martijn Vochteloo *(1)*

1. Institute for Life Science & Technology, Hanze University of Applied Science, Groningen, The Netherlands.

## License  

This project is licensed under the GNU GPL v3 license - see the [LICENSE.md](LICENSE.md) file for details
