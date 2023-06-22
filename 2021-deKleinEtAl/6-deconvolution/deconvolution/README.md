# Deconvolution  

## Introduction
This directory contains code written for the deconvolution of the MetaBrain bulk RNA-seq dataset.

## Prerequisites  

This program is developed in Pycharm 2020.3.3 (Professional Edition), performance on other system is not guaranteed.

Code in this directory uses the following packages:  

 * [adjustText](https://pypi.org/project/adjustText/) (v0.7.3; UNKNOWN)
 * [colour](https://pypi.org/project/colour/) (v0.1.5; BSD (3-clause))
 * [functools](https://pypi.org/project/functools/) (v0.5; BSD) 
 * [matplotlib](https://pypi.org/project/matplotlib/) (v3.3.4; PSF)  
 * [numpy](https://pypi.org/project/numpy/) (v1.19.5; BSD)  
 * [openpyxl](https://pypi.org/project/openpyxl/) (v3.0.6; MIT)
 * [pandas](https://pypi.org/project/pandas/) (v1.2.1; BSD (3-clause))
 * [requests](https://pypi.org/project/requests/) (v2.25.1; Apache 2.0)
 * [scikit-learn](https://pypi.org/project/scikit-learn/) (v0.24.1; new BSD)
 * [scipy](https://pypi.org/project/scipy/) (v1.6.0; BSD)
 * [seaborn](https://pypi.org/project/seaborn/) (v0.11.1; BSD (3-clause))
 * [sklearn](https://pypi.org/project/sklearn/) (v0.0; UNKNOWN)
 * [statsmodels](https://pypi.org/project/statsmodels/) (v0.12.1; BSD License)    
 * [sympy](https://pypi.org/project/sympy/) (v1.7.1; BSD)
 * [UpSetPlot](https://pypi.org/project/UpSetPlot/) (v0.4.1; BSD (3-clause))  
  

Install the packages using:
```console  
pip install -r requirements.txt
```  

## Contents
  
 * **custom_interaction_analyser/** code used to test for interacting eQTLs using an f-test.  
 * **matrix_preparation/** code used to pre-process the MetaBrain data. Conists of several steps including reordering correcting for cohort effects, filtering and reorder of the input data, predict cell type proportions, etc.  
 * **paper_scripts/** code used in the process of writing the manuscript.  
 * **partial_deconvolution/** code to perform partial deconvolution om the MetaBrain bulk RNA-seq dataset. Allows the dynamically change a lot of settings/inputs/comparisons in order to evaluate performance.     
 * **presentation_scripts/**  code for creating a plot for presentation purposes.  
 * **r_scripts/** code written in R.   
 * **sn_scripts/** code related to the single-nucleus data analysis.   
 * **test_scripts/** divers code made to perform a single task (e.g. test/plot something).  
  
## Author  

Martijn Vochteloo *(1)*

1. University Medical Centre Groningen, Groningen, The Netherlands

## License  

This project is licensed under the GNU GPL v3 license - see the [LICENSE.md](LICENSE.md) file for details
