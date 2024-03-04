# MetaBrain Whole Genome Sequencing QC pipeline
This pipeline is developed for automated quality control of WGS data. The following steps are performed in this pipeline:

- General filtering (minor allele frequency, call-rate, missingness, Hardy-Weinberg P)
- Sample filtering based on a sex check
- Sample filtering based on a heterozygosity check
- Sample filtering based on a relatedness check
- Sample filtering based on z-score outliers
- Population based filtering of outliers

## System requiremets
This pipeline is developed to run on an HPC with the following software installed:
- Nextflow (v23.10)
- Apptainer (v1.1.3)
- Java (v11 or higher)
- Bash (v3.2 or higher)

## Installation
Please follow the steps below to set up the pipeline.

### 1. Clone the repository
Run the following command to clone the repository:

```
git clone https://github.com/molgenis/metabrain.git
```

### 2. Pull the container
This pipeline comes with an Apptainer container that has all required software installed to run the pipeline. First, create a cache directory and set it as `APPTAINER_CACHEDIR`:

```
mkdir apptainer_cache
export APPTAINER_CACHEDIR=apptainer_cache
```

In the code above, `apptainer_cache` can be replaced with any desired name. Then, use the following command to pull the container:

```
apptainer pull docker://quay.io/eqtlgen/popassign:v0.7
```

A download process will start. After the process is finished a new file called `popassign_v0.7.sif` will be created.

### 3. Configure parameters
In `nextflow.config`, set the `process.container` parameter with the path to the image file created in the previous step, and set the `singularity.cacheDir` parameter with the path to the cache directory, also created in the previous step. Configure the remaining parameters as described in the next section.

## Input
The pipeline expects certain required inputs, and there are some optional inputs. These inputs should be configured in `nextflow.config` or can be passed as command line parameters.
### Required inputs

The following parameters are required:

`--bindFolder` Folder that the singularity container should bind to.

`--sampleFile` Text file containing a list of input VCF samples. This file can be created with the instructions in the Usage section of this README.

`--outDir` Path to the desired output directory

`--refPath` Path to the reference panel in plink format (bed/bim/bam) without file extensions

`--refPop` Path to a text file containing information about the populations of the 1000 Genomes reference panel in the following format:
  `SampleID	Population	Superpopulation`

`--refAfs` Path to the gzipped allele frequencies of the reference panel

`--sexFile` Path to a file containing sex information about the samples in the following format: `SampleId,F` (female) or `SampleId,M` (male)

### Optional inputs

The following are optional:

`--maf` Threshold for minor allele frequency filtering (default: 0.01)

`--geno` Threshold for maximum missing genotype rate (default: 0.05)

`--mind` Threshold for the maximum proportion of missing genotypes for individuals (default: 0.05)

`--hwe` Threshold for Hardy-Weinberg equilibrium (default: 1e-6)

`--kingTableFilter` Threshold for kinship in the relatedness check (default: 2^-4.5)

`--populationOutlierThreshold` Threshold for filtering population outliers (default: 0.4)

## Usage
### 1. Create sample file
To create an input sample file, you can use the `index_dir.py` script. This script expects the following parameters:

`--input_dir` The directory in which the sample files are stored (REQUIRED)

`--out` Desired prefix of the output file (REQUIRED)

`--file_substring` Only files containing this substring are indexed (OPTIONAL)

### 2. Run the pipeline
After the installation and configuration as described in the previous parts of this README, you can simply run the pipeline with the following command:

```
nextflow run main.nf
```

Numerous things can cause the pipeline to crash. If this happens, the pipeline can easily be resumed since Nextflow automatically caches progress of all processes. To resume the pipeline, just extend the run command with the `-resume` flag:

```
nextflow run main.nf -resume
```

## Acknowledgements
The code for this pipeline is built upon the work of Urmo Võsa, Robert Warmerdam, Harm-Jan Westra, and Martijn Vochteloo, to whom we want to express our sincere gratitude.

### Citations
[P. DI Tommaso, M. Chatzou, E. W. Floden, P. P. Barja, E. Palumbo, and C. Notredame, “Nextflow enables reproducible computational workflows,” Nature Biotechnology 2017 35:4, vol. 35, no. 4, pp. 316–319, Apr. 2017, doi: 10.1038/nbt.3820.](https://www.nature.com/articles/nbt.3820)

[P. Danecek et al., “Twelve years of SAMtools and BCFtools,” Gigascience, vol. 10, no. 2, pp. 1–4, Jan. 2021, doi: 10.1093/GIGASCIENCE/GIAB008.](https://pubmed.ncbi.nlm.nih.gov/33590861/)

[S. Purcell et al., “PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses,” Am J Hum Genet, vol. 81, no. 3, p. 559, 2007, doi: 10.1086/519795.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/)

[F. Prive, H. Aschard, A. Ziyatdinov, and M. G. B. Blum, “Efficient analysis of large-scale genome-wide data with two R packages: bigstatsr and bigsnpr,” Bioinformatics, vol. 34, no. 16, pp. 2781–2787, Aug. 2018, doi: 10.1093/BIOINFORMATICS/BTY185.](https://pubmed.ncbi.nlm.nih.gov/29617937/)

[A. Auton et al., “A global reference for human genetic variation,” Nature 2015 526:7571, vol. 526, no. 7571, pp. 68–74, Sep. 2015, doi: 10.1038/nature15393.](https://www.nature.com/articles/nature15393)

## Contact
For any questions or comments about this pipeline, contact: j.m.bakker@umcg.nl
