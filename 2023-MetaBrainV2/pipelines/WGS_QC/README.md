# MetaBrain Whole Genome Sequencing QC pipeline

## System requiremets
This pipeline is developed to run on an HPC with the following software installed:
- Nextflow (v23.10)
- Apptainer (v1.1.3)
- Java (v11 or higher)
- Bash (v3.2 or higher)

## Input
The pipeline expects certain required inputs, and there are some optional inputs. These inputs should be configured in `nextflow.config` or can be passed as command line parameters.
### Required inputs

The following parameters are required:

`--bindFolder` Folder that the singularity container should bind to.

`--sampleFile` Text file containing a list of input VCF samples

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
