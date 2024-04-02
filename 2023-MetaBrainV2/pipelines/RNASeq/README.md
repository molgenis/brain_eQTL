# RNA-Seq Quantification Pipeline
*Interns: Orfeas Gkourlias & Joost Bakker*<br>
*Suprevisor: Harm-Jan Westra* <br>
*UMCG Group: Functional Genomics*<br>
o.gkourlias@umcg.nl <br>
j.bakker@umcg.nl <br>
h.j.westra@umcg.nl

This repository consists of a newly developed Nextflow pipeline. The pipeline
is based on previous works created by Niek, also under suprevision of Harm-Jan.
The pipeline consists of several steps, known as processes in Nextflow, which come together to perform alignment and quantification on RNA-Seq samples.

The pipeline consists of steps, which are called processes in Nextflow. These processes have distinct functionalities and all of them handle the input or proceses files in some way. A step by step description can be found below.

1. convertBAMToFASTQ
   - Input: A singular unsorted .BAM file or one (two if paired) .faq file(s)
   - Output: A FastQ file.
   - Function: Handles .BAM to .faq conversion in case the input is in .BAM format. Always uses SAMtools sort before conversion.
2. fastqcQualityControl
   - Input: .FastQ
   - Output: Zipped .FastQ file.
   - Function: Takes FastQ file and performs FastQC on it. It then returns it.
3. alignWithSTAR
   - Input: .FastQ file.
   - Output: .Bam file & Alignment log files.
   - Function: Align the input FastQ file against the reference panel.
4. sortBAM
   - Input: .BAM file.
   - Output: Sorted .Bam file.
   - Function: Sorted the 
5. markDuplicates
   - Input: .BAM file.
   - Output: Filtered .BAM file & marked duplicate log files.
6. QCwithRNASeqMetrics
   - Get random reads in both forward and reverse orientation for the polishing step. Amount can be set in the config.
7. QCwithMultipleMetrics
   - Returns the variants in reads compared to the raw centroid draft.
8. identifyAlternativeSplicingSitesrMATS
   - Create the consensus sequence from the reference sequence and variants.
9. identifyAlternativeSplicingSitesLeafCutter
   - The consensus sequence is compared to different types of vaccins. If the percentage identity is 98% or higher, the consensus sequence is typed as this vaccine stem.  
10. convertBAMToCRAM
    - Identify and censor samples with low coverage. Returns 'N's as consensus sequence if the assigned reads is lower than the threshold.

## System requirements

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
apptainer pull docker://ogkourlias/rnaseq_metabrain:latest
```

A download process will start. After the process is finished a new file called `popassign_v0.7.sif` will be created.

## Input
The pipeline expects certain required inputs, and there are some optional inputs. These inputs should be configured in `nextflow.config` or can be passed as command line parameters.

`--refFlat` Path to the refflat file 

`--referenceGenome` Path to a reference genome FASTA file

`--gtfAnnotationFile` Path to the gene annotation GTF file

`--sampleFile` Text file containing a list of input BAM or FASTQ samples. This file can be created with the instructions in the Usage section of this README.

`--outDir` Path to the desired output directory

`--refDir` Path to the directory containing the STAR reference

`--bindFolder` Folder that the apptainer container should bind to.

`--ribosomalIntervalList` Path to the ribosomal interval list

## Usage
### 1. Create sample file
To create an input sample file, you can use the `index_dir.py` script. This script expects the following parameters:

`--input_dir` The directory in which the sample BAM or FASTQ files are stored (REQUIRED)

`--out` Desired prefix of the output file (REQUIRED)

After running the script, fill in the path to the resulting sample file in the `nextflow.config` file. 

### 2. Run the pipeline
After the installation and configuration as described in the previous parts of this README, you can simply run the pipeline with the following command:

```
nextflow run main.nf
```

Numerous things can cause the pipeline to crash. If this happens, the pipeline can easily be resumed since Nextflow automatically caches progress of all processes. To resume the pipeline, just extend the run command with the `-resume` flag:

```
nextflow run main.nf -resume
```

## Output

After succesful completion of the pipeline, the output directory has the following content:

```
output_dir
    ├── cram
    ├── fastqc
    ├── leafcutter
    ├── mark_duplicates
    ├── multiple_metrics
    ├── rmats
    ├── rna_seq_metrics
    ├── star
```

Every folder in the above directory tree contains output of the corresponding tool.

## Citations

[P. Danecek et al., “Twelve years of SAMtools and BCFtools,” Gigascience, vol. 10, no. 2, pp. 1–4, Jan. 2021, doi: 10.1093/GIGASCIENCE/GIAB008.](https://pubmed.ncbi.nlm.nih.gov/33590861/)

[Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21 (2013).](https://pubmed.ncbi.nlm.nih.gov/23104886/)

[Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data.](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Picard Tools - By Broad Institute.](http://broadinstitute.github.io/picard/)

[Wang, Y. et al. rMATS-turbo: an efficient and flexible computational tool for alternative splicing analysis of large-scale RNA-seq data. Nature Protocols 2024 1–22 (2024) doi:10.1038/s41596-023-00944-2.](https://www.nature.com/articles/s41596-023-00944-2)

[Li, Y. I. et al. Annotation-free quantification of RNA splicing using LeafCutter. Nature Genetics 2017 50:1 50, 151–158 (2017).](https://www.nature.com/articles/s41588-017-0004-9)

## Contact
Orfeas Gkourlias: o.gkourlias@umcg.nl <br>
Joost Bakker: j.bakker@umcg.nl <br>
Harm-Jan Westra: h.j.westra@umcg.nl <br>

## License
The code in this direcotory and its subdirectories is licensed under the BSD 3-Clause License. See the `LICENSE` file in this directory for more information. The license applies only to the current directory and all its subdirectories; it does not extend to higher directories.
