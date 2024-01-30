# RNA-Seq Quantification Pipeline.
*Interns: Orfeas Gkourlias & Joost Bakker*<br>
*Suprevisor: Harm-Jan Westra* <br>
*UMCG Group: Functional Genomics*<br>
o.gkourlias@umcg.nl <br>
j.bakker@umcg.nl <br>
h.j.westra@umcg.nl

## Introduction
This repository consists of a newly developed Nextflow pipeline. The pipeline
is based on previous works created by Niek, also under suprevision of Harm-Jan.
The pipeline consists of several steps, known as processes in Nextflow, which come together to perform alignment and quantification on RNA-Seq samples.

## Requirements & Configuration
There are two hard software requirements to perform the pipeline. These are:  
1. The ability to build docker images.
2. A Unix based system with Nextflow installed.

It is highly recomended to run the pipeline
on a HPC cluster. 

There are input requirements and configuartions the user needs to tweak to their use case. The paths to all the files and directories above need to be in the nextflow.config file. The parameters present there may be changed.

- `sampleFile`: A text file containg paths to the input RNA-Seq sample files.
- `referenceGenome`: A reference genome, like Grch38: .fa.
- `gtfAnnotationFile`: A corresponding GTF file.
- `ref_dir`: Index files created with STAR version 3.7.11a. This is important. The config file expects a directory containing all index files.


## Usage
Once the parameters in the config file have been configured, the pipeline
is ready to run. This can be done in a couple of ways, depending on the system. The recommended approach is to run the pipeline using singularity/apptainer on a HPC cluster. To do this, use the following command, assuming the working directory is where main.nf is located:

```bash
Nextflow run --with-singularity rnaseq_metabrain.sif main.nf
```

Replace the "--with-singularity" option if applicable.

## Pipeline Description
The pipeline consists of steps, which are called processes in Nextflow. These processes have distinct functionalities and all of them handle the input or proceses files in some way. A step by step description can be found below in addition to a DAG.

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


![DAG](archive/dag.png)

## Contact
Orfeas Gkourlias: o.gkourlias@umcg.nl <br>
Joost Bakker: j.bakker@umcg.nl <br>
Harm-Jan Westra: h.j.westra@umcg.nl <br>
