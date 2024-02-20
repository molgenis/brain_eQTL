# RNA-Seq Based Genotyping Pipeline.
This pipeline program is dedicated to the genotyping of RNA-Seq based STAR alignments using GATk4.4.0.0.
The pipeline uses Nextflow and is built for HPC cluster utilizing the SLURM scheduler.
By taking RNA-Seq based alignments, this pipeline is able to call variants and genotypes based on the alignment input.

## Installation
All of the required software except for Nextflow itself is present in the docker container. This container can either be built through useage of the Dockerfile or through Dockerhub. You may simply call upon the following command to pull it from Dockerhub.

```bash
docker pull ogkourlias/genotyping:latest
```
Once the docker container is installed and running, you are ready to use the nextflow script.

## Usage
Before using the pipeline, you will need to configure the parameters in the Nextflow.config file.
The following parameters are required.

- process.container: Path to the container you have built or pulled.
- cram_dir: Path to a directory containing all of the .cram format files you want to call genotypes for.
- referenceGenome: Path to the reference genome file (.fa format) that was used for the alignment of your input file(s) 
- knownSites: Path to a known sites file which contains known variant IDs, such as GCF_000001405.40.gz
- knownSitesIndex: Path to the index file for the knownSites file in TBI format (e.g. GCF_00001405.40.gz.tbi)
- outDir: Output directory for the final .vcf.gz files containg the variants and genotypes.
- n_thread: The amount of threads you wish to dedicate to the pipeline.

There are also singularity settings which may be interesting if you are using singularity on your HPC enviroment.
After assigning values to the parameters in Nextflow.config, make sure you have the Nextflow module installed and loaded on your machine.
You make now invoke the pipeline and start the run using the following command when the bin directory is your working directory.

```bash
nextflow run main.nf
```

Once this command has been invoked in your CLI, Nextflow will start creating processes.
The structure and run order of the pipeline processes is as follows:

1. cramToBam
2. indexBam
3. encodeConvert
4. splitNCigarReads
5. AddOrReplaceReadGroups
6. baseRecalibrator
7. applyBQSR
8. haplotypeCaller
9. indexGvcf
10. jointGenotype

After completing of the jointGenotype step, the output VCF files will be present in the user defined directory.
## License

[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)