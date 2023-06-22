#!/usr/bin/env bash

# Download and make quantification pipeline
# This script clones the Molgenis Compute pipeline for alignin FastQ files and
# modifies the protocols and parameter files for analysis of Brain eQTL data.
# Makes samplesheet dependent parameter files.

# NOTE: This script changes paths in a parameters.csv that are specific to our clusters
#       If you want to use this somewhere else you have to change that in the script"


set -e
set -u

module load Python

cohort=
project_dir=
tmpdir=
samplesheet=
RNAseqDir=
script_dir=$(dirname "$0")

main(){
    parse_commandline "$@"
    # move to project dir if it exists.
    # I prefer making the project dir outside of this script in case wrong directory is given
    if [ -d "$project_dir" ];
    then
        echo "Changing directory to $project_dir"
        cd $project_dir
    else
        echo "ERROR: $project_dir does not exist"
        exit 1;
    fi
    clone_pipelines
    adjust_workflows
    change_protocols
    change_parameter_files
    make_samplesheets
    change_prepare_scripts
    cohort_specific_steps
    make_pipeline_scripts
    cd -
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -c cohort -p project_directory -t tmp04 -s samplesheet.txt -r /path/to/fastq/"
    echo "  -c      provide the cohort to prepare pipeline for (mandatory)"
    echo "  -p      Base of the project_dir where pipeline scripts will be put (mandatory)"
    echo "  -t      tmpdir of the cluster (e.g. tmp03 for boxy, tmp04 for calculon). Will use /groups/<tmp>/.../"
    echo "  -s      Samplesheet to use to get sample info from"
    echo "  -r      Location of directory that contains RNAseq files, usually either fastq or bam files"
    echo "  -h      display help"
    exit 1
}

clone_pipelines(){
    # Clone the relevant github repositories and move some files around
    echo "Cloning pipeline...."

    # In case the program has already run in current direcotry
    rm -rf Public_RNA-seq_QC
    rm -rf Public_RNA-seq_quantification
    git clone https://github.com/harmjanwestra/molgenis-pipelines.git
    # Move the relevant pipeline directories to the main project_dir
    mv molgenis-pipelines/compute5/Public_RNA-seq_QC/ .
    mv molgenis-pipelines/compute5/Public_RNA-seq_quantification/ .
    rm -rf molgenis-pipelines
}

adjust_workflows(){
    # Adjust the molgenis compute workflows to only run the steps necesarry for this project and to add new steps
    echo "Adjusting workflows..."

    # Remove unused workflows to make workflow directory less messy
    mv Public_RNA-seq_QC/workflows/workflowSTAR.csv /tmp/workflowSTAR.csv
    rm Public_RNA-seq_QC/workflows/*
    mv /tmp/workflowSTAR.csv Public_RNA-seq_QC/workflows/

    # Remove all the steps we don't need from the workflow
    sed -i '/VerifyBamID/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i '/VariantEval/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i '/GatkUnifiedGenotyper/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i '/AddOr/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv

    # Change/add steps that are not in by default
    # Have the collect metrics depend on STARMapping. Leave in dependence on CreateCramFiles in for a bit, will be removed later
#    sed -i 's/SortBam/CreateCramFiles/' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i 's;Kallisto,../protocols/Kallisto.sh,;Sailfish,../protocols/Sailfish.sh,;' Public_RNA-seq_quantification/workflows/workflow.csv
}



change_protocols(){
    # Adjust the molgenis compute protocols based on needs of the brain eQTL project
    echo "Changing protocols..."

    # add a line to delete the unfiltered BAM and sorted BAM after cramming
    echo "rm \${sortedBamDir}/\${uniqueID}.bam" >> Public_RNA-seq_QC/protocols/CreateCramFiles.sh

    # give createcramfiles a bit more mem
    sed -i 's;mem=8gb;mem=10gb;' Public_RNA-seq_QC/protocols/CreateCramFiles.sh
    
    # give CollectMultipleQcMetrics a few less cores
    sed -i 's;ppn=4;ppn=2;' Public_RNA-seq_QC/protocols/CollectMultipleQcMetrics.sh
    sed -i 's;ParallelGCThreads=4;ParallelGCThreads=2;' Public_RNA-seq_QC/protocols/CollectMultipleQcMetrics.sh
    sed -i 's;mem=4gb;mem=5gb;' Public_RNA-seq_QC/protocols/CollectMultipleQcMetrics.sh

    # Because we first convert to cram before running the collectMetrics jobs, change this for all Collect*sh scripts
    # This has changed for the new cohorts like Brainseq, keeping this in temporarily but will be removed later
#    sed -i 's;#string sortedBam;#string cramFileDir;' Public_RNA-seq_QC/protocols/Collect*sh
#    sed -i 's;#string sortedBai;#string uniqueID;' Public_RNA-seq_QC/protocols/Collect*sh
#    sed -i 's;${sortedBam};${cramFileDir}${uniqueID}.cram;' Public_RNA-seq_QC/protocols/Collect*sh

    # Did not add readgroup information at this point, so remove line METRIC_ACCUMULATION_LEVEL
    # Otherwise, tries to do it per readgroud. Now does it per file
    sed -i '/METRIC_ACCUMULATION_LEVEL/d' Public_RNA-seq_QC/protocols/CollectRnaSeqQcMetrics.sh

    # Original SAM to BAM conversion is done in seperate step, but this step is removed from the workflow
    # Add the conversion to the STAR alignment script
    echo "echo \"convert SAM to BAM\"" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "module load SAMtools/1.16.1-GCCcore-11.3.0" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "mkdir -p ${project_dir}/results/\${unfilteredBamDir}/" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "samtools view -h -b \${alignmentDir}/\${uniqueID}.sam > \${unfilteredBamDir}/\${uniqueID}.bam" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "rm \${alignmentDir}/\${uniqueID}.sam" >> Public_RNA-seq_QC/protocols/STARMapping.sh

    # Removed the filteredBam step, so change this in the SortBam protocol
    sed -i 's;filteredBam;unfilteredBam;' Public_RNA-seq_QC/protocols/SortBam.sh

    # remove unfilteredBam after sorting
    echo "rm \${unfilteredBam}" >> Public_RNA-seq_QC/protocols/SortBam.sh


    # Since we already converted to cram, change the bam part in HTSeq to cram
    sed -i 's;${bam};${cramFileDir}${uniqueID}.cram;' Public_RNA-seq_quantification/protocols/HtseqCount.sh
    sed -i 's;#string bam;#string cramFileDir;' Public_RNA-seq_quantification/protocols/HtseqCount.sh

    # Same for Sailfish, already converted to Cram so need to convert to fastq first
    # Since it's lot of lines it it is put in modified protocol
    rsync -P $script_dir/modified_protocols/Sailfish.sh Public_RNA-seq_quantification/protocols/
}

change_parameter_files(){
    # Change those parameters in parameters.csv that are old (e.g. old path to resource files) or cluster specific (e.g. tmp03/tmp04)
    echo "Changing parameter files..."

    # Change the alignment pipeline parameter file
    # NOTE: The replacemets are cluster specific, but don't want to make this command line options because that would make too many of them
    #       Either change below code or change the parameters.csv file
    sed -i 's;group,umcg-wijmenga;group,umcg-biogen;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i "s;resDir,/groups/umcg-wijmenga/tmp04/resources/;resDir,/groups/umcg-biogen/tmp01/apps/data;" Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i "s;resDir,/groups/umcg-wijmenga/tmp03/resources/;resDir,/groups/umcg-biogen/tmp01/apps/data;" Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i "s;projectDir,\${root}/\${group}/\${tmp}/projects/umcg-ndeklein/\${project}/;projectDir,${project_dir}/results/;" Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;alignmentDir,${projectDir}/hisat/;alignmentDir,${projectDir}/star;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;fastqExtension,.gz;fastqExtension,.fq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;goolf-1.7.20;foss-2018b;g' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;foss-2018b;foss-2018b;g' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;fastqExtension,.fq.gz;fastqExtension,.fastq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;1.102-Java-1.7.0_80;1.119-Java-1.7.0_80;' Public_RNA-seq_QC/parameter_files/parameters.csv


    sed -i 's;genesRefFlat,${resDir}/Ensembl/release-${ensemblVersion}/gtf/homo_sapiens/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.refflat;genesRefFlat,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.collapsedGenes.refflat;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;STARindex,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/STAR/${starVersion}/;STARindex,${resDir}/UMCG/STAR_index/gencode_32/STAR_genome_GRCh38_gencode.v32.primaryAssembly_oh100/;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;onekgGenomeFasta,${resDir}/${genomeBuild}/indices/human_g1k_v${human_g1k_vers}.fasta;onekgGenomeFasta,${resDir}//ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;rRnaIntervalList,${resDir}//picard-tools/Ensembl${ensemblVersion}/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.rrna.interval_list;rRnaIntervalList,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.collapsedGenes.rRNA.interval_list;' Public_RNA-seq_QC/parameter_files/parameters.csv

    sed -i 's;genomeBuild,b37;genomeBuild,b38;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;genomeGrchBuild,GRCh37;genomeGrchBuild,GRCh38;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;human_g1k_vers,37;human_g1k_vers,38;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;ensemblVersion,75;ensemblVersion,?;' Public_RNA-seq_QC/parameter_files/parameters.csv

    HOSTNAME=$(hostname)
    if [ "$HOSTNAME" = "gearshift" ];
    then
        echo "Change some parameters specifically for gearshift because different versions of tools are installed there"
	sed -i s';Molgenis-Compute/v16.04.1-Java-11-LTS;Molgenis-Compute/v19.01.1-Java-8-LTS;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
#        sed -i s';picardVersion,2.18.26-Java-1.8.0_74;picardVersion,2.20.5-Java-11-LTS;' Public_RNA-seq_QC/parameter_files/parameters.csv
#        sed -i s';iolibVersion,1.14.6-foss-2018b;iolibVersion,1.14.11-GCCcore-7.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv
#        sed -i s';samtoolsVersion,1.5-foss-2018b;samtoolsVersion,1.5-GCCcore-7.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv

	# edits 01-feb-2023:
       sed -i s';samtoolsVersion,1.5-foss-2018b;samtoolsVersion,1.16.1-GCCcore-11.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';samtoolsVersion,1.5-foss-2015b;samtoolsVersion,1.16.1-GCCcore-11.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';starVersion,2.6.1c-foss-2015b;starVersion,2.6.1c-foss-2018b;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';pythonVersion,3.6.4-intel-2018a;pythonVersion,3.7.4-GCCcore-7.3.0-bare;' Public_RNA-seq_QC/parameter_files/parameters.csv

       sed -i s';bedtoolsVersion,2.25.0-foss-2015b;bedtoolsVersion,2.30.0-GCCcore-11.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';cutadaptVersion,1.8.1-foss-2015b-Python-2.7.9;cutadaptVersion,4.2-GCCcore-11.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';picardVersion,2.18.26-Java-1.8.0_74;picardVersion,2.20.5-Java-11-LTS;' Public_RNA-seq_QC/parameter_files/parameters.csv
#	can't seem to find iolib; I hope this will work.. :/
       sed -i s';iolibVersion,1.14.6-foss-2015b;iolibVersion,1.14.11-GCCcore-7.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv
	
       sed -i s';tabixVersion,0.2.6-foss-2015b;tabixVersion,1.16-GCCcore-11.3.0;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';jdkVersion,1.8.0_25;jdkVersion,Java/8-LTS;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';fastqcVersion,0.11.3-Java-1.7.0_80;fastqcVersion,0.11.9-Java-8-LTS;' Public_RNA-seq_QC/parameter_files/parameters.csv
       sed -i s';RVersion,3.2.1-foss-2015b;RVersion,3.6.1-foss-2018b-bare;' Public_RNA-seq_QC/parameter_files/parameters.csv

#       sed -i s';;;' Public_RNA-seq_QC/parameter_files/parameters.csv
	
    fi


}


make_samplesheets(){
    # TODO: added -s option for giving samplesheet as input, but currently most cohorts still have hardcoded paths.
    #       same for newly added option -r
    # Have to update this
    # Make the samplesheets. How the samplesheet is made is dependent om the cohort
    echo "Making samplesheets..."
    rm Public_RNA-seq_QC/samplesheet1.csv
    mkdir -p Public_RNA-seq_QC/samplesheets/

    samplesheet_script_dir=$script_dir/samplesheet_scripts/
    if [[ "$cohort" == "TargetALS" ]];
    then
        python $samplesheet_script_dir/make_samplesheet_TargetALS.py $project_dir/../sample_annotation/UMG-Target_ALS_RNA_Clinical_Data_06122018.txt \
                                                                    $project_dir/../data/ \
                                                                     Public_RNA-seq_QC/samplesheets/samplesheet_TargetALS_RNA.

        python $samplesheet_script_dir/make_samplesheet_TargetALS.py $project_dir/../sample_annotation/RNA_Metadata_TALS_2_5July2018.txt \
                                                                     $project_dir/../data/ \
                                                                     Public_RNA-seq_QC/samplesheets/samplesheet_TargetALS_RNA.samplesheet5july2018_

        # give this an extra prefix because tehre is double data in this file and the Clinical_Data file. Since jobs have already been run for the
        # Clinical_Data file, I don't want those to be overwritten by jobs based on RNA_Metadata. Have to check afterwards by hand which are double
        python $samplesheet_script_dir/make_samplesheet_TargetALS.py $project_dir/../sample_annotation/Target_ALS_RNA_Metadata_03062019.txt \
                                                                     $project_dir/../data/ \
                                                                     Public_RNA-seq_QC/samplesheets/samplesheet_TargetALS_RNA.03062019_
    elif [[ "$cohort" == "CMC" ]];
    then
        python $samplesheet_script_dir/make_samplesheet_CMC.py /groups/umcg-biogen/tmp03/input/CMC/CMC_RNAseq_samplesheet.txt \
                                                                /groups/umcg-biogen/tmp03/input/CMC/pipelines/results/fastq/
    elif [[ "$cohort" == "CMC_HBCC" ]];
    then
        # here RNAseqDir is the directory containing the BAM files, did not want to make an extra variable for it
        echo "python $samplesheet_script_dir/make_samplesheet_CMC_HBCC.py $samplesheet $RNAseqDir $project_dir/results/fastq/ $project_dir"
        python $samplesheet_script_dir/make_samplesheet_CMC_HBCC.py $samplesheet $RNAseqDir $project_dir/results/fastq/
    elif [[ "$cohort" == "Braineac" ]];
    then
        python $samplesheet_script_dir/make_samplesheet_Braineac.py /groups/umcg-biogen/tmp03/input/ucl-upload-biogen/data/fastq/ \
                                                                    data/original_fastq_files.txt \
                                                                    SampleInfoBiogen.csv
    elif [[ "$cohort" == "Brainseq" ]];
    then
        python $samplesheet_script_dir/make_samplesheet_Brainseq.py /groups/umcg-biogen/tmp03/input/rawdata/Brainseq/RNAseq_merged/ \
                                                                    /groups/umcg-biogen/tmp03/input/rawdata/Brainseq/phenotype_data/phenotypeFile_LIBD_szControl.csv
    elif [[ "$cohort" == "NABEC" ]];
    then
        python $samplesheet_script_dir/make_samplesheet_NABEC.py /groups/umcg-biogen/tmp03/input/rawdata/NABEC/phenotypes/NABEC_phenotypes.txt /groups/umcg-biogen/tmp03/input/rawdata/NABEC/fastq/
    elif [[ "$cohort" == "BipSeq" ]];
    then
        # psychEncode has multiple datasets, easiest is to have separate script for creating samplesheet
        # samplesheet can be downloaded from Synapse
        python $samplesheet_script_dir/make_samplesheet_BipSeq.py /groups/umcg-biogen/tmp04/input/rawdata/psychEncode/metadata/BipSeq/LIBD-JHU_BipSeq_R01MH105898_Metadata_RNAseq_DLPFC_Limbic.csv \
                                                                  /groups/umcg-biogen/tmp04/input/rawdata/psychEncode/RNAseq/BipSeq/DLPFC_merged/
    elif [[ "$cohort" == "BrainGVEx" ]];
    then
        # psychEncode has multiple datasets, easiest is to have separate script for creating samplesheet
        # samplesheet can be downloaded from Synapse
        python $samplesheet_script_dir/make_samplesheet_BrainGVEx.py /groups/umcg-biogen/tmp04/input/rawdata/psychEncode/metadata/BrainGVEX/UIC-UChicago-U01MH103340_BrainGVEX_RNAseqMetadata_August2016Release.csv \
                                                                     /groups/umcg-biogen/tmp04/input/rawdata/psychEncode/RNAseq/BrainGVEX/merged
    elif [[ "$cohort" == "EpiGABA" ]];
    then
        # psychEncode has multiple datasets, easiest is to have separate script for creating samplesheet
        # samplesheet can be downloaded from Synapse
        echo 'ERROR: Not implemented yet'
        exit 1
        python $samplesheet_script_dir/make_samplesheet_EpiGABA.py
    elif [[ "$cohort" == "UCLA_ASD" ]];
    then
        # psychEncode has multiple datasets, easiest is to have separate script for creating samplesheet
        # samplesheet can be downloaded from Synapse
        python $samplesheet_script_dir/make_samplesheet_UCLA_ASD.py $samplesheet \
                                                                    $RNAseqDir
    elif [[ "$cohort" == "MSBB" ]];
    then
        python $samplesheet_script_dir/make_samplesheet_MSBB.py $samplesheet \
                                                                $RNAseqDir
    elif [[ "$cohort" == "WUSTL" ]];
    then
        python $samplesheet_script_dir/make_samplesheet_WUSTL.py $samplesheet \
                                                                $RNAseqDir
    elif [[ "$cohort" == "ENA" ]];
    then
        echo "ERROR: need to get genotypes as well for the ENA samples. Because the pipeline needs to be set up quite differently, use make_alignmentQuantificationAndGenotype_pipeline.sh instead"
        exit 1;
    else
        echo "No code written for cohort $cohort"
        exit 1;
    fi
}

change_prepare_scripts(){
    # Change the molgenis compute prepare scirpts
    echo "Chaning prepare scripts..."

    # Do general changes for the alignment pipeline
    sed -i 's;/path/to/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;parameters.converted.csv;${project_dir}/Public_RNA-seq_QC/parameter_files/parameters.converted.csv;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;workflow.csv;${project_dir}/Public_RNA-seq_QC/workflows/workflowSTAR.csv;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;-rundir /groups/umcg-wijmenga/tmp04/umcg-ndeklein/rundirs/QC/;-rundir ${project_dir}/alignment//alignmentDir;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/samplesheets/;${project_dir}/Public_RNA-seq_QC/samplesheets/;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh

    # Do general changes for the quantification pipeline
    sed -i "s;workflows/workflow.csv;${project_dir}/Public_RNA-seq_quantification/workflows/workflow.csv;" Public_RNA-seq_quantification/prepare_quantification.sh
    sed -i "s;parameter_files/;${project_dir}/Public_RNA-seq_quantification/parameter_files/;" Public_RNA-seq_quantification/prepare_quantification.sh


    # Do specific changes per samplesheet batch
    mkdir Public_RNA-seq_QC/prepare_scripts/
    mkdir Public_RNA-seq_quantification/prepare_scripts/
    for samplesheet in `ls Public_RNA-seq_QC/samplesheets/samplesheet*txt`;
    do
        name=$(echo ${samplesheet} | awk -F"." '{print $2}')
        # Change it for the alignment pipeline per batch
        echo "making prepare script using $samplesheet: Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh"
        rsync -P Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
        sed -i "s;samplesheet.csv;samplesheet_${cohort}_RNA.$name.txt;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
        sed -i "s;alignmentDir;$name;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh

        # Change it for the quantificaiton pipeline per batch
        rsync -P Public_RNA-seq_quantification/prepare_quantification.sh Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;samplesheet.csv;${project_dir}/Public_RNA-seq_QC/samplesheets/samplesheet_${cohort}_RNA.$name.txt;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;alignmentDir;$name;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;-rundir results/;-rundir ${project_dir}/quantification/$name;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
    done
}


make_pipeline_scripts(){
    # Make the pipeline scripts using Molgenis Compute
    echo "Making pipeline scripts..."
    echo $PWD

    # convert the parameter files to wide format now instead of in the change_parameters function because the cohort_specific_steps function might have changed some of them
    bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_QC/parameter_files/parameters.csv Public_RNA-seq_QC/parameter_files/parameters.converted.csv
    bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_quantification/parameter_files/parameters.csv Public_RNA-seq_quantification/parameter_files/parameters.converted.csv
    cd Public_RNA-seq_QC/prepare_scripts;
    for f in *sh;
    do
        bash $f;
    done
    cd -

#    cd Public_RNA-seq_quantification/prepare_scripts
#    for f in *sh;
#    do
#        bash $f;
#    done
#    cd -
}


cohort_specific_steps(){
    # Collection of cohort specific steps
    echo "Running cohort specific steps..."

    if [[ "$cohort" == "TargetALS" ]];
    then
#    rsync -P STARMappingTwoPass.sh Public_RNA-seq_QC/protocols/
#   echo "STARMappingTwoPass,../protocols/STARMappingTwoPass.sh,CreateCramFiles" >> Public_RNA-seq_QC/workflows/workflowSTAR.csv
        echo "TargetALS specific methods..."
        echo "HtseqCountTwoPass,../protocols/HtseqCountTwoPass.sh," >> Public_RNA-seq_quantification/workflows/workflow.csv
        rsync -P $script_dir/modified_protocols/HtseqCountTwoPass.sh Public_RNA-seq_quantification/protocols/
        rsync -P $script_dir/modified_protocols/HtseqCountOneSampleTwoPass.sh Public_RNA-seq_quantification/protocols/HtseqCount.sh

        # the project structure is different for TargetALS, so change (after changing directory, see project structure)
        for prepare_script in Public_RNA-seq_QC/prepare_scripts/* Public_RNA-seq_quantification/prepare_scripts/*;
        do
            sed -i 's;-rundir /groups/umcg-biogen/tmp04/biogen/input/TargetALS/;-rundir /groups/umcg-biogen/tmp04/biogen/input/TargetALS/pipelines/DEXSEQ_test;' $prepare_script
        done
        # same also for resultsdir
        sed -i 's;projectDir,${root}/${group}/${tmp}/biogen/input/TargetALS/results/;projectDir,${root}/${group}/${tmp}/biogen/input/TargetALS/pipelines/results/DEXSEQ_test;' Public_RNA-seq_quantification/parameter_files/parameters.csv
        sed -i 's;projectDir,${root}/${group}/${tmp}/biogen/input/TargetALS/results/;projectDir,${root}/${group}/${tmp}/biogen/input/TargetALS/pipelines/results/DEXSEQ_test;' Public_RNA-seq_QC/parameter_files/parameters.csv

        echo "unfilteredTwoPassBamDir,/groups/umcg-biogen/tmp04/biogen/input/TargetALS/pipelines/results/DEXSEQ_test/unfilteredTwoPassBam/" >> Public_RNA-seq_quantification/parameter_files/parameters.csv
    fi
    if [[ "$cohort" == "CMC_HBCC" ]];
    then
        rsync -P $script_dir/modified_protocols/ConvertBamToFastq.sh Public_RNA-seq_QC/protocols/
        sed -i '2iConvertBamToFastq,../protocols/ConvertBamToFastq.sh,' Public_RNA-seq_QC/workflows/workflowSTAR.csv
        sed -i 's;STARMapping,../protocols/STARMapping.sh,;STARMapping,../protocols/STARMapping.sh,ConvertBamToFastq;' Public_RNA-seq_QC/workflows/workflowSTAR.csv
        sed -i 's;FastQC,../protocols/Fastqc.sh,;FastQC,../protocols/Fastqc.sh,ConvertBamToFastq;' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    fi
}

parse_commandline(){
    # Parse the command line arguments
    echo "Parsing commandline..."

    # Check to see if at least one argument is given
    if [ $# -eq 0 ]
    then
        echo "ERROR: No arguments supplied"
        usage
        exit 1;
    fi

    while [[ $# -ge 1 ]]; do
        case $1 in
            -c | --cohort )         shift
                                    cohort=$1
                                    ;;
            -p | --project_dir )    shift
                                    project_dir=$1
                                    ;;
            -t | --tmpdir )         shift
                                    tmpdir=$1
                                    ;;
            -s | --samplesheet )    shift
                                    samplesheet=$1
                                    ;;
            -r | --RNAseqDir )      shift
                                    RNAseqDir=$1
                                    ;;
            -h | --help )           usage
                                    exit
                                    ;;
            * )                     echo "ERROR: Undexpected argument: $1"
                                    usage
                                    exit 1
        esac
        shift
    done
    echo "Cohort = $cohort"
    echo "Project directory = $project_dir"

    # if -z tests if variable is empty. Make sure the relevant variables are set
    if [ -z "$RNAseqDir" ];
    then
        echo "ERROR: -r/--RNAseqDir not set!"
        usage
        exit 1;
    fi
    if [ -z "$cohort" ];
    then
        echo "ERROR: -c/--cohort not set!"
        usage
        exit 1;
    fi
    if [ -z "$project_dir" ];
    then
        echo "ERROR: -p/--project_dir not set!"
        usagae
        exit 1;
    fi
    if [ -z "$tmpdir" ];
    then
        echo "ERROR: -t/--tmpdir not set!"
        usagae
        exit 1;
    fi
    if [ -z "$samplesheet" ];
    then
        echo "ERROR: -s/--samplesheet not set!"
        usagae
        exit 1;
    fi
    echo "Making pipeline..."
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
main "$@"
