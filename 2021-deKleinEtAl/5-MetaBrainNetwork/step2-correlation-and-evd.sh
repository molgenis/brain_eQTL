#!/bin/bash
# This runs the complete GeneNetwork processing pipeline

set -e
set -u

module load Python

#### Only lines that should be changed ####
# TMPDIR for writing files
TMPDIR=
# original expression table
expression_file=
# Directory to write output to
output_dir=
# Directory to write configuration files to (will be put in $project_dir/configs/)
project_dir=
# directory containing V13.jar
jar_dir=
# file with samples to analyze
sample_file=
# github dir
github_dir=
# covariate table
covar_table=
# GTF file
gtf=
# number of threads to use for correlation and PCA step
threads=
# memory to use when doing normalization, PCA etc
mem=
# qos to run sbatch jobs in (default = regular)
qos=
# Name of the project (e.g. GeneNetwork, MetaBrain, KidneyNetwork, etc)
name=
# the previous step which has the correctly selected samples. e.g. step1b, will use this to find right input directory
prev_step=
####

main(){
    parse_commandline "$@"
    print_command_arguments
    cd $github_dir/GeneNetwork
    if [ ! -f RunV13.jar ];
    then
        echo "RunV13.jar does not exist, download now"
        wget https://github.com/npklein/brain_eQTL/releases/download/2020-03-02/RunV13.jar
    fi
    md5=$(md5sum RunV13.jar | awk '{print $1}')
    if [ $md5 != "9b292956792206b7cdb8a4e7fee5a06b" ];
    then
        echo "ERROR: md5sum of RunV13.jar not correct"
        exit 1;
    fi
    cd -
    mkdir -p $output_dir
    echo "0_remove_genes"
    0_remove_genes
    echo "1_select_samples"
    1_select_samples
    echo "2_remove_duplicate_samples"
    2_remove_duplicate_samples
    echo "4_DeseqNormalizedData"
    4_DeseqNormalizedData
    echo "5_RemoveCovariates"
    5_RemoveCovariates
    echo "6_CorrelationMatrix"
    6_CorrelationMatrix
    echo "7_evd_on_correlation"
    7_evd_on_correlation
    echo "Finished correlation and eigenvalue decomposition, check output/7_evd_on_correlation_matrix/eigenvalues.png how many eigenvalues to use."
    echo "Then run the next step."
}


print_command_arguments(){
    echo "-------------------------------------------"
    echo "Processing expression data for GeneNetwork."
    echo "-------------------------------------------"
    echo "Starting program with:"
    echo "expression_file=$expression_file"
    echo "output_dir=$output_dir"
    echo "project_dir=$project_dir"
    echo "jar_dir=$jar_dir"
    echo "sample_file=$sample_file"
    echo "github_dir=$github_dir"
    echo "covar_table=$covar_table"
    echo "gtf=$gtf"
    echo "mem=$mem"
    echo "qos=$qos"
    echo "name=$name"
    echo "-------------------------------------------"
    echo ""
    echo ""
}


0_remove_genes(){
    # Step 1. Select samples from expression file
    outfile_step0="${output_dir}/0_remove_genes/$(basename ${expression_file%.txt.gz}).genesRemoved.txt.gz"
    if [ ! -f ${outfile_step0} ]
    then
        echo "${outfile_step0} does not exist, start step 0"
        python $github_dir/scripts_per_step/0_remove_duplicate_sameCounts_scaffolded_genes.py \
            -e $expression_file \
            -o $TMPDIR/0_remove_genes/$(basename $outfile_step0) \
            -g $gtf
        mkdir -p ${output_dir}/0_remove_genes/
        mv $TMPDIR/0_remove_genes/* ${output_dir}/0_remove_genes/
        rmdir $TMPDIR/0_remove_genes/
    else
        echo "$outfile_step0 exists, go to step 1"
    fi
}

1_select_samples(){
    # Step 1. Select samples from expression file
    new_output_file_step1=$output_dir/1_selectSamples_$prev_step/$(basename ${outfile_step0%.txt.gz})_extractedColumnsnoVarianceRowsRemoved.txt.gz
    echo $new_output_file_step1
    if [ ! -f ${new_output_file_step1} ]
    then
        echo "start step 1"
        bash $github_dir/scripts_per_step/1_select_samples.sh \
            -e $outfile_step0 \
            -p $project_dir \
            -o $TMPDIR/1_selectSamples_$prev_step/ \
            -c $github_dir/config_file_templates/ \
            -g $github_dir/ \
            -s $sample_file
        # input of next file expect postfix of extractedColumnsnoVarianceRowsRemoved.txt.gz but is extractedColumns_noVarRemoved.txt.gz
        # change this
        f1=$output_dir/1_selectSamples_$prev_step/$(basename ${outfile_step0%.txt.gz})_extractedColumns_noVarRemoved.txt.gz
        mkdir -p $output_dir/1_selectSamples_$prev_step
        echo "mv $TMPDIR/1_selectSamples/* $output_dir/"
        mv $TMPDIR/1_selectSamples_$prev_step/* $output_dir/1_selectSamples_$prev_step/
        rmdir $TMPDIR/1_selectSamples_$prev_step/
        mv $f1 $new_output_file_step1
        if [ ! -f ${new_output_file_step1} ];
        then
            echo "${new_output_file_step1} does not exist"
            exit 1;
        fi
    fi
}

2_remove_duplicate_samples(){
    # We have some ENA samples in our data, so run for duplicates. Add the relevant input files to the config file
    # Step 2. remove duplicates
    output_file_step2="${output_dir}/2_removeDuplicates_$prev_step/$(basename ${expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumns.txt.gz"
    if [ ! -f ${output_file_step2} ];
    then
        echo "start step 2"
        bash $github_dir/scripts_per_step/2_remove_duplicate_samples.sh \
            -p $project_dir \
            -e $new_output_file_step1 \
            -o $TMPDIR/$(basename $output_file_step2) \
            -c $github_dir/config_file_templates/ \
            -g $github_dir/
        echo "done"
        mkdir -p $(dirname $output_file_step2)
        mv $TMPDIR/$(basename $output_file_step2) $output_file_step2
    fi
}

4_DeseqNormalizedData(){
    output_file_step4="$output_dir/4_deseqNormalized/$(basename ${output_file_step2%.txt.gz}noVarianceRowsRemoved.DESeqNorm.log2.txt)"
    if [ ! -f $output_file_step4 ];
    then
        # Step 4. Do deseq normalisation on the original counts
        bash $github_dir/scripts_per_step/4_DeseqNormalizedData.sh \
            -p $project_dir \
            -e $output_file_step2 \
            -f $TMPDIR/4_deseqNormalized/$(basename ${output_file_step2%.txt.gz}noVarianceRowsRemoved.DESeqNorm.txt.gz) \
            -o $TMPDIR/4_deseqNormalized/ \
            -g $github_dir/ \
            -z $TMPDIR/4_deseqNormalized/PCA_corrected_expression/ \
            -m $mem
        mkdir -p $output_dir/4_deseqNormalized/
        mv $TMPDIR/4_deseqNormalized/* $output_dir/4_deseqNormalized/
        rmdir $TMPDIR/4_deseqNormalized/
    fi
}

5_RemoveCovariates(){
    output_file_step5=$output_dir/5_covariatesRemoved/$(basename ${output_file_step4%.txt}.CovariatesRemoved.txt.gz)
    echo "$output_file_step5"
    if [ ! -f $output_file_step5 ];
    then
        # Step 5. Remove covariates from deseq normalized data
        bash $github_dir/scripts_per_step/5_RemoveCovariates.sh \
            -p $project_dir \
            -e $output_file_step4 \
            -o $TMPDIR/5_covariatesRemoved \
            -c $github_dir/config_file_templates/ \
            -g $github_dir/ \
            -z $covar_table \
            -m $mem
        mkdir -p $output_dir/5_covariatesRemoved/
        echo "Removing rows with NaN values"
        zcat $TMPDIR/5_covariatesRemoved/$(basename $output_file_step5) | grep -v 'NaN' > $TMPDIR/5_covariatesRemoved/$(basename ${output_file_step5}).tmp
        mv $TMPDIR/5_covariatesRemoved/$(basename ${output_file_step5}).tmp $TMPDIR/5_covariatesRemoved/$(basename ${output_file_step5%.gz})
        rm $TMPDIR/5_covariatesRemoved/$(basename $output_file_step5)
        gzip $TMPDIR/5_covariatesRemoved/$(basename ${output_file_step5%.gz})
        mv $TMPDIR/5_covariatesRemoved/* $output_dir/5_covariatesRemoved/
        rmdir $TMPDIR/5_covariatesRemoved
    else
        echo "$output_file_step5 exists, skip step"
    fi
}

6_CorrelationMatrix(){
    output_file_step6="$output_dir/6_correlation_matrix/$name.deseqNorm.covarCorrected.correlation.txt"
    if [ ! -f $output_file_step6 ]  && [ ! -f ${output_file_step6}.gz ];
    then
        # Step 6. Make correlation matrix
        bash $github_dir/scripts_per_step/6_CorrelationMatrix.sh \
            -p $project_dir \
            -e $output_file_step5 \
            -o ${output_file_step6} \
            -c $github_dir/config_file_templates/ \
            -g $github_dir/ \
            -t $threads \
            -m $mem \
            -q $qos
    fi
}

7_evd_on_correlation(){
    if [ ! -f $output_file_step6 ] && [ -f ${output_file_step6}.gz ];
    then
        echo "zcat ${output_file_step6}.gz > ${output_file_step6}"
        zcat ${output_file_step6}.gz > ${output_file_step6}
    fi

    if [ ! -f ${output_file_step5%.gz} ] && [ -f ${output_file_step5} ];
    then
        echo "zcat ${output_file_step5} > ${output_file_step5%.gz}"
         zcat ${output_file_step5} > ${output_file_step5%.gz}
    fi

    output_file_step_7=$output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.txt.gz
    if [ ! -f $output_file_step_7 ];
    then
        # step 7. Run evd on correlation matrix
        mkdir -p  $output_dir/7_evd_on_correlation_matrix/
        rsync -vP $github_dir/scripts_per_step/7_evd_on_correlation.sh $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACECORMATRIX;$output_file_step6;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACECOVCORRECTEDEXPRESSION;${output_file_step5%.gz};" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACEOUTDIR;$output_dir/;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACETHREADS;$threads;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACESCRIPTDIR;$github_dir/scripts_per_step/;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACESVDSOLVER;auto;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACEMEM;$mem;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACEQOS;$qos;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        sed -i "s;REPLACENAME;$name;" $output_dir/7_evd_on_correlation_matrix/7_evd_on_correlation.sh
        cd $output_dir/7_evd_on_correlation_matrix
        echo "submit eigenvector decomposition to cluster. Can take a while. Will wait untill finished"
        jobid=$(sbatch --parsable 7_evd_on_correlation.sh)
        echo "job submitted (jobid=$jobid), check if job is running, if so sleep for 10 minutes"
        startDate=`date`

        while [ ! $(squeue -j $jobid 2> /dev/null | wc -l) -le 1 ];
        do
            echo "--------------------------------"
            echo "sbatch starting time: $startDate"
            echo "current time: `date`"
            echo "Job still running. Sleep 10 minutes"
            sleep 600;
        done
        if [ ! -f ${output_file_step_7} ];
        then
            echo "ERROR: Outputfile ${output_file_step_7} does not exist, check if job exited with an error"
            exit 1;
        fi
        echo "done"
        cd -
    fi
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -t TMPDIR -e expression_file -o output_dir -p project_dir"
    echo "                    -j jar_dir -s sample_file -g github_dir -z covar_table -v threads"
    echo "                    -r conbach_alpha -n name -m mem -y prev_step [-q qos]"
    echo "  -t      TMPDIR where files will be written during runtime"
    echo "  -e      Expression file"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output directory where results will be written"
    echo "  -j      Location of eqtl-mapping-pipeline.jar"
    echo "  -s      File with samples to include"
    echo "  -g      Location of git cloned https://github.com/npklein/brain_eQTL/ directory"
    echo "  -z      Covariate table"
    echo "  -a      GTF file"
    echo "  -v      Number of threads to use for correlation step and PCA step"
    echo "  -m      Memory to use for some steps"
    echo "  -q      qos to run sbatch jobs in"
    echo "  -n      Name that will be used in output file"
    echo "  -y      the previous step which has the correctly selected samples. e.g. step1b, will use this to find right input directory"
    echo "  -h      display help"
    exit 1
}

parse_commandline(){
    # Check to see if at least one argument is given
    if [ $# -eq 0 ]
    then
        echo "ERROR: No arguments supplied"
        usage
        exit 1;
    fi

    while [[ $# -ge 1 ]]; do
        case $1 in
            -t | --TMPDIR )                 shift
                                            TMPDIR=$1
                                            ;;
            -p | --project_dir )            shift
                                            project_dir=$1
                                            ;;
            -e | --expression_file )        shift
                                            expression_file=$1
                                            ;;
            -o | --output_dir )             shift
                                            output_dir=$1
                                            ;;
            -j | --jar_dir )                shift
                                            jar_dir=$1
                                            ;;
            -s | --sample_file )            shift
                                            sample_file=$1
                                            ;;
            -g | --github_dir )             shift
                                            github_dir=$1
                                            ;;
            -a | --gtf )                    shift
                                            gtf=$1
                                            ;;
            -v | --threads )                shift
                                            threads=$1
                                            ;;
            -z | --covar_table )            shift
                                            covar_table=$1
                                            ;;
            -m | --mem )                    shift
                                            mem=$1
                                            ;;
            -q | --qos )                    shift
                                            qos=$1
                                            ;;
            -n | --name )                   shift
                                            name=$1
                                            ;;
            -y | --prev_step )              shift
                                            prev_step=$1
                                            ;;
            -h | --help )                   usage
                                            exit
                                            ;;
            * )                             echo "ERROR: Undexpected argument: $1"
                                            usage
                                            exit 1
        esac
        shift
    done

    # if -z tests if variable is empty. Make sure the relevant variables are set
    if [ -z "$project_dir" ];
    then
        echo "ERROR: -p/--project_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$expression_file" ];
    then
        echo "ERROR: -e/--expression_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$output_dir" ];
    then
        echo "ERROR: -o/--output_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$jar_dir" ];
    then
        echo "ERROR: -j/--jar_dir not set!"
        usage
        exit 1;
    fi
        if [ -z "$sample_file" ];
    then
        echo "ERROR: -s/--sample_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$github_dir" ];
    then
        echo "ERROR: -g/--github_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$gtf" ];
    then
        echo "ERROR: -a/--gtf not set!"
        usage
        exit 1;
    fi
    if [ -z "$covar_table" ];
    then
        echo "ERROR: -z/--covar_table not set!"
        usage
        exit 1;
    fi
    if [ -z "$threads" ];
    then
        echo "ERROR: -v/--threads not set!"
        usage
        exit 1;
    fi
    if [ -z "$mem" ];
    then
        echo "ERROR: -m/--mem not set!"
        usage
        exit 1;
    fi
    if [ -z "$qos" ];
    then
        echo "ERROR: -q/--qos not set!"
        usage
        exit 1;
    fi
    if [ -z "$name" ];
    then
        echo "ERROR: -n/--name not set!"
        usage
        exit 1;
    fi
    if [ -z "$prev_step" ];
    then
        echo "ERROR: -y/--prev_step not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
main "$@";






