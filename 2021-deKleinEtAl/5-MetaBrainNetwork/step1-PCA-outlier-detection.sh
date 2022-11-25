#!/bin/bash
# This runs the complete GeneNetwork processing pipeline

# TODO: because after PCA and smaple removal this has to be rerun, the removing duplicate
#       samples is in double for the second round. Can organize this better

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
# directory containing RunV13.jar
jar_dir=
# github dir
github_dir=
# GTF file
gtf=
# File with samples to use
sample_file=
# memory to use when doing normalization, PCA etc
mem=
# Current step
step=
####

main(){
    parse_commandline "$@"
    print_command_arguments
    mkdir -p $output_dir
    echo "0_remove_genes"
    0_remove_genes
    echo "1_select_samples"
    1_select_samples
    echo "2_remove_duplicate_samples"
    2_remove_duplicate_samples
    echo "3_quantileNormalized"
    3_quantileNormalized
    echo "Finished quantile normalizing, check the output PCA to see if any outliers need to be removed, otherwise go to next step"
}


print_command_arguments(){
    echo "Processing expression data for GeneNetwork."
    echo "Starting program with:"
    echo "expression_file=$expression_file"
    echo "output_dir=$output_dir"
    echo "project_dir=$project_dir"
    echo "jar_dir=$jar_dir"
    echo "github_dir=$github_dir"
    echo "gtf=$gtf"
    echo "TMPDIR=$TMPDIR"
    echo "sample_file=$sample_file"
    echo "mem=$mem"
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
        mv $TMPDIR/0_remove_genes ${output_dir}/
    else
        echo "{$outfile_step0} exists, go to step 1"
    fi
}

1_select_samples(){
    # Step 1. Select samples from expression file
    new_output_file_step1=$output_dir/1_selectSamples_$step/$(basename ${outfile_step0%.txt.gz})_extractedColumnsnoVarianceRowsRemoved.txt.gz
    echo $new_output_file_step1
    if [ ! -f ${new_output_file_step1} ]
    then
        echo "start step 1"
        bash $github_dir/scripts_per_step/1_select_samples.sh \
            -e $outfile_step0 \
            -p $project_dir \
            -o $TMPDIR/1_selectSamples_$step/ \
            -c $github_dir/config_file_templates/ \
            -g $github_dir/ \
            -s $sample_file
        # input of next file expect postfix of extractedColumnsnoVarianceRowsRemoved.txt.gz but is extractedColumns_noVarRemoved.txt.gz
        # change this
        f1=$output_dir/1_selectSamples_$step/$(basename ${outfile_step0%.txt.gz})_extractedColumns_noVarRemoved.txt.gz
        echo  "mv $TMPDIR/1_selectSamples_$step/ $output_dir/"
        mv $TMPDIR/1_selectSamples_$step/ $output_dir/
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
    output_file_step2="${output_dir}/2_removeDuplicates_$step/$(basename ${expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumns.txt.gz"
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

3_quantileNormalized(){
    # Step 3. calculate PCs of quantile normalized data. Do visual inspection, if samples need to be removed, rerun with different sample file as input
    # and remake the PC plot. Repeat until no outliers are left over. Then go to next step of GeneNetwork
    output_file_step3=$output_dir/3_quantileNormalized_$step/$(basename ${output_file_step2%.txt.gz}.QuantileNormalized.txt.gz)
    if [ ! -f $output_dir/3_quantileNormalized_$step/pca.png ];
    then
       bash $github_dir/scripts_per_step/3_PCA_on_quantNormalizedData.sh \
           -p $project_dir \
           -e $output_file_step2 \
           -o $TMPDIR/3_quantileNormalized_$step/ \
           -c $github_dir/config_file_templates/ \
           -j $jar_dir \
           -m $mem

       # The export.sh file has hardcoded paths for running PCA, change these
#       rsync -vP $github_dir/GeneNetwork/scripts_per_step/run_PCA.sh $TMPDIR/3_quantileNormalized/run_PCA.sh
#       sed -i "s;REPLACEGENECOVARIANCE;$TMPDIR/3_quantileNormalized/gene_covariance.txt;" $TMPDIR/3_quantileNormalized/run_PCA.sh
#       sed -i "s;REPLACEOUT;$TMPDIR/3_quantileNormalized//;" $TMPDIR/3_quantileNormalized/run_PCA.sh
#       sed -i "s;REPLACEPRECOR;$TMPDIR/3_quantileNormalized/pre_Correlation_Or_Covariance.txt;" $TMPDIR/3_quantileNormalized/run_PCA.sh
#       sbatch --wait $TMPDIR/3_quantileNormalized/run_PCA.sh
#       wait
        mv $TMPDIR/3_quantileNormalized_$step/ $output_dir/
    fi
}


usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname \\"
    echo "                    -t TMPDIR -e expression_file\\"
    echo "                    -s sample_file -q step \\"
    echo "                    -o output_dir -p project_dir \\"
    echo "                    -j jar_dir -g github_dir -m mem \\"
    echo "                    -d GeneNetworkDir -m memory -a gtf"
    echo "  -t      TMPDIR where files will be written during runtime"
    echo "  -e      Expression file"
    echo "  -s      File with list of samples to include (exclude all others)"
    echo "  -q      Current step, e.g. 1a or 1b. This way can distinguish between multiple PCA rounds"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output directory where results will be written"
    echo "  -j      Location of RunV13.jar and RunV13.jar"
    echo "  -g      Github GeneNetwork directory"
    echo "  -a      GTF file"
    echo "  -m      Memory to give to the eqtlgen jar file when running, will be appendended to -Xmx, so 8g will be java -jar -Xmx8g"
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
            -c | --config_template_dir )    shift
                                            config_template_dir=$1
                                            ;;
            -j | --jar_dir )                shift
                                            jar_dir=$1
                                            ;;
            -g | --github_dir )             shift
                                            github_dir=$1
                                            ;;
            -a | --gtf )                    shift
                                            gtf=$1
                                            ;;
            -s | --sample_file )            shift
                                            sample_file=$1
                                            ;;
            -q | --step )                   shift
                                            step=$1
                                            ;;
            -m | --mem )                    shift
                                            mem=$1
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
    if [ -z "$github_dir/config_file_templates/" ];
    then
        echo "ERROR: -c/--config_template_dir not set!"
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
    if [ -z "$sample_file" ];
    then
        echo "ERROR: -s/--sample_file not set!"
        usage
    exit 1;
    fi
    if [ -z "$step" ];
    then
        echo "ERROR: -q/--step not set!"
        usage
    exit 1;
    fi
    if [ -z "$mem" ];
    then
        echo "ERROR: -m/--mem not set!"
        usage
    exit 1;
    fi
}

main "$@"





