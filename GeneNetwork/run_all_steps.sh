
# This runs the complete GeneNetwork processing pipeline

set -e
set -u

#### Only lines that should be changed ####
# TMPDIR for writing files
TMPDIR=
# original expression table
input_expression_file=
# Directory to write output to
output_dir=
# Directory to write configuration files to (will be put in $project_dir/configs/)
project_dir=
# directory with config file templates
config_template_dir=
# directory containing V13.jar
jar_dir=
# file with samples to analyze
sample_file=
# github dir
github_dir=
# quant normalize only necesarry to select samples. we already know which are the good ones, so skip
quant_norm=false
# covariate table
covar_table=
####

main(){
    parse_commandline "$@"
    mkdir -p $output_dir
    1_remove_duplicate_samples
    2_select_samples
    3_quantileNormalized
    4_DeseqNormalizedData
    5_RemoveCovariates
    6_CorrelationMatrix
    7_PCA_on_correlation
}
1_remove_duplicate_samples(){
    # We have some ENA samples in our data, so run for duplicates. Add the relevant input files to the config file
    # Step 1. remove duplicates
    outfile_step1="${output_dir}/1_removeDuplicates/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved.txt.gz"
    if [ ! -f ${outfile_step1} ];
    then
        echo "start step 1"
        bash $github_dir/GeneNetwork/scripts_per_step/1_remove_duplicate_samples.sh \
            -p $project_dir \
            -e $input_expression_file \
            -o $TMPDIR/$(basename $outfile_step1) \
            -c $config_template_dir \
            -j $jar_dir
        echo "done"
        echo ""
        echo ""
        echo "Do: mv $TMPDIR/$(basename $outfile_step1) $outfile_step1"
        mkdir -p $(dirname $outfile_step1)
        mv $TMPDIR/$(basename $outfile_step1) $outfile_step1
    fi
}

2_select_samples(){
    # Step 2. Select samples from expression file
    outfile_step2="${output_dir}/2_selectSamples/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumns.txt.gz"
    new_outfile_step2=$output_dir/2_selectSamples/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumnsnoVarianceRowsRemoved.txt.gz
    if [ ! -f ${new_outfile_step2} ]
    then
        echo "start step 2"
        bash $github_dir/GeneNetwork/scripts_per_step/2_select_samples.sh \
            -p $project_dir \
            -e $outfile_step1 \
            -o $TMPDIR/2_selectSamples/ \
            -c $config_template_dir \
            -j $jar_dir \
            -s $sample_file
        # input of next file expect postfix of extractedColumnsnoVarianceRowsRemoved.txt.gz but is extractedColumns_noVarRemoved.txt.gz
        # change this
        f1=$TMPDIR/2_selectSamples/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumns_noVarRemoved.txt.gz
        mv $f1 $TMPDIR/2_selectSamples/$(basename $new_outfile_step2)
        mv $TMPDIR/2_selectSamples/ $output_dir/
    fi
}

3_quantileNormalized(){
    if [ "quant_norm" = true ];
    then
        # Step 3. make PCA of quantile normalized data. Do visual inspection, if samples need to be removed, run from step 2 again (selecting only correct samples)
        # and then skip step 1 and 3 and go to step 4. I think this step is nececarry for step 4, not sure why though.
        outfile_step3=$output_dir/3_quantileNormalized/$(basename ${new_outfile_step2%.txt.gz}.QuantileNormalized.txt.gz)
        if [ ! -f $outfile_step3 ];
        then
            bash $github_dir/GeneNetwork/scripts_per_step/3_PCA_on_quantNormalizedData.sh \
                -p $project_dir \
                -e $new_outfile_step2 \
                -o $TMPDIR/3_quantileNormalized/ \
                -c $config_template_dir \
                -j $jar_dir

            # The export.sh file has hardcoded paths for running PCA, change these
            rsync -vP $github_dir/GeneNetwork/scripts_per_step/run_PCA.sh $TMPDIR/3_quantileNormalized/3_run_PCA.sh
            sed -i "s;REPLACEGENECOVARIANCE;$TMPDIR/3_quantileNormalized/gene_covariance.txt;" $TMPDIR/3_quantileNormalized/3_run_PCA.sh
            sed -i "s;REPLACEOUT;$TMPDIR/3_quantileNormalized//;" $TMPDIR/3_quantileNormalized/3_run_PCA.sh
            sed -i "s;REPLACEPRECOR;$TMPDIR/3_quantileNormalized/pre_Correlation_Or_Covariance.txt;" $TMPDIR/3_quantileNormalized/3_run_PCA.sh
            bash $TMPDIR/3_quantileNormalized/3_run_PCA.sh
            mv $MPDIR/3_quantileNormalized/ $output_dir/
        fi
    fi
}

4_DeseqNormalizedData(){
    outfile_step4=$output_dir/4_deseqNormalized/$(basename ${new_outfile_step2%.txt.gz}noVarianceRowsRemoved.DESeqNorm.txt.gz)
    if [ ! -f $outfile_step4 ];
    then
        # Step 4. Do deseq normalisation on the original counts
        bash $github_dir/GeneNetwork/scripts_per_step/4_DeseqNormalizedData.sh \
            -p $project_dir \
            -e $new_outfile_step2 \
            -o $TMPDIR/4_deseqNormalized/ \
            -c $config_template_dir \
            -j $jar_dir \
            -z $TMPDIR/4_deseqNormalized/PCA_corrected_expression/

        mv $TMPDIR/4_deseqNormalized/ $output_dir/
    fi
}

5_RemoveCovariates(){
    outfile_step5=$output_dir/5_covariatesRemoved/$(basename ${outfile_step4%.txt.gz}.CovariatesRemoved.txt.gz)
    if [ ! -f $outfile_step5 ];
    then
        # Step 5. Remove covariates from deseq normalized data
        bash $github_dir/GeneNetwork/scripts_per_step/5_RemoveCovariates.sh \
            -p $project_dir \
            -e $outfile_step4 \
            -o $TMPDIR/5_covariatesRemoved \
            -c $config_template_dir \
            -j $jar_dir \
            -z $covar_table
        mv $TMPDIR/5_covariatesRemoved $output_dir/
    fi
}

6_CorrelationMatrix(){
    outfile_step6="$output_dir/6_correlation_matrix/MetaBrain.deseqNorm.covarCorrected.correlation.txt.gz"
    if [ ! -f $outfile_step6 ];
    then
        # Step 6. Make correlation matrix
        bash $github_dir/GeneNetwork/scripts_per_step/6_CorrelationMatrix.sh \
            -p $project_dir \
            -e $outfile_step5 \
            -o $TMPDIR/$(basename $outfile_step6) \
            -c $config_template_dir \
            -j $jar_dir
        mkdir -p $(dirname $outfile_step6)
        mv $TMPDIR/$(basename $outfile_step6) $outfile_step6
    fi
}

7_PCA_on_correlation(){
    if [ ! -f $output_dir/7_PCA_on_correlation_matrix/MetaBrain.pc-scores.txt ];
    then
        # step 7. Run PCA on correlation matrix
        rsync -vP $github_dir/GeneNetwork/scripts_per_step/7_PCA_on_correlation.sh $TMPDIR/7_PCA_on_correlation.sh
        sed -i "s;REPLACEOUTDIR;$output_dir/;" $TMPDIR/7_PCA_on_correlation.sh
        sed -i "s;REPLACEOUTFILE;$outfile_step6;" $TMPDIR/7_PCA_on_correlation.sh
        sed -i "s;REPLACEPCASETUP;$github_dir/GeneNetwork/PCA_setup_calculon.sh;" $TMPDIR/7_PCA_on_correlation.sh
        sbatch $TMPDIR/7_PCA_on_correlation.sh
    fi
}


usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -t TMPDIR -e input_expression_file -o output_dir -p project_dir -c config_template_dir -j jar_dir -s sample_file -g github_dir -z covar_table -q quant_norm (default: false)"
    echo "  -t      TMPDIR where files will be written during runtime"
    echo "  -e      Expression file to remove duplciates from"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -o      Output file that will be written"
    echo "  -j      Location of V13 jar file"
    echo "  -s      File with samples to include"
    echo "  -g      Github GeneNetwork directory"
    echo "  -z      Covariate table"
    echo "  -q      true/false wether quntile normalization should be done"
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
            -t | --TMPDIR               shift
                                        project_dir=$1
                                        ;;
            -p | --project_dir )        shift
                                        project_dir=$1
                                        ;;
            -e | --expression_file )    shift
                                        expression_file=$1
                                        ;;
            -o | --outfile )            shift
                                        outfile=$1
                                        ;;
            -c | --config_templates )   shift
                                        config_templates=$1
                                        ;;
            -j | --jardir )             shift
                                        jardir=$1
                                        ;;
            -s | --sample_file )        shift
                                        sample_file=$1
                                        ;;
            -g | --github_dir )         shift
                                        github_dir=$1
                                        ;;
            -z | --covar_table )        shift
                                        covar_table=$1
                                        ;;
            -q | --quant_norm )         shift
                                        quant_norm=$1
                                        ;;
            -h | --help )               usage
                                        exit
                                        ;;
            * )                         echo "ERROR: Undexpected argument: $1"
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
    if [ -z "$outfile" ];
    then
        echo "ERROR: -o/--outfile not set!"
        usage
        exit 1;
    fi
    if [ -z "$jardir" ];
    then
        echo "ERROR: -j/--jardir not set!"
        usage
        exit 1;
    fi
    if [ -z "$config_templates" ];
    then
        echo "ERROR: -c/--config_templates not set!"
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
    if [ -z "$config_templates" ];
    then
        echo "ERROR: -c/--config_templates not set!"
        usage
        exit 1;
    fi
    if [ -z "$covar_table" ];
    then
        echo "ERROR: -z/--covar_table not set!"
        usage
        exit 1;
    fi    
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






