#!/bin/bash
# This runs the complete GeneNetwork processing pipeline

set -e
set -u

module load Python
module load RPlus

#### Only lines that should be changed ####
# TMPDIR for writing files
TMPDIR=
# Directory to write output to
output_dir=
# Directory to write configuration files to (will be put in $project_dir/configs/)
project_dir=
# directory containing V13.jar
jar_dir=
# github dir
github_dir=
# GTF file
gtf=
# number of threads to use for correlation and evd step
threads=
# gene_network_dir that contains GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar and PathwayMatrix/
gene_network_dir=
# memory to use when doing normalization, evd etc
mem=
# Name of the project (e.g. GeneNetwork, MetaBrain, KidneyNetwork, etc)
# This needs to be same as previous step
name=
# Number of eigenvectors to keep
n_eigenvectors=
# qos when submitting to cluster
qos=regular
# ENSG to hgnc gzipped file (first column ENSG, 3rd column HNCG
ensg_hgnc_file=
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
    echo "8_center_scale"
    8_center_scale
    echo "9_correlate_eigenvectors"
    9_correlate_eigenvectors
    echo "10_GeneNetwork_predictions"
    10_GeneNetwork_predictions
    echo "11_WebsiteMatrixCreator"
    11_WebsiteMatrixCreator
}


print_command_arguments(){
    echo "-------------------------------------------"
    echo "Processing expression data for GeneNetwork."
    echo "-------------------------------------------"
    echo "Starting program with:"
    echo "output_dir=$output_dir"
    echo "project_dir=$project_dir"
    echo "jar_dir=$jar_dir"
    echo "github_dir=$github_dir"
    echo "gtf=$gtf"
    echo "gene_network_dir=$gene_network_dir"
    echo "mem=$mem"
    echo "name=$name"
    echo "qos=$qos"
    echo "n_eigenvectors=$n_eigenvectors"
    echo "-------------------------------------------"
    echo ""
    echo ""
}

8_center_scale(){
    output_file_step8="$output_dir/8_CenterScaledColumnsRows/$name.${n_eigenvectors}_eigenvectors.columnsRowsCenterScaled.txt.gz"
    if [ ! -f ${output_file_step8} ];
    then
        echo "${output_file_step8}.gz does not exist yet, start step 8"
        if [ ! -f ${output_dir}/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt.gz ];
        then
            nforcut=$(expr $n_eigenvectors + 1)
            zcat < $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.txt.gz | cut -f1-$nforcut > $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt
            gzip $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt
        fi
        bash $github_dir/GeneNetwork/scripts_per_step/8_CenterScaleColumnsRows.sh \
            -i $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt \
            -o ${output_file_step8}
    else
        if [ ! -f ${output_dir}/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt.gz ];
        then
            nforcut=$(expr $n_eigenvectors + 1)
            zcat < $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.txt.gz | cut -f1-$nforcut > $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt
            gzip $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt
        fi
        echo "${output_file_step8} exists, go to step 9"
    fi
}

9_correlate_eigenvectors(){
    output_file_step9="$output_dir/9_eigenvector_correlation_matrix/$name.${n_eigenvectors}_eigenvectors.columnsRowsCenterScaled.correlation.txt"
    if [ ! -f ${output_file_step9}.gz ];
    then
        echo "${output_file_step9}.gz does not exist yet, start step 9"
        # Step 8. Make correlation matrix of eigenvectors
        bash $github_dir/GeneNetwork/scripts_per_step/9_CorrelateEigenvectors.sh \
            -p $project_dir \
            -e $output_file_step8 \
            -o ${output_file_step9} \
            -c $github_dir/GeneNetwork/config_file_templates/ \
            -j $jar_dir \
            -t $threads \
            -q $qos \
            -m $mem \
            -g $github_dir/GeneNetwork/
    else
        echo "${output_file_step9}.gz already exists, go to step 10"
    fi
}

10_GeneNetwork_predictions(){
    if [ ! -f $github_dir/GeneNetwork/GeneNetworkBackend/GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar ]
    then
        echo "Download GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar"
        cd $github_dir/GeneNetwork/GeneNetworkBackend/
        wget  'https://molgenis50.gcc.rug.nl/jenkins/job/systemsgenetics/473/nl.systemsgenetics$GeneNetworkBackend/artifact/nl.systemsgenetics/GeneNetworkBackend/1.0.7-SNAPSHOT/GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar'
        cd -
    fi

    python $github_dir/GeneNetwork/table_scripts/remove_gene_version.py $output_dir/7_evd_on_correlation_matrix/${name}.eigenvectors.${n_eigenvectors}_eigenvectors.txt

    mkdir -p  $output_dir/10_GeneNetwork_predictions/scripts/
    job_ids=()
    outfiles=()
    # if all output files already made we have to keep track, so a later step can be skipped
    at_least_1_job_submitted=false
    cd $output_dir/10_GeneNetwork_predictions/scripts/
    for f in $github_dir/GeneNetwork/GeneNetworkBackend/PathwayMatrix/*matrix.txt;
    do
        echo "Start predicting with $f"
        matrix_name=$(basename ${f%_matrix.txt})
        matrix_name=${matrix_name%.matrix.txt}
        mkdir -p $output_dir/10_GeneNetwork_predictions/scripts/
        outfile="$output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt"
        echo $matrix_name
        echo $outfile
        if [ ! -f $outfile ] && [ ! -f ${outfile}.gz ];
        then
            at_least_1_job_submitted=true
            rsync -vP $github_dir/GeneNetwork/scripts_per_step/10_GeneNetwork_predictions.sh $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACENAME;${matrix_name}.${n_eigenvectors}_eigenvectors.predictions;g" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACEGENENETWORKDIR;${gene_network_dir};" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACEIDENTITYMATRIX;$f;g" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACEEIGENVECTORS;$output_dir/7_evd_on_correlation_matrix/${name}.eigenvectors.${n_eigenvectors}_eigenvectors.txt;g" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACEOUT;$outfile;" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACEBACKGROUND;${f%matrix.txt}genesInPathways.txt;g" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACEMEM;${mem};g" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACEQOS;${qos};" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            sed -i "s;REPLACENEIG;${n_eigenvectors};" $output_dir/10_GeneNetwork_predictions/scripts/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh
            job_ids+=("$(sbatch --parsable ${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.sh)")
            outfiles+=("$outfile")
        else
            echo "$output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt exists, skip"
        fi
    done

    if [ "$at_least_1_job_submitted" = true ];
    then
        echo "sumbitted ${#job_ids[@]} jobs, with job ids ${job_ids}"
        echo "will loop over all jobs, and wait for each until it is finished"
        startDate=`date`
        for jobid in ${job_ids[@]};
        do
            echo "Checking if job $jobid is running, if so sleep for 10 minutes"
            while [ ! $(squeue -j $jobid 2>/dev/null | wc -l ) -le 1 ];
            do
                echo "--------------------------------"
                echo "sbatch starting time: $startDate"
                echo "current time: `date`"
                echo "Job $jobid still running. Sleep 2 minutes"
                sleep 120;
            done
            echo "Done!"
        done
        for outfile in ${outfiles[@]};
        do
            if [ ! -f $outfile ];
            then
                echo "ERROR: Outputfile $outfile does not exist, check if job exited with an error"
                exit 1;
            fi
        done
    else
        echo "All output files exist"
    fi
    cd -

    for f in $github_dir/GeneNetwork/GeneNetworkBackend/PathwayMatrix/*terms.txt;
    do
        matrix_name=$(basename ${f%_terms.txt})
        matrix_name=${matrix_name%.terms.txt}
        auc="${output_dir}/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.txt"
        if [ ! -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt ] && [ ! -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt.gz ];
        then
            if [ -f ${auc}.gz ] && [ ! -f ${auc} ];
            then
                echo "zcat ${auc}.gz to ${auc}"
                zcat ${auc}.gz > ${auc}
            fi
            Rscript ${github_dir}/GeneNetwork/table_scripts/calculate_bonferonni.R \
                                          -i $auc \
                                          -o $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt \
                                          -t $f
        fi
        outfile=$output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt
        if [ ! -f $outfile ] && [ ! -f ${outfile}.gz ];
        then
            if [ -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt.gz ] && [ ! -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt ];
            then
                echo "zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt.gz > zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt"
                zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt
            fi

            if [ -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt.gz ] && [ ! -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt ];
            then
                echo "zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt"
                zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt
            fi
            bash ${github_dir}/GeneNetwork/scripts_per_step/11_select_columns.sh -g ${github_dir}/GeneNetwork/ \
                                                                             -a $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt \
                                                                             -c $github_dir/GeneNetwork/config_file_templates/ \
                                                                             -p $project_dir \
                                                                             -q $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.txt \
                                                                             -o $outfile
        fi
    done

}


11_WebsiteMatrixCreator(){
    echo "# Copy the commented files to your local machine where you cloned molgenis-app-genenetwork in a separate data/ directory, then run this code" > $output_dir/11_ImportToWebsite/populate_database.sh
    for f in $github_dir/GeneNetwork/GeneNetworkBackend/PathwayMatrix/*.matrix.txt;
    do
        matrix_name=$(basename ${f%.matrix.txt})
        matrix_name=${matrix_name%_matrix.txt}
        echo $output_dir/11_ImportToWebsite/${matrix_name}.${n_eigenvectors}_eigenvectors_gnInputFormat.txt
        if [ ! -f $output_dir/11_ImportToWebsite/${n_eigenvectors}/${matrix_name}.${n_eigenvectors}_eigenvectors.matrix_gnInputFormat.txt ];
        then
            if [ ! -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt ] && [ -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt.gz ];
            then
                echo "zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt.gz > output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt"
                zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt
            fi

            if [ ! -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt ] && [ -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt.gz ];
            then
                echo "zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt"
                zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt
            fi

            if [ ! -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt ] && [ -f $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt.gz ];
            then
                echo "zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt"
                zcat $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt.gz > $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt
            fi

            bash ${github_dir}/GeneNetwork/scripts_per_step/12_GeneNetwork_WebsiteMatrixCreator.sh  -i $f \
                                                                                                    -t $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt \
                                                                                                    -a $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt \
                                                                                                    -z $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt \
                                                                                                    -o $output_dir/11_ImportToWebsite/${n_eigenvectors}/ \
                                                                                                    -g $github_dir/GeneNetwork/ \
                                                                                                    -c $github_dir/GeneNetwork/config_file_templates/ \
                                                                                                    -p $project_dir
        fi

        echo "# $github_dir/GeneNetwork/GeneNetworkBackend/PathwayMatrix/${matrix_name}.matrix.txt" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "# $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "# $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "# $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "# $output_dir/11_ImportToWebsite/*" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
    done
    today=$(date +"%Y-%m-%d")
    if [ ! -f $output_dir/${today}-list_of_genes.txt ]
    then
        echo "extract list of genes"
        awk '{print $1}' $output_dir/10_GeneNetwork_predictions/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt | tail -n+2 > $output_dir/${today}-list_of_genes.${n_eigenvectors}_eigenvectors.txt
    fi
    echo "# $output_dir/list_of_genes.txt" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh

    echo "# $ensg_hgnc_file" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
    echo "# $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt.gz" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
    echo "# ${output_file_step8}" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
    echo "set -e" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
    echo "set -u" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh

    for f in $github_dir/GeneNetwork/GeneNetworkBackend/PathwayMatrix/*${n_eigenvectors}_eigenvectors.terms.txt;
    do
        matrix_name=$(basename ${f%_terms.txt})
        matrix_name=${matrix_name%.terms.txt}
        echo "### $matrix_name ###" >> $output_dir/11_ImportToWebsite/populate_database.sh
        echo "bash data_scripts/populate.sh \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -d /data/$name/ \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -g data/list_of_genes.${n_eigenvectors}_eigenvectors.txt \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -e data/$(basename $ensg_hgnc_file) \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -b ${matrix_name} \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -t data/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonSigTerms.txt \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -i data/${matrix_name}.${n_eigenvectors}_eigenvectors.matrix.txt \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -z data/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.bonSigOnly.txt \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -a data/${matrix_name}.${n_eigenvectors}_eigenvectors.predictions.AUC.bonferonni.txt \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -c data/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt.gz \\" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
        echo "  -n $n_eigenvectors" >> $output_dir/11_ImportToWebsite/populate_database.${n_eigenvectors}_eigenvectors.sh
    done

    echo "Done making all backend files. copy $output_dir/11_ImportToWebsite/populate_database.sh to your local machine to make the webserver"
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -t TMPDIR -o output_dir -p project_dir"
    echo "                    -j jar_dir -g github_dir -v threads -e ensg_hgnc_file"
    echo "                    -d GeneNetworkDir -n name -m mem -z n_eigenvectors [-q qos]"
    echo "  -t      TMPDIR where files will be written during runtime"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output directory where results will be written"
    echo "  -j      Location of eqtl-mapping-pipeline.jar"
    echo "  -g      Location of git cloned https://github.com/npklein/brain_eQTL/ directory"
    echo "  -a      GTF file"
    echo "  -v      Number of threads to use for correlation step and PCA step"
    echo "  -d      GeneNetwork directory (with backend data for predictions"
    echo "  -m      Memory to use for some steps"
    echo "  -n      Name that will be used in output file. Needs to be same as for the previous script!"
    echo "  -z      Number of eigenvectors to use"
    echo "  -q      qos to use when submitting to cluster"
    echo "  -e      file with 1st column ensg ID, 3rd column HNGC symbol"
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
            -o | --output_dir )             shift
                                            output_dir=$1
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
            -v | --threads )                shift
                                            threads=$1
                                            ;;
            -d | --gene_network_dir )       shift
                                            gene_network_dir=$1
                                            ;;
            -m | --mem )                    shift
                                            mem=$1
                                            ;;
            -n | --name )                   shift
                                            name=$1
                                            ;;
            -z | --n_eigenvectors )         shift
                                            n_eigenvectors=$1
                                            ;;
            -q | --qos )                    shift
                                            qos=$1
                                            ;;
            -e | --ensg_hgnc_file )         shift
                                            ensg_hgnc_file=$1
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
    if [ -z "$gene_network_dir" ];
    then
        echo "ERROR: -d/--gene_network_dir not set!"
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
    if [ -z "$name" ];
    then
        echo "ERROR: -n/--name not set!"
        usage
        exit 1;
    fi
    if [ -z "$n_eigenvectors" ];
    then
        echo "ERROR: -z/--n_eigenvectors not set!"
        usage
        exit 1;
    fi
    if [ -z "$ensg_hgnc_file" ];
    then
        echo "ERROR: -e/--ensg_hgnc_file not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
main "$@";






