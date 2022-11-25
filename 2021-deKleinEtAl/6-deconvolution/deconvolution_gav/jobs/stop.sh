#!/bin/bash

JOB_PREFIX=$1
LOGFILE=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/jobs/$JOB_PREFIX/start.txt

re='^[0-9]+$'
while IFS= read -r LINE; do
	JOBID=$(echo $LINE | cut -f5 -d' ')
	if [[ $JOBID =~ $re ]] ; then
		echo "scancel $JOBID"
		scancel $JOBID
	fi
done < "$LOGFILE"