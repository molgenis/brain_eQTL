#!/bin/bash

JOBS_DIR=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/
LOGFILE=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/start.txt
JOB_PREFIX=$1
START=$2
END=$3

echo $(date) > $LOGFILE

for I in $(seq $START $END); do
	JOBNAME="${JOB_PREFIX}${I}"
	FILENAME="${JOBS_DIR}${JOBNAME}.sh"
	if [ -f $FILENAME ]; then
		echo -e "$JOBNAME:\t$(sbatch $FILENAME)" >> $LOGFILE
	else
		echo -e "error: no such file $FILENAME" >> $LOGFILE
	fi
done