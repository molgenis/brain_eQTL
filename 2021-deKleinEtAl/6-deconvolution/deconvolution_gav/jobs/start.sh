#!/bin/bash

JOB_PREFIX=$1
START=$2
END=$3
QOS_INPUT=$4

QOS_OPTIONS=("leftover" "regular" "priority" "panic" "ds", "dev")
QOS="regular"
if [[ " ${QOS_OPTIONS[@]} " =~ " $QOS_INPUT " ]]; then
    if [[ $QOS_INPUT = "panic" ]]; then
      QOS_INPUT="panic mode"
    fi
    QOS=$QOS_INPUT
fi

JOBS_DIR=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/jobs/$JOB_PREFIX
OUTPUT_DIR=output
LOGFILE=$JOBS_DIR/start.txt

for I in $(seq $START $END); do
	JOBNAME="${JOB_PREFIX}${I}"
	JOBFILE="${JOBS_DIR}/${JOBNAME}.sh"
	OUTFILE="${JOBS_DIR}/${OUTPUT_DIR}/${JOBNAME}.out"
	if [ -f $OUTFILE ]; then
    PANIC=$(cat $OUTFILE | grep 'Panic!!!' | wc -l)
		if [ $PANIC = 1 ]; then
      # echo sbatch --qos=$QOS $JOBFILE
      echo -e "$JOBNAME:\t$(sbatch --qos=$QOS $JOBFILE)" >> $LOGFILE
		fi
	else
	  if [ -f $JOBFILE ]; then
	    # echo sbatch --qos=$QOS $JOBFILE
		  echo -e "$JOBNAME:\t$(sbatch --qos=$QOS $JOBFILE)" >> $LOGFILE
	  else
		  echo -e "error: no such file $JOBFILE" >> $LOGFILE
	  fi
	fi
done