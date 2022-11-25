#!/bin/bash

LOG_DIR=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/
JOB_PREFIX=$1
START=$2
END=$3

TOTAL_DONE=0
TOTAL_EQTLS=0

echo -e "JOBNAME\tJOBID\tNODE\tELAPSED\tERRORS\tEQTLS_DONE\tN_EQTLS\tPROGRESS\tPANIC\tFINISHED"
for I in $(seq $START $END); do
	JOBNAME="${JOB_PREFIX}${I}"
	FILENAME="${LOG_DIR}${JOBNAME}.out"
	if [ -f $FILENAME ]; then
		JOBID=$(tail $FILENAME | grep 'Resources consumed' | cut -f5 -d' ')
		NODE=$(tail $FILENAME | grep 'Resources consumed' | cut -f13 -d' ' | sed 's/.$//')
    ELAPSED=$(tail $FILENAME | grep "${JOBID}.batch" | awk -F' ' '{print $2}')
		ERRORS=$(cat $FILENAME | grep -i 'error' | wc -l)
		DONE_EQTLS=$(cat $FILENAME | grep 'Processing' | tail -n 1 | grep -oP '[0-9]{1,2}' | sed -n 1p)
		N_EQTLS=$(cat $FILENAME | grep 'Processing' | tail -n 1 | grep -oP '[0-9]{1,2}' | sed -n 2p)
    PROGRESS=$(cat $FILENAME | grep 'Processing' | tail -n 1 | grep -oP '[0-9]{1,3}%')
    PANIC=$(cat $FILENAME | grep 'Panic!!!' | wc -l)
		FINISHED=$(tail $FILENAME | grep 'Shutting down' | wc -l)
		echo -e "$JOBNAME\t$JOBID\t$NODE\t$ELAPSED\t$ERRORS\t$DONE_EQTLS\t$N_EQTLS\t$PROGRESS\t$PANIC\t$FINISHED"

		TOTAL_DONE=$((TOTAL_DONE + DONE_EQTLS))
		TOTAL_EQTLS=$((TOTAL_EQTLS + N_EQTLS))
	fi
done

TOTAL_PROGRESS=$(echo "scale=2; (${TOTAL_DONE}/${TOTAL_EQTLS})*100" | bc)

echo "-------------------------------------------------------------------------"
echo "Total progress: $TOTAL_DONE/$TOTAL_EQTLS [$TOTAL_PROGRESS%]"