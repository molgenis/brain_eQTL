#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=10g
#SBATCH --cpus-per-task=8
#SBATCH -J JOBNAME
#SBATCH -o LOGPREFIX.log
#SBATCH -e LOGPREFIX.err

set -e
set -u

#ml Java/11-LTS
# ml Java/11.0.2
ml Java/11.0.16

# CHROM, BATCHFILE, OUTPREFIX
# EXP, GTE, GENOTYPE
threads=8
nice -n20 java -Xmx8g \
	-Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads \
	-Dmaximum.threads=$threads -Dthread.pool.size=$threads \
        -jar MbQTL-1.3-SNAPSHOT-jar-with-dependencies.jar \
        -m mbqtl \
	-a ANNOTATION \
	-e EXPRESSION \
	-g GTE \
	-v GENOTYPE \
	--chr CHROM \
	--perm 1000 \
	--minobservations 30 \
	-gl BATCHFILE \
	--replacemissinggenotypes \
	-o OUTPREFIX \
        --expgroups GROUPS



