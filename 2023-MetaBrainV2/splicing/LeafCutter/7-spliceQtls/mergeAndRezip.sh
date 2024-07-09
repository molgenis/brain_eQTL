
ml Python
ml Java
ml HTSlib

indir=$1
threads=2
for chr in {1..22}
do

nice -n20 python mergeAllEffects.py \
                ${indir} \
                ${chr}

nice -n20 java -Xmx4g \
        -Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads \
        -Dmaximum.threads=$threads -Dthread.pool.size=$threads \
        -jar MbQTL-1.4.5-SNAPSHOT-jar-with-dependencies.jar \
        -m sortfile \
        --input ${indir}/chr${chr}-AllEffects.txt.gz \
        --out ${indir}/chr${chr}-AllEffects-sorted.txt.gz

nice -n20 zcat ${indir}/chr${chr}-AllEffects-sorted.txt.gz | nice -n20 bgzip -c > ${indir}/chr${chr}-AllEffects-sorted.bgz.txt.gz
nice -n20 tabix -f -s 3 -b 4 -e 4 -S 1 ${indir}/chr${chr}-AllEffects-sorted.bgz.txt.gz

done
