for tissue in Basalganglia-EUR Cerebellum-EUR Cortex-AFR Cortex-EAS Hippocampus-EUR Spinalcord-EUR 
do
        echo $tissue
        cd output/$tissue/

        for dataset in $(find * | grep -oP '\K.*(?=_meta)' | sort | uniq);
        do
                echo $dataset
                zcat ${dataset}_meta-brain_perind.counts.gz.phen_chr1 | head -n 1 > ${dataset}_merged_fractions.txt
                
		for chrFile in ${dataset}_meta*;
		do
			zcat $chrFile | tail -n +2 >> ${dataset}_merged_fractions.txt
		done
		gzip ${dataset}_merged_fractions.txt
                echo 'merged'
        done
        cd ../../
done

