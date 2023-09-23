\


for tissue in Basalganglia-EUR Cerebellum-EUR Cortex-AFR Cortex-EAS Cortex-EUR Hippocampus-EUR Spinalcord-EUR 
do

	python count.py \
		${tissue} \
		./output/${tissue}-PSI-unfiltered.txt.gz \
		../2-annotation/output/annotation_leafcutter.txt.gz
done
