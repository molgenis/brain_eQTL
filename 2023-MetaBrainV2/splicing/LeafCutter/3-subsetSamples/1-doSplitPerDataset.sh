module purge
module load Python/3.6.3-foss-2015b

mkdir -p output

for tissue in Basalganglia-EUR Cerebellum-EUR Cortex-AFR Cortex-EUR Hippocampus-EUR Spinalcord-EUR
do

	python splitPerDataset-v2.py \
		../1-createClusters/output/meta-brain_perind.counts.gz \
		./cramToRNAseqId.txt.gz \
		/groups/umcg-biogen/tmp04/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-${tissue}-SampleToDataset.txt \
		/groups/umcg-biogen/tmp04/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-${tissue}-SelectedRNASeqSamples.txt \
		./output/${tissue} 2 &
	

done

for tissue in Cortex-EAS
do
	python splitPerDataset-v2.py \
		../1-createClusters/output/meta-brain_perind.counts.gz \
		./cramToRNAseqId.txt.gz \
		/groups/umcg-biogen/tmp04/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-${tissue}-SampleToDataset.txt \
		/groups/umcg-biogen/tmp04/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-${tissue}-SelectedRNASeqSamples.txt \
		./output/${tissue} 1 &

done

