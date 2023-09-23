#ml Anaconda3
#source /home/umcg-hbrugge/.bashrc
#source activate numpyenv

for tissue in Basalganglia-EUR Cerebellum-EUR Cortex-AFR Cortex-EAS Hippocampus-EUR Spinalcord-EUR
do
        echo $tissue
        mkdir ./output/${tissue}
        python ./splitPerDataset.py \
                ../1-createClusters/output/ \
                /groups/umcg-biogen/tmp04/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-${tissue}-SampleToDataset.txt \
                /groups/umcg-biogen/tmp04/output/2021-FreezeThree/2021-02-18-splicing/2022-01-25-SampleLists-MetaBrainLatest/2020-05-26-${tissue}-SelectedRNASeqSamples.txt \
                ./output/${tissue}/
done

