#module load Python/3.10.4-GCCcore-11.3.0

module load Python/2.7.16-GCCcore-11.3.0-bare
python leafcutter_cluster_regtools_fix.py \
    -j junction_files.txt \
    -o output/outprefix \
    -m 50 -l 500000 \
    –checkchrom


