
mkdir -p ./output/
python createQtlAnnotationLeafcutter2.py \
        "../1-createClusters/output/meta-brain_perind.counts.gz" \
        "./output/annotation_leafcutter.txt.gz" \
	"/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
