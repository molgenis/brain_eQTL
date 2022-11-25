module load Java
#module load parallel
module load Python
module load BWA

if [ $# -lt 2 ]
  then
    echo "Usage: indir threads"
    exit 1
fi

indir=$1
threads=$2
dateprefix=$(date '+%Y-%m-%d')

# calculate fdr if not done yet
if [ ! -f "$indir/eQTLsFDR0.05.txt.gz" ]; then
	lnct=$(zcat $indir/eQTLs.txt.gz | wc -l)
	lnct=$(($lnct - 1))
	java -Xmx10g -jar /groups/umcg-biogen/tmp01/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar \
		--mode util \
		--fdr --fdrmethod full \
		--in $indir \
		--threshold 0.05 \
		--nreqtls $lnct \
		--perm 10
fi

# select combos
zcat $indir/eQTLsFDR0.05.txt.gz | awk -v OFS='\t' '{print $2,$3,$4,$5}' > $indir/eQTLsFDR0.05-combosForCrossMap.txt

workdir="$indir/crossmap/"
readdir="$workdir/reads-10mbwindow-35bp/"
mkdir -pv $readdir

# make sequences
java -Xmx8g -jar /groups/umcg-biogen/tmp01/tools/SequenceCrossmap.jar \
        makeseq \
        $indir/eQTLsFDR0.05-combosForCrossMap.txt \
	/groups/umcg-biogen/tmp01/annotation/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa \
        /groups/umcg-biogen/tmp01/annotation/gencode/gencode.v32.primary_assembly.annotation.collapsedGenes.gtf.gz \
        $readdir \
        10000000 35 2 true false
# exit 1

# index
#parallel -j $threads < "$readdir/index.sh"
bash "$readdir/index.sh"

# align
#parallel -j $threads < "$readdir/align.sh"
bash < "$readdir/align.sh"

# convert to samse
#parallel -j $threads < "$readdir/samse.sh"
bash "$readdir/samse.sh"

# quantify crossmapping reads per eqtl
java -Xmx8g -jar /groups/umcg-biogen/tmp01/tools/SequenceCrossmap.jar \
        quantifyIndividualEQTL \
        $indir/eQTLsFDR0.05.txt.gz \
        $readdir \
        "$workdir/$dateprefix-eQTLCrossMap-5mb-35bp.txt.gz"
		
# cleanup: remove reads and alignments
rm -rv $readdir

# select non-crossmapping eqtls
python /groups/umcg-biogen/tmp01/tools/brain_eQTL/eqtlfile/transcrossmap-selectNonCrossMapping.py \
		$indir/eQTLs.txt.gz \
		"$workdir/$dateprefix-eQTLCrossMap-5mb-35bp.txt.gz" \
		0.05 \
		$indir/eQTLs-crossMappingEQTLsRemoved.txt.gz
		
# copy files to temporary directory
fdrdir="$indir/fdrtmp/"
mkdir -p $fdrdir
cp $indir/eQTLs-crossMappingEQTLsRemoved.txt.gz $fdrdir/eQTLs.txt.gz
cp $indir/PermutedEQTLsPermutationRound* $fdrdir

lnct=$(zcat $fdrdir/eQTLs.txt.gz | wc -l)
lnct=$(($lnct - 1))

# repeat FDR calculation
java -Xmx10g -jar /groups/umcg-biogen/tmp01/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar \
	--mode util \
	--fdr \
	--fdrmethod full \
	--skipqqplot \
	--in $fdrdir \
	--threshold 0.05 \
	--perm 10 \
	--nreqtls $lnct

# copy result back to indir
cp $fdrdir/eQTLsFDR0.05.txt.gz $indir/eQTLs-crossMappingEQTLsRemoved-FDR0.05.txt.gz

# cleanup: remove temparary FDR dir
rm -rv $fdrdir

