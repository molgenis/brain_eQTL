ml Java
ml Python


if [ "$#" -ne 3 ]; then
    echo "Usage: dumpdir outprefix fdrreference-FDR0.05.txt.gz"
    exit
fi

dumpdir=$1
outprefix=$2
fdrref=$3

#if [ 1 -eq 2 ]; then
# calculate MAF
java -Xmx10g -jar /groups/umcg-biogen/tmp01/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar --mode util --getmaffromqclog \
        --in $dumpdir/SNPQCLog.txt.gz \
        --out $dumpdir/SNPQCLog-MAF.txt.gz

# sort by Z
java -Xmx10g -jar /groups/umcg-biogen/tmp01/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar --mode util --sortfile \
	--in $dumpdir/eQTLDump.txt.gz \
	--out $dumpdir/eQTLDump-sortZ.txt.gz

# apply FDR
java -Xmx10g -jar /groups/umcg-biogen/tmp01/tools/BinaryMetaAnalyzer-1.0.14-SNAPSHOT-jar-with-dependencies.jar --applyfdr \
	--in $dumpdir/eQTLDump-sortZ.txt.gz \
	--out $dumpdir/eQTLDump-sortZ-FDR \
	--eqtls $fdrref


mkdir $dumpdir/splitByChr/
# sort by position
java -Xmx10g -jar /groups/umcg-biogen/tmp01/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar --mode util --sorteqtlfilebypos \
        --in $dumpdir/eQTLDump-sortZ-FDR.txt.gz \
        --out $dumpdir/eQTLDump-sortPos-FDR.txt.gz

#fi

# split by chr
java -Xmx10g -jar /groups/umcg-biogen/tmp01/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar --mode util --spliteqtlfilebychr \
        --in $dumpdir/eQTLDump-sortPos-FDR.txt.gz \
        --out $dumpdir/splitByChr/$outprefix

#fi


# reformat biogen
mkdir -pv $dumpdir/biogenformat/
mkdir -pv $dumpdir/smr/
mkdir -pv $dumpdir/biogenformat/fdr005/

