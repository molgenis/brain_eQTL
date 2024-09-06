ml Java
ml plink
ml Python/2.7.12-foss-2015b

dsname=$1
dsdir=$2

snpmapref=/groups/umcg-biogen/tmp03/annotation/SNPMappings-dbsnp151-b38.txt.gz
kgbim=/groups/umcg-biogen/tmp03/annotation/1kgp3v5a-20190312-b38/ALL-merged-common.bim
kgpedmap=/groups/umcg-biogen/tmp03/annotation/1kgp3v5a-20190312-b38/ALL-merged-common


vcfutils=/groups/umcg-biogen/tmp03/tools/VCFUtils.jar
ttvcfconverter=/groups/umcg-biogen/tmp03/tools/ConvertTTToVCF.jar
updatersscript=/groups/umcg-biogen/tmp03/tools/genotypeqc/updaters.py
filteranddedupscript=/groups/umcg-biogen/tmp03/tools/liftover/filteranddedup.py
missingnessplotscript=/groups/umcg-biogen/tmp03/tools/lmiss-hist.Rscript

workdir=$dsdir/qcout/
mkdir -p $workdir
if [ -d $workdir ]; then
	rm $workdir/*
fi

rm $workdir/*

# convert to VCF
java -Xmx4g -jar $ttvcfconverter \
	$dsdir \
	$workdir \
	0.01 0.5 0.0001

#read -p "Press enter to continue1"
# convert to plink
plink --vcf $workdir/ConvertedFromTriTyper.vcf.gz \
	--make-bed --out $workdir/$dsname


# replace SNP ids to 1kg
python $updatersscript $workdir/$dsname.bim $snpmapref

# missingness
plink --bfile $workdir/$dsname \
	--missing \
	--out $workdir/$dsname

#read -p "Press enter to continue2

# R script for missingness
Rscript $missingnessplotscript $workdir/$dsname.lmiss $workdir/$dsname.lmiss.pdf

#read -p "Press enter to continue3"
# heterozygosity
plink --bfile $workdir/$dsname --het --out $workdir/$dsname

#read -p "Press enter to continue4"
# compare to 1kg
python $filteranddedupscript \
       $kgbim \
       $workdir/$dsname.bim \
       $workdir/snpintersect.txt

#read -p "Press enter to continue5"
plink --bfile $workdir/$dsname --extract $workdir/snpintersect.txt --make-bed --out $workdir/$dsname\-1kgsnps

#read -p "Press enter to continue6"
plink --bfile $kgpedmap \
	--extract $workdir/snpintersect.txt --make-bed --out $workdir/1kg-snpintersect

#read -p "Press enter to continue7"
plink --bfile $workdir/1kg-snpintersect \
	--bmerge $workdir/$dsname\-1kgsnps --make-bed --out $workdir/$dsname\-1kgmerged
#read -p "Press enter to continue8"

if [ ! -f $workdir/$dsname\-1kgmerged.bim ]; then
        plink --bfile $workdir/$dsname\-1kgsnps --make-bed \
		--out $workdir/$dsname\-1kgsnps-filter \
		--exclude $workdir/$dsname\-1kgmerged-merge.missnp
#read -p "Press enter to continue9"
        plink --bfile $workdir/1kg-snpintersect --make-bed \
		--out $workdir/1kg-snpintersect-filter \
		--exclude $workdir/$dsname\-1kgmerged-merge.missnp
#read -p "Press enter to continue10"
        plink --bfile $workdir/1kg-snpintersect-filter \
		--bmerge $workdir/$dsname\-1kgsnps-filter \
		--make-bed --out $workdir/$dsname\-1kgmerged
#read -p "Press enter to continue11"


fi

plink --bfile $workdir/$dsname\-1kgmerged --indep-pairwise 50 5 0.2 --out $workdir/$dsname\-1kgmerged
#read -p "Press enter to continue12"
plink --bfile $workdir/$dsname\-1kgmerged --extract $workdir/$dsname\-1kgmerged.prune.in --pca 4 --out $workdir/$dsname\-pca
#read -p "Press enter to continue13"


plink --bfile $workdir/$dsname --indep-pairwise 50 5 0.2 --out $workdir/$dsname\-fulldata

# IBS
plink --bfile $workdir/$dsname \
	--extract $workdir/$dsname\-fulldata.prune.in \
	--recode vcf --out $workdir/$dsname\-pruned
	
gzip -v $workdir/$dsname\-pruned.vcf

#read -p "Press enter to continue14"
java -Xmx20g -jar $vcfutils --similarity \
        -l $workdir/$dsname\-fulldata.prune.in \
        -i $workdir/$dsname\-pruned.vcf.gz \
        -o $workdir/$dsname\-geneticsimilarity

# MDS
plink --bfile $workdir/$dsname --genome \
	--extract $workdir/$dsname\-fulldata.prune.in \
	--out $workdir/$dsname\-genome
	
plink --bfile $workdir/$dsname --read-genome $workdir/$dsname\-genome.genome \
	--extract $workdir/$dsname\-fulldata.prune.in \
	--cluster --mds-plot 4 --out $workdir/$dsname\-mds

#rm $workdir/*.bed
#rm $workdir/*.bim
#rm $workdir/*.fam
#rm $workdir/*.vcf.gz
#rm $workdir/*.ped
#rm $workdir/*.map
#rm $workdir/ConversionLog.txt.gz
#rm $workdir/*.prune.in
#rm $workdir/*.prune.out
#rm $workdir/*.nosex
#rm $workdir/*.log

gzip -v $workdir/*
