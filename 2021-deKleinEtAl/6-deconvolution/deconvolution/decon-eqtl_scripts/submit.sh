dir=$1
yourfilenames=`ls $dir/*.sh`

for entry in $yourfilenames
do
	echo $entry
	sbatch $entry
done