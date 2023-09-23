import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate leafcutter files.')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('study', help="Study of the cramfiles")
parser.add_argument('--sample_split', help='Character to split file name on and take [0] as sample name', default='.cram')

# refGenome = "/groups/umcg-biogen/tmp01/annotation/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
#refGenome = "/apps/data///ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa"
refGenome = "/groups/umcg-biogen/tmp01/annotation/GeneReference/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa" # Gearshift
regtools = "/groups/umcg-biogen/tmp01/tools/leafcutter/regtools/build/regtools"
# clusteroutput = "/groups/umcg-biogen/tmp04/output/2019-11-08-FreezeTwoDotOne/2020-03-31-leafcutter-output/clusteroutput/"
clusteroutput = "/groups/umcg-biogen/tmp01/output/2021-FreezeThree/2021-02-18-splicing/2023-06-27-WUSTL-Leafcutter/0-Leafcutter/output/"
tmpdir = "/groups/umcg-biogen/tmp01/output/2021-FreezeThree/2021-02-18-splicing/2023-06-27-WUSTL-Leafcutter/0-Leafcutter/tmp/"

args = parser.parse_args()

cram_files = glob.glob(args.cram_base_directory+'*cram', recursive=True)

print('found ',len(cram_files),'cram files')

outdir = args.output_directory
job_base_dir = args.jobs_directory

if not os.path.exists(outdir): 
    os.mkdir(outdir)
if not os.path.exists(job_base_dir):
    os.mkdir(job_base_dir)

def make_jobs(template):
    x = 0
    studies = set([])
    for cram in cram_files:
        if not cram.endswith('.cram') and not cram.endswith('.bam') or '/BPD/' in cram:
            continue
        study = args.study
        studies.add(study)
        jobs_dir = job_base_dir +study+'/'
        if not os.path.exists(jobs_dir):
            os.mkdir(jobs_dir)
        x += 1
        results_dir = outdir+'/'+study
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)

        sample = cram.split('/')[-1].split(args.sample_split)[0]
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACETMP', tmpdir)
        new_template = new_template.replace('REPLACELOG', clusteroutput+sample)
        new_template = new_template.replace('REPLACEJUNCOUT', results_dir+'/'+sample+".junc")
        new_template = new_template.replace('REPLACECRAM', cram)
        new_template = new_template.replace('REPLACEREFGENOME', refGenome)
        new_template = new_template.replace('REPLACEREGTOOLS', regtools)
        new_template = new_template.replace('REPLACEBAM', cram.replace('cram','bam'))

        with open(jobs_dir+'/'+sample+'.sh','w') as out:
            out.write(new_template)
    print('Added: '+','.join(studies))


template = '''#!/bin/bash
#SBATCH --job-name=leafcutter_REPLACENAME
#SBATCH --output=REPLACELOG.out
#SBATCH --error=REPLACELOG.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
#module load SAMtools
module load SAMtools/1.17-GCCcore-11.3.0 # gearshift
export TMPDIR=REPLACETMP
echo $TMPDIR
echo "converting cram to bam"
BAMFILE="$TMPDIR/$(basename REPLACEBAM)"
samtools view -@ 1 -b -T REPLACEREFGENOME REPLACECRAM > $BAMFILE

# test if it is single end or paired end.
# 1. get the header out of the BAM file (samtools view -H). This contains the use command that was used to make BAM file
# 2. get the line out of the header that contains the user command (grep)
# 3. split the line on "user command line" (awk) and print everything that is after
# 4. in this string, count the number of occurences of fastq.gz and fq.gz
FASTQINPUTFILE=$(samtools view -H $BAMFILE | grep "user command line" | awk -F"user command line:" '{ print $2}' | grep -o ".fastq.gz\s\|.fq.gz\s\|.fastq\s\|.fqs\s\ " | wc -l)

t=
if [ "$FASTQINPUTFILE" -eq 1 ];
then
    t="single"
    echo "BAM file is single-end"
elif [ "$FASTQINPUTFILE" -eq 2 ];
then
    t="paired"
    echo "BAM file is paired-end"
else
    echo "ERROR: PAIRED was $PAIRED, should have been 1 or 2";
    exit 1;
fi

# 5. determine read length
READLEN=$(samtools view $BAMFILE |head -n1 |awk '{print $10}'|tr -d "\\n" |wc -m)
echo "Read length: $READLEN"

# index bam
samtools index $BAMFILE

ml purge
#ml CMake/3.10.2-GCC-6.4.0-2.28
ml CMake/3.11.4-GCCcore-7.3.0 # gearshift

# run regtools
REPLACEREGTOOLS junctions extract -s 0 -a 8 -m 50 -M 500000 $BAMFILE -o REPLACEJUNCOUT

touch REPLACEJUNCOUT.finished

if [ $? -eq 0 ];
then
    N=$(wc -l < REPLACEJUNCOUT;)
    if [[ "$N" -eq 1 ]];
    then
        echo "ERROR: only header written (possibly running paired end on single end files or vice versa?)"
        exit 1
    fi
    echo "succes!"
    echo "returncode: $?"
    rm $TMPDIR/$(basename REPLACEBAM)
    rm $TMPDIR/$(basename REPLACEBAM).bai
else
    echo "FAIL!"
    echo "returncode: $?"
    rm $TMPDIR/$(basename REPLACEBAM)
    rm $TMPDIR/$(basename REPLACEBAM).bai
    exit 1;
fi


'''


make_jobs(template)
