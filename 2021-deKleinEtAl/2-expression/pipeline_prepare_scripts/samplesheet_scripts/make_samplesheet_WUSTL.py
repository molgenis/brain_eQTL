
# loop through all subdirs to find all fastq files
import os.path
import glob
import os
import re
import argparse 

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for CMC.')
parser.add_argument('samplesheet', help='CMC samplesheet from synapse')
parser.add_argument('fastq_dir', help='path to fastq file dir')
parser.add_argument('outdir',help='Directory where output is written',
                             nargs='?',
                             default = 'Public_RNA-seq_QC/samplesheets/')

args = parser.parse_args()

individual_per_sample = {}
samples_per_batch = {}
batch_count = {}

samples = {}

with open(args.samplesheet) as input_file:
    for line in input_file:
        felems = line.strip().split("/")
        orig = felems[-1]
        sample = felems[-1].replace("_R1.fastq.gz","")
        sample = sample.replace("_R2.fastq.gz","")
        sample = sample.replace(".bam","")
        samples[sample] = [args.fastq_dir+'/'+sample+"_R1.fastq.gz",args.fastq_dir+'/'+sample+"_R2.fastq.gz",args.fastq_dir+'/'+sample+".bam"]		
        print("Added sample: "+sample+" --> "+orig)

with open(args.outdir+'samplesheet_WUSTL_RNA.1.txt','w') as out:
    out.write("internalId,project,sampleName,reads1FqGz,reads2FqGz,alignedBamOrCram\n")
    for sample in samples.keys():
        data = samples.get(sample)
        out.write(sample+",WUSTL,"+sample+","+data[0]+","+data[1]+","+data[2]+"\n") 
        print("Written "+sample)
