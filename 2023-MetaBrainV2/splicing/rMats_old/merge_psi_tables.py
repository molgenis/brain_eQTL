# merge psi (percent spliced-in)  tables together
import glob
import argparse
import os

parser = argparse.ArgumentParser(description='Merge rMATs files.')
parser.add_argument('psi_base_dir', help='Base dir which contains all subdirectories from which files need to be merged')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('file_prefix', help='Prefix for the file names. Will create files in output_directory with <prefix>.<event>.txt')
parser.add_argument('gtf', help='GTF file')


args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

print('read gtf')
exons = {}
with open(args.gtf) as input_file:
    for line in input_file:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] != 'exon':
            continue

        exon_location = line[0]+'_'+line[3]+'_'+line[4]

        exon_id = line[8].split('exon_id "')[1].split('"')[0]
        if exon_location not in exons:
            exons[exon_location] = set([])
        exons[exon_location].add(exon_id)
events = ['A3SS','A5SS','MXE','RI','SE']


incl_per_sample_per_event = {'A3SS':{},'A5SS':{},'MXE':{},'RI':{},'SE':{}}
exons_per_event = {}
samples = set([])
genes = {}

def get_exon(chr, start, stop):
    start = str(int(start)+1)
    location = chr+'_'+start+'_'+stop
    try:
        return exons[location]
    except KeyError:
        return []
event_header = {}
samples = []
print('read data');
exon_per_splice_event = {}
for event in events:
    print(event)
    if event not in exons_per_event:
        exons_per_event[event] = set([])
    event_files = glob.glob(args.psi_base_dir+'/**/'+event+'*txt', recursive=True)
    n_event_files = len(event_files)
    print('parsing through',n_event_files,'files')
    for index, f in enumerate(event_files):
        if index % 100 == 0:
            print(str(index)+'/'+str(n_event_files))
        sample = f.split('/')[-2]
        if sample not in samples:
            samples.append(sample)
        with open(f) as input_file:
            header = input_file.readline().strip().split('\t')
            for line in input_file:
                line = line.strip().split('\t')
                gene = line[0]
                chr = line[2]
                splice_event = line[4]+'_'+line[5]+'_'+line[6]+'_'+line[7]+'_'+line[8]+'_'+line[9]
              
                exon_per_splice_event[splice_event] = [get_exon(chr, line[4], line[5]),
                                                      get_exon(chr, line[6], line[7]),
                                                      get_exon(chr, line[8], line[9])]
                if event == 'MXE':
                    exon_per_splice_event[splice_event].append(get_exon(chr, line[10], line[11]))
                    
                genes[splice_event] = gene
                exons_per_event[event].add(splice_event)
                incl_level = line[-1]
                if incl_level == 'NA':
                    incl_level = ''
                else:
                    incl_level = str(round(float(incl_level),4))
                if splice_event not in incl_per_sample_per_event[event]:
                    incl_per_sample_per_event[event][splice_event] = {}
                incl_per_sample_per_event[event][splice_event][sample] = incl_level


for event in events:
    splice_event_sorted = sorted(incl_per_sample_per_event[event])
    outfile = args.output_directory+'/'+args.file_prefix+'.'+event+'.txt'
    with open(outfile, 'w') as out:
        out.write('gene\tsplice_event')
        if event == 'A3SS' or event == 'A5SS':
            out.write('\tlongExon\tshort\tflanking')
        elif event == 'MXE':
            out.write('\t1stExon\t2ndExon\tupstream\tdownstream')
        elif event == 'RI':
            out.write('\triExon\tupstream\tdownstread')
        elif event == 'SE':
            out.write('\texon\tupstream\tdownstream')

        for sample in samples:
            out.write('\t'+sample)
        out.write('\n')

        for splice_event in splice_event_sorted:
            out.write(genes[splice_event]+'\t'+splice_event+'\t')
            out.write(','.join(exon_per_splice_event[splice_event][0])+'\t')
            out.write(','.join(exon_per_splice_event[splice_event][1])+'\t')
            out.write(','.join(exon_per_splice_event[splice_event][2]))
            if event == 'MXE':
                out.write('\t'+','.join(exon_per_splice_event[splice_event][3]))
            for sample in samples:
                try:
                    out.write('\t'+incl_per_sample_per_event[event][splice_event][sample])
                except:
                    out.write('\t')
            out.write('\n')
    print('file written to',outfile)    
