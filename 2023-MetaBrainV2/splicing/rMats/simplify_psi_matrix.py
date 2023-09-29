import glob
import argparse
import gzip
import os
import multiprocessing
import datetime
import sys

parser = argparse.ArgumentParser(description='Simplify psi matrix')
parser.add_argument('psi_dir', help='Directory with psi files')
parser.add_argument('outdir', help='outdir to write to')
parser.add_argument('--multiprocess', help='if multiprocessing should be used', action='store_true')
    
args = parser.parse_args()



if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def simplify_matrix(f):

    none_samples = set(['rowname','gene','splice_event','short','longExon','flanking','downstream','upstream','exon','riExon','downstread','1stExon','2ndExon'])

    prev_identifier = None
    print(f)
    sys.stdout.flush()
    with openfile(f,'rt') as input_file, openfile(args.outdir+'/'+f.split('/')[-1].replace('.txt','.simplified.txt'),'wt') as out:
        out.write('exon')
        header = input_file.readline().strip('\n').split('\t')
        samples = []
        for element in header:
            if element not in none_samples:
                out.write('\t'+element)
                samples.append(element)
        out.write('\n')
        current_identifier_psi_values = {}
        for line_index, line in enumerate(input_file):
            if line_index % 1000 == 0:
                currentDT = datetime.datetime.now()
                print(str(currentDT)+' - '+f,'-',line_index,'lines processed')
                sys.stdout.flush()
            line = line.strip('\n').split('\t')
            if 'MXE' in f:
                identifier = line[3]+'_'+line[4]
            else:
                identifier = line[3]
            if prev_identifier and identifier != prev_identifier:
                out.write(prev_identifier)
                for sample in samples:
                    if len(current_identifier_psi_values[sample]) > 0:
                        psi = 0
                        for value in current_identifier_psi_values[sample]:
                            psi += float(value)
                        psi /= len(current_identifier_psi_values[sample])
                        out.write('\t'+str(psi))
                    else:
                        out.write('\t')
                out.write('\n')

                current_identifier_psi_values = {}

            for index, element in enumerate(line):
                if header[index] in none_samples:
                    continue
                else:
                    if header[index] not in current_identifier_psi_values:
                        current_identifier_psi_values[header[index]] = set([])
                    if element != '':
                        current_identifier_psi_values[header[index]].add(element)
            prev_identifier = identifier

        # have to do here as well to include last identifier
        out.write(identifier)
        for sample in samples:
            if len(current_identifier_psi_values[sample]) > 0:
                psi = 0
                for value in current_identifier_psi_values[sample]:
                    psi += float(value)
                    psi /= len(current_identifier_psi_values[sample])
                    out.write('\t'+str(psi))
            else:
               out.write('\t')
        out.write('\n')
    return

input_files = glob.glob(args.psi_dir+'/*gz')
input_files_to_run = []
for f in input_files:
    outfile = args.outdir+'/'+f.split('/')[-1].replace('.txt','.simplified.txt')
#    if not os.path.isfile(outfile):
    input_files_to_run.append(f)

if args.multiprocess:
    pool = multiprocessing.Pool(len(input_files_to_run))
    pool.map(simplify_matrix, input_files_to_run)
    pool.close()
else:
    for f in input_files_to_run:
        simplify_matrix(f)
