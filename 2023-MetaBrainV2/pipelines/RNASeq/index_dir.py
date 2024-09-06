import os
import argparse
import sys

FASTQ_EXTENSIONS = ('fq', 'fastq', 'fq.gz', 'fastq.gz')

def index_dir(path, samples):
    '''Function that recursively loops through a directory and indexes all BAM and FASTQ'''
    global output_dict
    
    # Loop through all entries in a directory. And entry could be a file or another directory
    for entry in os.listdir(path):
        entry_path = os.path.join(path, entry)

        # If the entry is a file
        if os.path.isfile(entry_path):

            # If the file is gzipped, remove the gz extension
            if entry.split('.')[-1] == 'gz':
                file = '.'.join(entry.split('.')[:-1])
            else:
                file = entry

            # Get file name without extension
            extension = file.split('.')[-1]
            file_name = '.'.join(file.split('.')[:-1])

            # If the file is a .bam file, add it to the output dict
            if extension == 'bam' and (samples is None or any(substring in entry_path for substring in samples)):
                output_dict[file_name] = entry_path

            # If the file is a fastq file
            if extension in FASTQ_EXTENSIONS:

                # If the file name ends with '_1', it is assumed that it is one of two paired end files
                if file_name.endswith('_1') or file_name.lower().endswith('_r1') or file_name.endswith('_R1_001') or file_name.lower().endswith('.r1'):
                    if file_name.lower().endswith('.r1'):
                        sample_name = file_name[:-3]
                    
                    if file_name.endswith('_1'):
                        sample_name = file_name[:-2]
                    
                    if file_name.lower().endswith('_r1'):
                        sample_name = file_name[:-3]

                    if file_name.endswith('_R1_001'):
                        sample_name = file_name[:-7]
                    
                    
                    if sample_name in output_dict and (samples is None or any(substring in entry_path for substring in samples)):
                        output_dict[sample_name]['1'] = entry_path
                    else:
                        output_dict[sample_name] = {'1': entry_path}


                # If the file name ends with '_2', it is assumed that it is one of two paired end files
                elif file_name.endswith('_2') or file_name.lower().endswith('_r2') or file_name.endswith('_R2_001') or file_name.lower().endswith('.r2'):
                    if file_name.lower().endswith('.r2'):
                        sample_name = file_name[:-3]
                    
                    if file_name.endswith('_2'):
                        sample_name = file_name[:-2]
                    
                    if file_name.lower().endswith('_r2'):
                        sample_name = file_name[:-3]
                    if file_name.endswith('_R2_001'):
                        sample_name = file_name[:-7]

                    if sample_name in output_dict and (samples is None or any(substring in entry_path for substring in samples)):
                        output_dict[sample_name]['2'] = entry_path
                    else:
                        output_dict[sample_name] = {'2': entry_path}
                
                # If the filename does not end with '_1'/'_r1 or '_2'/'_r2', it is assumed to be single end and added to the output dict
                else:
                    if (samples is None or  any(substring in entry_path for substring in samples)):
                        output_dict[file_name] = entry_path

        # If the entry is a directory, recursively call the check_dir function again for that directory
        if os.path.isdir(entry_path):
            index_dir(entry_path, samples)


def write_lines_to_file(output_dict, out_file):
    '''Function that writes all paths indexed by the index_dir function to the out file'''
    output_list = []

    # Create a list of lines to write to file
    for k, v in output_dict.items():
        if type(v) == str:
            output_list.append(k.split('.')[0] + ',' + v)
        if type(v) == dict and '1' in v and '2' in v:
            output_list.append(k + ',' + v['1'] + ';' + v['2'])
        
    # Write lines to file    
    f = open(f'{out_file}.txt', 'w')

    for line in output_list:
        f.write(line)

        # Don't write a new line character after the last line
        if output_list.index(line) != len(output_list) - 1:
            f.write('\n')

    f.close()

def chunk_file(out_file, chunk_size):
    fh = open(f'{out_file}.txt','r')
    lines = fh.readlines()
    fh.close()
    ctr = 0 # line ctr
    wctr = 0 # lines written counter
    cctr = 1 # chunk counter
    fho = open(f'{out_file}_{cctr}.txt','w')
    while(ctr < len(lines)):
        if ctr % chunk_size == 0 and wctr > 0: # if ctr == 0, this will start a new chunk, unless wctr > 0
            fho.close()
            cctr += 1
            fho = open(f'{out_file}_{cctr}.txt','w')
            wctr = 0
        fho.write(lines[ctr])
        wctr += 1
        ctr += 1
    if wctr > 0:
       fho.close()

if __name__ == '__main__':
    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', help='The directory that needs to be indexed (REQUIRED)')
    parser.add_argument('--sample_file', help='''A text file containing a list of samples that need to be indexed. 
                        The sample names should be split in lines. If no sample list is provided, all samples will be indexed (OPTIONAL)''')
    parser.add_argument('--out', help='Prefix of the outfile.')
    args = parser.parse_args()

    # Check if the input directory is a valid directory
    input_directory = args.input_dir
    if not os.path.isdir(input_directory):
        raise ValueError('The input directory is not a valid directory. Please check the path.')
            
    # Get the samples that need indexing from the sample file
    if args.sample_file:
        sample_file = args.sample_file

        try:
            with open(sample_file, 'r') as sample_file:
                samples = sample_file.readlines()
                samples = [sample.strip() for sample in samples]
        except OSError:
            print('An error occured while reading the sample file')
            sys.exit()
    else:
        samples = None

    
    # Check if out parameter is set
    out_file = args.out
    if not out_file:
        raise ValueError('Please provide a prefix for the output file')
    
    output_dict = {}
    index_dir(input_directory, samples)
    write_lines_to_file(output_dict, out_file)
    chunk_size = 300
    chunk_file(out_file,chunk_size)
