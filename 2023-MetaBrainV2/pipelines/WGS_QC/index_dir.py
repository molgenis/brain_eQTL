import os
import argparse
import sys


def index_dir(path, file_substring = None):
    global output_list
    
    # Loop through all entries in a directory. An entry could be a file or another directory
    for entry in os.listdir(path):
        entry_path = os.path.join(path, entry)

        # If the entry is a file
        if os.path.isfile(entry_path):

            # Check if path is a (gzipped) vcf file
            if entry_path.endswith(('.vcf', '.vcf.gz')):

                # If the substring parameter is set, check if the file contains the specified substring
                # print(file_substring)
                if file_substring is not None:
                    if file_substring in entry_path:
                        output_list.append(entry_path)
                # Add path to output list if no substring parameter is set
                else:
                    output_list.append(entry_path)
            

        # If the entry is a directory, recursively call the check_dir function again for that directory
        if os.path.isdir(entry_path):
            index_dir(entry_path, file_substring)


def write_lines_to_file(output_list, out_prefix):
        
    # Write lines to file    
    f = open(f'{out_prefix}.txt', "w")

    for line in output_list:
        f.write(line)

        # Don't write a new line character after the last line
        if output_list.index(line) != len(output_list) - 1:
            f.write('\n')

    f.close()


if __name__ == '__main__':
    # Create argument parsers object and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", help="The directory that needs to be indexed (REQUIRED)")
    parser.add_argument("--file_substring", help="Substring that should be present in all files to select (OPTIONAL)")
    parser.add_argument("--out", help="Prefix of the output file (REQUIRED)")
    args = parser.parse_args()

    # Check if the input directory is a valid directory
    input_directory = args.input_dir
    if input_directory is None or not os.path.isdir(input_directory):
        raise ValueError('The input directory is not a valid directory. Please check the path.')
            
    
    # Check if the expected file substring is not empty
    file_substring = args.file_substring
    if file_substring == '':
        raise ValueError('The file substring is not a valid string.')


    # Check if out parameter is set
    out_prefix = args.out
    if out_prefix is None or out_prefix is '':
        raise ValueError('Please provide a prefix for the output file')


    output_list = []
    index_dir(input_directory, file_substring)
    write_lines_to_file(output_list, out_prefix)
        