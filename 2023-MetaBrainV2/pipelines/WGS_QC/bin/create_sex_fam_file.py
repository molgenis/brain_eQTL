#!/opt/conda/envs/eQTLGenPopAssign/bin/python

import sys

def read_file(file_path):
    with open(file_path, 'r') as file:
        lines = [l.strip() for l in file.readlines()]

    return lines

def create_sex_fam_file(fam_lines, sex_lines, out_prefix):
    sex_sample_dict = {l.split(',')[0]:l.split(',')[1] for l in sex_lines}
    sex_trans_dict = {'M': '1', 'F': '2'}

    with open(f'{out_prefix}.fam', 'w') as out_fam:
        for line in fam_lines:
            line_parts = line.split('\t')
            sample_sex = sex_sample_dict[line_parts[1]]
            line_parts[4] = sex_trans_dict[sample_sex]
            out_line = '\t'.join(line_parts)
            out_fam.write(out_line + '\n')


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("Usage: python create_sex_fam_file.py input.fam sex_file.txt out")
        print("Input fam file")
        print("Sex file")
        print("Prefix for the out file")
        sys.exit(-1)


    fam_path = sys.argv[1]
    sex_path = sys.argv[2]
    out_prefix = sys.argv[3]

    fam_lines = read_file(fam_path)
    sex_lines = read_file(sex_path)
    create_sex_fam_file(fam_lines, sex_lines, out_prefix)