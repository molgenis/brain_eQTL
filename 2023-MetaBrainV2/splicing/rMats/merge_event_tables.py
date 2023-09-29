import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description='Merge (gzipped) matrices by column')
parser.add_argument('outfile', help='File to write output to')
parser.add_argument('type', help='type of event')
parser.add_argument('input_files', help='files to merge', nargs='+')

args = parser.parse_args()

def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)


value_per_column_per_row = {}
row_set = set([])
columns = []

extra_cols = [0,1, 2, 3, 4]
if args.type == 'MXE':
    extra_cols = [0, 1,  2, 3, 4, 5]

extra_cols_names = []

extra_info = {}
for f in args.input_files:
    print(f)
    column_index = {}
    with openfile(f,'rt') as input_file:
        header = input_file.readline().rstrip('\n').split('\t')
        for index in extra_cols:
            if header[index] not in extra_cols_names:
                extra_info[header[index]] = {}
                extra_cols_names.append(header[index])
        for index, element in enumerate(header):
            if index not in extra_cols:
                columns.append(element)

        for line in input_file:
            line = line.rstrip('\n').split('\t')
            rowname = line[0]+'_'+line[1]
            row_set.add(rowname)
            if rowname not in value_per_column_per_row:
                value_per_column_per_row[rowname] = {}
            for index, element in enumerate(line):
                value_per_column_per_row[rowname][header[index]] = element
            for index in extra_cols:
                extra_info[header[index]][rowname] = line[index]

sorted_rownames = sorted(row_set)
print('done reading, start writing')
with openfile(args.outfile,'wt') as out:
    out.write('rowname')
    for extra_cols in extra_cols_names:
        out.write('\t'+extra_cols)
    for colname in columns:
        out.write('\t'+colname)
    out.write('\n')
    for rowname in sorted_rownames:
        out.write(rowname)
        for extra_cols in extra_cols_names:
            try:
                out.write('\t'+extra_info[extra_cols][rowname])
            except:
                print(extra_cols)
                print(rowname)
                print(extra_info[extra_cols])
                print(extra_info[extra_cols][rowname])
        for colname in columns:
            if colname in value_per_column_per_row[rowname]:
                out.write('\t'+value_per_column_per_row[rowname][colname])
            else:
                out.write('\t')
        out.write('\n')


