from impyute.imputation.cs import mice
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Use mice imputation on table.')
parser.add_argument('input_file', help='input file')
parser.add_argument('output_file', help='outputf file')


args = parser.parse_args()

df = pd.read_csv(args.input_file,sep="\t", index_col=0)   

# start the MICE training
print('start training')
imputed = mice(df)
print('done, write to file')

if args.output_file.endswith('gz'):
    imputed.write_csv(args.output_file,
                        compression='gzip')
else:
    imputed.write_csv(args.output_file)
