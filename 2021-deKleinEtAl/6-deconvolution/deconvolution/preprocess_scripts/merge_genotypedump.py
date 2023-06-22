#!/usr/bin/env python3

import sys
import glob
import os
import pandas as pd

indir = sys.argv[1]


def main():
    print(indir)

    df_list = []
    for i in range(1, 23):
        fpath = os.path.join(indir, "chr{}".format(i), "GenotypeData.txt.gz")
        if not os.path.exists(fpath):
            continue

        df = pd.read_csv(fpath, sep="\t", header=0, index_col=0)
        df_list.append(df)
    df = pd.concat(df_list, axis=0)
    print(df)

    df.to_csv(os.path.join(indir, "GenotypeData.txt.gz"), sep="\t", header=True, index=True, compression="gzip")


if __name__ == '__main__':
    main()
