import gzip
import sys
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


if len(sys.argv) < 3:
    print("Usage: metricsfile.txt.gz outfile.png")
    sys.exit()

filters = { "RNASEQ_METRICS_PCT_CODING_BASES": 0.1,
            "ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED": 0.6,
            "STAR_LOG_Uniquely_mapped_reads_PCT": 60
            }

infile  = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile,sep='\t', header=0, index_col=0, )


# assume samples on columns
cols = df.columns
if "RNASEQ_METRICS_PCT_CODING_BASES" not in cols:
    print("Metrics not on columns. Transposing...")
    df = df.transpose()
df["X"] = range(0,df.shape[0])
print(df.head)
df.drop(["INSERT_METRICS_PAIR_ORIENTATION"],axis=1,inplace=True)

# for column in df.columns:
#     print(column)

df = df.astype('float')
print(df[["X","ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED"]] )
print(df[["X","ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED"]].dtypes )

nrcols = 1
nrrows = 3

my_dpi = 100
figheight = 2000
figwidth =600
# plt.figure(figsize=(figheight/my_dpi, figwidth/my_dpi), dpi=my_dpi)


sns.set(rc={'figure.figsize': (nrcols * 10, nrrows * 4)})
sns.set_style("ticks")

# fig, axes = plt.subplots(nrows=nrrows, ncols=nrcols, sharex='col', sharey='row', figsize=(figheight/my_dpi, figwidth/my_dpi))
fig, axes = plt.subplots(nrows=nrrows, ncols=nrcols, sharex='col', sharey='row')

# fig.suptitle("QC metric")
sns.scatterplot(data=df, x="X", y="RNASEQ_METRICS_PCT_CODING_BASES", ax=axes[0])
ax=axes[0]
ax.axhline(filters["RNASEQ_METRICS_PCT_CODING_BASES"], ls='--', color="#D7191C", alpha=0.3)
ax.set_title("RNASEQ_METRICS_PCT_CODING_BASES",
    fontsize=15,
    fontweight='bold')

sns.scatterplot(data=df, x="X", y="ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED", ax=axes[1])
ax=axes[1]
ax.axhline(filters["ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED"], ls='--', color="#D7191C", alpha=0.3)
ax.set_title("ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED",
    fontsize=15,
    fontweight='bold')


sns.scatterplot(data=df, x="X", y="STAR_LOG_Uniquely_mapped_reads_PCT", ax=axes[2])
ax=axes[2]
ax.axhline(filters["STAR_LOG_Uniquely_mapped_reads_PCT"], ls='--', color="#D7191C", alpha=0.3)
ax.set_title("STAR_LOG_Uniquely_mapped_reads_PCT",
    fontsize=15,
    fontweight='bold')
for ax in axes:
    ax.set_ylabel("")
    ax.set_xlabel("Sample #")
    
plt.tight_layout()

fig.savefig(outfile)
plt.close(fig)

print(df.shape)
rslt_df = df.loc[df['RNASEQ_METRICS_PCT_CODING_BASES'] > filters['RNASEQ_METRICS_PCT_CODING_BASES']]
print(rslt_df.shape)
rslt_df = df.loc[df['ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED'] > filters['ALIGNMENT_METRICS_PCT_PF_READS_ALIGNED']]
print(rslt_df.shape)
rslt_df = rslt_df.loc[df['STAR_LOG_Uniquely_mapped_reads_PCT'] > filters['STAR_LOG_Uniquely_mapped_reads_PCT']]
print(rslt_df.shape)



fho = open("passingQC.txt",'w')
for sample in rslt_df.index:
    fho.write(sample+"\n")
fho.close()


