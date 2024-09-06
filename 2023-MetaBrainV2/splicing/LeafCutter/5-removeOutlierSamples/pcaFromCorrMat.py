import sys
import gzip
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats
import numpy
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from sklearn.decomposition import PCA

matplotlib.use('Agg')

if len(sys.argv) < 3:
    print("Usage: infile outprefix npcs [threads:1] [subsample:1.0]")
    sys.exit(0)

file = sys.argv[1]
outprefix = sys.argv[2]
npcs = int(sys.argv[3])
threads = 1
subsample = 1
if len(sys.argv) > 3:
    threads = int(sys.argv[4])
if len(sys.argv) > 4:
    subsample = float(sys.argv[5])


def corr(x, y):
    xtmp = []
    ytmp = []

    for i in range(0,len(x)):
        if not numpy.isnan(x[i]) and not numpy.isnan(y[i]):
            xtmp.append(x[i])
            ytmp.append(y[i])

    if len(xtmp) == 0:
        return 0
    result = scipy.stats.pearsonr(xtmp,ytmp)
    return result.statistic

def processI(i, df):
    x = df.iloc[:, i].values
    corrs = [0] *  len(df.columns)
    for j in range(i+1,len(df.columns)):
        y = df.iloc[:, j].values
    
        r = corr(x,y)
        # print(f"{i}\t{j}\t{r}")
        corrs[j] = r
    return corrs


print("parsing "+file)
sys.stdout.flush()
# pandas weirdness..
# fh = gzip.open(file,'rt')
# header = fh.readline().strip().split("\t")
# fh.close()

# df = pd.read_csv(file, sep='\t',header = None,skiprows=[0])
df = pd.read_csv(file, sep='\t',index_col=0)
print("Read matrix: {} x {}".format(df.shape[0], df.shape[1]))




cormat = df
print(cormat)

print("Performing decomposition.")
pca = PCA(n_components=npcs)
pca.fit(cormat)
components = pca.components_
explVar = pca.explained_variance_ratio_
#			print(explVar)
pcadf = pd.DataFrame(components)
pcadf.columns = cormat.columns
#			pcadf.index.name = 'Component'
pcadf = pcadf.transpose()
colnames = []
for comp in range(1,npcs+1):
    colnames.append("PC"+str(comp))
pcadf.columns = colnames
#			print(pcadf)

#pcadfmelt = pcadf.melt()
#			print(pcadfmelt)
pcadf.to_csv(outprefix+"_PCs.txt", sep='\t')
fig, ax = plt.subplots()
sns.scatterplot(data=pcadf, x="PC1", y="PC2")
fig.savefig(outprefix+"_PC1and2.png")
plt.close(fig)

sns.scatterplot(data=pcadf, x="PC1", y="PC4")
fig.savefig(outprefix+"_PC1and4.png")
plt.close(fig)
