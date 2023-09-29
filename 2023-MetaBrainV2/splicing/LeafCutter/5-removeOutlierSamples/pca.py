import sys
import gzip
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

matplotlib.use('Agg')

if len(sys.argv) < 3:
    print("Usage: infile outprefix npcs")
    sys.exit(0)

file = sys.argv[1]
outprefix = sys.argv[2]
npcs = int(sys.argv[3])

print("parsing "+file)
# pandas weirdness..
# fh = gzip.open(file,'rt')
# header = fh.readline().strip().split("\t")
# fh.close()

# df = pd.read_csv(file, sep='\t',header = None,skiprows=[0])
df = pd.read_csv(file, sep='\t',index_col=0)
print("Read matrix: {} x {}".format(df.shape[0], df.shape[1]))
# df.set_index('-')
print("Correlating..")
cormat = df.corr()
cormat = cormat.fillna(0)

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
