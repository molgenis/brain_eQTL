import glob
# import rglob
import sys
import os
import subprocess

template = "readmetemplate31.txt"
treecmd = "tree -L 1  --dirsfirst  -F |sed 's/*//g' "

if len(sys.argv) < 3:
    print("Usage: indir template.txt")
    sys.exit()

indir = sys.argv[1]
template = sys.argv[2]



def fast_scandir(dirname):
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders

def generate(dir):
    print(dir)
    tree = subprocess.run(['tree', '-L','1','--dirsfirst','-F',dir], stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")
    print()
    tree = tree[0:len(tree)-3]
    foldername = dir.split("/")[-1]
    relativefoldername = dir.replace(indir,"")
    for i in range(0,len(tree)):
        line = tree[i]
        line = line.replace(dir,".")
        line = line.replace('*',"")
        line = line + " - "
        tree[i] = line
        # print(line)
    tree = "\n".join(tree)

    readmeout = tpl.replace("TREEMAP",tree)
    readmeout = readmeout.replace("FOLDER",relativefoldername)
    # print(readmeout)
    outf = dir+"/README.txt"
    print(outf)
    fho = open(outf,'w')
    fho.write(readmeout)
    fho.close()
    # print(output)
    # sys.exit()

# for dir in os.walk(indir):
#     print(dir)
fh = open(template,'r')
tpl = fh.readlines()
for i in range(0,len(tpl)):
    tpl[i] = tpl[i].strip()
tpl = "\n".join(tpl)
fh.close()

from glob import glob
dirs = fast_scandir(indir)
for dir in dirs:
    generate(dir)

generate(indir)