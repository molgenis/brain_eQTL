import glob
import sys
import gzip



if len(sys.argv) < 3:
    print("Usage: indir outfile.txt.gz")
    sys.exit(-1)

indir = sys.argv[1]
outfile = sys.argv[2]

def getfh(file):
    if file.endswith(".gz"):
        return gzip.open(file,'rt')
    else:
        return open(file,'r')

def getfho(file):
    if file.endswith(".gz"):
        return gzip.open(file,'wt')
    else:
        return open(file,'w')

def getSample(file,removeStr):
    sample = file.split("/")[-1]
    sample = sample.replace(".gz","")
    sample = sample.replace(removeStr,"")
    return sample    

def getSamples(files,removeStr):
    print("{} files".format(len(files)))
    samples = set()
    for file in files:
        sample = getSample(file,removeStr)
        samples.add(sample)
    nrsamples = len(samples)
    print(f"{nrsamples} samples with {removeStr}")
    return samples

def toSortedArr(s):
    d = []
    for i in s:
        d.append(i)
    d.sort()
    return d

    
# inventorize samples
files = glob.glob(indir+"/collectMultipleMetrics_QC/*.multiplemetrics.alignment_summary_metrics")
samples = getSamples(files,".multiplemetrics.alignment_summary_metrics")
print("{} samples loaded sofar ".format(len(samples)))
files = glob.glob(indir+"/collectMultipleMetrics_QC/*.multiplemetrics.insert_size_metrics")
samples.update(getSamples(files,".multiplemetrics.insert_size_metrics"))
print("{} samples loaded sofar ".format(len(samples)))
files = glob.glob(indir+"/collectRnaSeqMetrics_QC/*.rna_metrics.log")
samples.update(getSamples(files,".rna_metrics.log"))
print("{} samples loaded sofar ".format(len(samples)))
files = glob.glob(indir+"/star/*.ReadsPerGene.out.tab")
samples.update(getSamples(files,".ReadsPerGene.out.tab"))
print("{} samples loaded sofar ".format(len(samples)))
files = glob.glob(indir+"/star/*.Log.final.out")
samples.update(getSamples(files,".Log.final.out"))
print("{} samples loaded total ".format(len(samples)))



# collection bins
data = {}
metrics = set()

## collectMultipleMetrics_QC
# alignment_summary_metrics - skip lines w/ #, skip empty lines
# header line starts with CATEGORY
# line with PAIR has average numbers
files = glob.glob(indir+"/collectMultipleMetrics_QC/*.multiplemetrics.alignment_summary_metrics")
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,".multiplemetrics.alignment_summary_metrics")
    fh = getfh(file)
    header = []
    sampledata = {}

    for line in fh:
        line = line.strip()
        if line.startswith("CATEGORY"):
            elems = line.strip().split("\t")
            for i in range(1,len(elems)):
                col = elems[i]
                if col != "SAMPLE" and col != "LIBRARY" and col != "READ_GROUP":
                    header.append("ALIGNMENT_METRICS_"+elems[i])
                    metrics.add("ALIGNMENT_METRICS_"+elems[i])
        elif line.startswith("PAIR"):
            elems = line.strip().split("\t")
            for i in range(1,len(elems)):
                sampledata[header[i-1]] = elems[i]
            break # no need to read the rest
    fh.close()    
    data[sample] = sampledata        

# insert_size_metrics - skip lines w/ #, skip empty lines; first line begins with; second line has stats 
# MEDIAN_INSERT_SIZE
files = glob.glob(indir+"/collectMultipleMetrics_QC/*.multiplemetrics.insert_size_metrics")
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,".multiplemetrics.insert_size_metrics")
    fh = getfh(file)
    header = []
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}

    for line in fh:
        line = line.strip()
        if line.startswith("MEDIAN_INSERT_SIZE"):
            elems = line.strip().split("\t")
            for i in range(0,len(elems)):
                col = elems[i]
                if col != "SAMPLE" and col != "LIBRARY" and col != "READ_GROUP":
                    dta = sampledata.get(col)
                    col = "INSERT_METRICS_"+col
                    header.append(col)
                    metrics.add(col)

            dataline = fh.readline()
            elems = dataline.strip().split("\t")
            for i in range(0,len(elems)):
                sampledata[header[i]] = elems[i]
            break
    fh.close()
    data[sample] = sampledata

## collectRnaSeqMetrics_QC
# *.rna_metrics.log
# header starts with ## METRICS CLASS; next line has column names, starting with PF_BASES
# line afterwards values
files = glob.glob(indir+"/collectRnaSeqMetrics_QC/*.rna_metrics.log")
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,".rna_metrics.log")
    fh = getfh(file)
    header = []
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}

    for line in fh:
        line = line.strip()
        if line.startswith("PF_BASES"):
            elems = line.strip().split("\t")
            for i in range(0,len(elems)):
                col = elems[i]
                if col != "SAMPLE" and col != "LIBRARY" and col != "READ_GROUP":
                    dta = sampledata.get(col)
                    col = "RNASEQ_METRICS_"+col
                    header.append(col)
                    metrics.add(col)

            dataline = fh.readline()
            elems = dataline.strip().split("\t")

            for i in range(0,len(elems)):
                sampledata[header[i]] = elems[i]
            break
    fh.close()
    data[sample] = sampledata

# STAR count metrics
files = glob.glob(indir+"/star/*.ReadsPerGene.out.tab")
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,".ReadsPerGene.out.tab")
    samples.add(sample)
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}
    fh = getfh(file)
    mapped = [0,0,0]
    
    for line in fh:
        elems = line.strip().split("\t")
        phen = elems[0]
        if not phen.startswith("ENSG"):
            for i in range(1,4):
                ct = elems[i]
                phenoname = elems[0]
                if i == 1:
                    phenoname = "STAR_TAB_"+phenoname + "_sum"
                elif i == 2:
                    phenoname = "STAR_TAB_"+phenoname + "_strandA"
                else:
                    phenoname = "STAR_TAB_"+phenoname + "_strandB"
                sampledata[phenoname] = ct
                metrics.add(phenoname)
        else:
            # sum up the counts for each column
            for i in range(1,4):
                ct = int(elems[i])
                mapped[i-1] += ct
    for i in range(0,len(mapped)): 
            if i == 0:
                phenoname = "STAR_TAB_N_mapped_sum"
            elif i == 1:
                phenoname = "STAR_TAB_N_mapped_strandA"
            else:
                phenoname = "STAR_TAB_N_mapped_strandB"
            sampledata[phenoname] = str(mapped[i])
            metrics.add(phenoname)

    fh.close()
    data[sample] = sampledata


# STAR final logs
allowedp = [
    "Number_of_input_reads",
    "Average_input_read_length",
    "Uniquely_mapped_reads_number",
    "Uniquely_mapped_reads_PCT",
    "Average_mapped_length",
    "Number_of_splices_Total",
    "Number_of_splices_Annotated_(sjdb)",
    "Number_of_splices_GT/AG",
    "Number_of_splices_GC/AG",
    "Number_of_splices_AT/AC",
    "Number_of_splices_Non-canonical",
    "Mismatch_rate_per_base_PCT",
    "Deletion_rate_per_base",
    "Deletion_average_length",
    "Insertion_rate_per_base",
    "Insertion_average_length",
    "Number_of_reads_mapped_to_multiple_loci",
    "PCT_of_reads_mapped_to_multiple_loci",
    "Number_of_reads_mapped_to_too_many_loci",
    "PCT_of_reads_mapped_to_too_many_loci",
    "PCT_of_reads_unmapped_too_many_mismatches",
    "PCT_of_reads_unmapped_too_short",
    "PCT_of_reads_unmapped_other",
    "Number_of_chimeric_reads",
    "PCT_of_chimeric_reads"
]
allowedpset = set()
for p in allowedp:
    allowedpset.add(p)

files = glob.glob(indir+"/star/*.Log.final.out")
for file in files:
    print("Parsing: "+file)
    sample = getSample(file,".Log.final.out")
    samples.add(sample)
    sampledata = data.get(sample)
    if sampledata is None:
        sampledata = {}

    fh = getfh(file)
    for line in fh:
        line = line.strip()
        if "|" in line:
            line = line.replace("|","").strip()
            
            elems = line.split("\t")
            if len(elems) == 2:
                p = elems[0].strip().replace(" ","_")
                p = p.replace(":","")
                p = p.replace(",","")
                p = p.replace("%","PCT")
                # print(p)
                if p in allowedpset:
                    v = elems[1].replace("%","")
                    p = "STAR_LOG_"+p
                    sampledata[p] = v
                    metrics.add(p)
    fh.close()
    # sys.exit()
    data[sample] = sampledata

print("{} metrics loaded ".format(len(metrics)))
print("{} samples loaded ".format(len(samples)))

metrics = toSortedArr(metrics)
samples = toSortedArr(samples)

# for sample in samples:
#     print("{}".format(sample))

print("Writing: "+outfile)
fho = getfho(outfile)
header = "-\t"+"\t".join(samples)+"\n"
fho.write(header)
for metric in metrics:
    outln = metric
    for sample in samples:
        sampledata = data.get(sample)
        if sample is None:
            outln+="\tnan"
        else:
            v = sampledata.get(metric)
            if v is None:
                print("metric {} is missing for sample: {}".format(metric, sample))
                # sys.exit(-1)
                outln +="\tnan"
            else:
                outln +="\t"+v
    
    fho.write(outln+"\n")
fho.close()