from datetime import datetime
import gzip
import wget
from time import strftime
import os
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description='Make GO matrices.')
parser.add_argument('uniprot_ID_mapping',
                    help='idmapping_selected.tab.gz from uniprot.org')
parser.add_argument('ordered_gene_list',
                    help='List with ordered gene IDs')



args = parser.parse_args()

today = datetime.now().strftime("%Y-%m-%d") 

if not os.path.exists(today+'-goa_human.gaf.gz'):
    print('Getting latest go annotation')
    url = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
    wget.download(url, out=today+'-goa_human.gaf.gz')


uniprot_to_geneID = {}
with gzip.open(args.uniprot_ID_mapping, 'rt') as input_file:
    for line in input_file:
        line = line.strip('\n').split('\t')
        uniprot_id = line[0]
        ensembl_ids = line[1].split(';')
        ensembl_ids = [x.strip() for x in ensembl_ids]
        uniprot_to_geneID[uniprot_id] = ensembl_ids


with gzip.open(today+'-goa_human.gaf.gz','rt') as input_file:
    for line in input_file:
        line = line.strip()
        if 'Date Generated by GO' in line:
            date_generated = line.replace('!','')
            break
go_date = date_generated.split('Date Generated by GOC: ')[1]

Path("PathwayMatrix/").mkdir(parents=True, exist_ok=True)
f = 'PathwayMatrix/'+go_date+'-goa_human_F_matrix.txt'
c = 'PathwayMatrix/'+go_date+'-goa_human_C_matrix.txt'
p = 'PathwayMatrix/'+go_date+'-goa_human_P_matrix.txt'
if os.path.exists(f) and os.path.exists(c) and os.path.exists(p):
    print('All matrices exist, done')
    exit()


gene_per_go = {'F':{}, 'C':{}, 'P':{}}
genes_included = {'F':set([]), 'C':set([]), 'P':set([])}
IDs_not_mapped = 0
IDs_mapped = 0
with gzip.open(today+'-goa_human.gaf.gz','rt') as input_file:
    for line in input_file:
        if line.startswith('!'):
            continue
        line = line.strip().split('\t')
        go_type = line[8]
        uniprot_id = line[1]
        go_id = line[4]
        if uniprot_id not in uniprot_to_geneID:
            IDs_not_mapped += 1
            continue
        else:
            IDs_mapped += 1
        for gene_id in  uniprot_to_geneID[uniprot_id]:
            if go_id not in gene_per_go[go_type]:
                gene_per_go[go_type][go_id] = set([])
            gene_per_go[go_type][go_id].add(gene_id)
            genes_included[go_type].add(gene_id)
print('Able to map',IDs_mapped,' uniprot IDs to ensembl')
print('Unable to map',IDs_not_mapped,' uniprot IDs to ensembl')
def write_matrix(go_type, outfile):
    with open(args.ordered_gene_list) as input_file, open(outfile,'w') as out, open(outfile.replace('matrix.txt','genesInPathways.txt'),'w') as out2:
        out.write(date_generated)
        for go_term in gene_per_go[go_type]:
            out.write('\t'+go_term)
        out.write('\n')
        genes_in_at_least_1_pathway = set()
        for gene in input_file.readline():
            gene = gene.strip()
            out.write(gene)
            for go_term in gene_per_go[go_type]:
                if gene in gene_per_go[go_type][go_term]:
                    out.write('\t1.0')
                    genes_in_at_least_1_pathway.add(gene)
                else:
                    out.write('\t0.0')
            out.write('\n')

        for gene in genes_in_at_least_1_pathway:
            out2.write(gene+'\n')
    print('outfile written to',outfile)
write_matrix('F', f)
write_matrix('P', p)
write_matrix('C', c)
