#!/opt/conda/envs/eQTLGenPopAssign/bin/python

import sys

# Read population assignment path input
population_assignment_file = sys.argv[1]

# Open population assignment file and read lines
with open(population_assignment_file, 'r') as pop_file:
    lines = [line.strip().split('\t') for line in pop_file.readlines()][1:]

# Create population translation dict
pop_dict = {
    0: 'EUR',
    1: 'EAS',
    2: 'AMR',
    3: 'SAS',
    4: 'AFR',
}

# Initialize res dict
res_dict = {
    'EUR': [],
    'EAS': [],
    'AMR': [],
    'SAS': [],
    'AFR': [],
}

# For every line in population assigment file, get lowest distance and assign population
for line in lines:
    sample = line[0]
    distances = [float(distance) for distance in line[1:]]
    population = pop_dict[distances.index(min(distances))]
    res_dict[population].append(sample)

# # Write summary to count file
with open('population_counts.txt', 'w') as count_file:
    for k, v in res_dict.items():
        count_file.write(f'{k}: {len(v)} samples\n')

with open('population_file.txt', 'w') as pop_file:
    for pop, samples in res_dict.items():
        for s in samples:
            pop_file.write(f'{s},{pop}\n')

# For every population, write out samples to file if n > 49
for k, v in res_dict.items():
    if len(v) > 1:
        with open(f'populations/{k}', 'w') as pop_out_file:
            for sample in v:
                pop_out_file.write(f'{sample}\n')

