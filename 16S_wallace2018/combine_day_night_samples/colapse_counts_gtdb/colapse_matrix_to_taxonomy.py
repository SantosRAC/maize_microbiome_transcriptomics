#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Collapse counts matrix to genus level')
parser.add_argument('--verbose', action="store_true", dest='verbose')
parser.add_argument('--gtdb_classification', type=str, metavar='gtdb_classification.tsv',
                    dest='gtdb_classification_file',
                    help='GTDB classification file using Qiime2 classifier', required=True)
parser.add_argument('--counts_matrix', type=str, metavar='otu_counts.tsv',
                    dest='counts_file',
                    help='OTU Counts table', required=True)
parser.add_argument('--genus_counts_matrix', type=str, metavar='genus_counts.tsv',
                    dest='genus_counts_file',
                    help='OTU Counts table (Genus)', required=False)
parser.add_argument('--family_counts_matrix', type=str, metavar='family_counts.tsv',
                    dest='family_counts_file',
                    help='OTU Counts table (Family)', required=False)

args = parser.parse_args()

classification_file = args.gtdb_classification_file
counts_file = args.counts_file

family_count_matrix = args.family_counts_file
genus_count_matrix = args.genus_counts_file

if not family_count_matrix and not genus_count_matrix:
    print('\nERROR: At least one of the output files must be specified\n')
    parser.print_help()
    exit(1)

record2genus = {}
record2family = {}

total_otus = 0
no_genus_no_family = 0
no_genus = 0

total_otu_counts_genus = 0
discarded_otu_counts_no_genus = 0

total_otu_counts_family = 0
discarded_otu_counts_no_family = 0

with open(classification_file, 'r') as file:
    for line in file:
        line = line.strip()
        feature_name, taxon_path, confidence_value = line.split('\t')
        taxon_path = taxon_path.split(';')
        for taxon in taxon_path:
            if taxon.startswith('g__'):
                genus = taxon.split('__')[1]
                if feature_name in record2genus.keys():
                    print(f'Something very weird is happening here.')
                else:
                    record2genus[feature_name] = genus
        for taxon in taxon_path:
            if taxon.startswith('f__'):
                family = taxon.split('__')[1]
                if feature_name in record2family.keys():
                    print(f'Something very weird is happening here.')
                else:
                    record2family[feature_name] = family

        if feature_name not in record2genus.keys() and\
            feature_name not in record2family.keys() and args.verbose:
            no_genus_no_family += 1
            print(f'WARNING: No genus/family found for {feature_name}: {taxon_path}')
        
        if feature_name not in record2genus.keys() and\
            feature_name in record2family.keys() and args.verbose:
            no_genus += 1
            print(f'WARNING: No genus found for {feature_name}: {taxon_path}')

genus2count = {}
family2count = {}
sample_names = []
first_column_name = ''

with open(counts_file, 'r') as file:
    header = file.readline().strip()
    sample_names = header.split('\t')
    first_column_name = sample_names.pop(0)
    for line in file:
        line = line.strip()
        sample_counts = line.split('\t')
        otu_name = sample_counts.pop(0)
        total_otus += 1
        for (sample_name, sample_count) in zip(sample_names, sample_counts):
            if otu_name in record2genus.keys():
                genus = record2genus[otu_name]
                if genus in genus2count.keys():
                    if sample_name in genus2count[genus].keys():
                        genus2count[genus][sample_name] += int(float(sample_count))
                        total_otu_counts_genus += int(float(sample_count))
                    else:
                        genus2count[genus][sample_name] = int(float(sample_count))
                        total_otu_counts_genus += int(float(sample_count))
                else:
                    genus2count[genus] = {}
                    genus2count[genus][sample_name] = int(float(sample_count))
                    total_otu_counts_genus += int(float(sample_count))
            else:
                discarded_otu_counts_no_genus += int(float(sample_count))
                total_otu_counts_genus += int(float(sample_count))
            
            if otu_name in record2family.keys():
                family = record2family[otu_name]
                if family in family2count.keys():
                    if sample_name in family2count[family].keys():
                        family2count[family][sample_name] += int(float(sample_count))
                        total_otu_counts_family += int(float(sample_count))
                    else:
                        family2count[family][sample_name] = int(float(sample_count))
                        total_otu_counts_family += int(float(sample_count))
                else:
                    family2count[family] = {}
                    family2count[family][sample_name] = int(float(sample_count))
                    total_otu_counts_family += int(float(sample_count))
            else:
                discarded_otu_counts_no_family += int(float(sample_count))
                total_otu_counts_family += int(float(sample_count))

if genus_count_matrix:
    with open(genus_count_matrix, 'w') as outfile:
        outfile.write(f'{first_column_name}')
        for sample_name in sample_names:
            outfile.write(f'\t{sample_name}')
        outfile.write('\n')
        for genus in genus2count.keys():
            outfile.write(f'{genus}')
            for sample_name in sample_names:
                outfile.write(f'\t{genus2count[genus][sample_name]}')
            outfile.write('\n')

if family_count_matrix:
    with open(family_count_matrix, 'w') as outfile:
        outfile.write(f'{first_column_name}')
        for sample_name in sample_names:
            outfile.write(f'\t{sample_name}')
        outfile.write('\n')
        for family in family2count.keys():
            outfile.write(f'{family}')
            for sample_name in sample_names:
                outfile.write(f'\t{family2count[family][sample_name]}')
            outfile.write('\n')

if args.verbose:
    print(f'\n\nColapsing SUMMARY:::....*****')
    print(f'Total OTUs: {total_otus}')
    print(f'Number of records without genus/family: {no_genus_no_family}')
    print(f'Number of records without genus: {no_genus}')
    print(f'Number of genera: {len(genus2count.keys())}')
    print(f'Number of families: {len(family2count.keys())}')
    if genus_count_matrix:
        print(f'Total Number of Counts: {total_otu_counts_genus}')
        print(f'Total Counts Discarded (no genus): {discarded_otu_counts_no_genus} ({round((discarded_otu_counts_no_genus / total_otu_counts_genus) * 100, 2)}%)')
        print(f'Created Genus count matrix: {genus_count_matrix}')
    if family_count_matrix:
        print(f'Total Number of Counts: {total_otu_counts_family}')
        print(f'Total Counts Discarded (no family): {discarded_otu_counts_no_family} ({round((discarded_otu_counts_no_family / total_otu_counts_family) * 100, 2)}%)')
        print(f'Created Family count matrix: {family_count_matrix}')