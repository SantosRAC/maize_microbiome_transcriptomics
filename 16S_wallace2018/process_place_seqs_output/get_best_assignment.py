#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Process place_seqs output and get best assignment and add to feature table')
parser.add_argument('-per_query_file', help='per_query.tsv', dest='per_query_tsv_file', required=True,
                    type=str, metavar='per_query.tsv file (output from Gappa)')
parser.add_argument('-feature_table', help='feature_table.tsv', dest='feature_table_file', required=True,
                    type=str, metavar='feature_table.tsv (output from Qiime2 feature table (.biom) converted to .tsv format using biom)')
args = parser.parse_args()

per_query_tsv_file = args.per_query_tsv_file
feature_table_file = args.feature_table_file

asv2taxon = {}

# Read per_query.tsv
with open(per_query_tsv_file, 'r') as file:

    _ = file.readline()

    for line in file:
        # Split the line into columns
        columns = line.strip().split('\t')
        
        # Assign columns to variables
        asv_id = columns[0]
        lwr_value = float(columns[1])
        taxopath = columns[5]
        
        if asv_id in asv2taxon.keys():
            if lwr_value in asv2taxon[asv_id].keys():
                asv2taxon[asv_id][lwr_value].append(taxopath)
            else:
                asv2taxon[asv_id][lwr_value] = [taxopath]
        else:
            asv2taxon[asv_id] = {lwr_value: [taxopath]}

asv2definitive_taxon = {}

# Print asv_ids
for asv_id in asv2taxon.keys():
    highest_lwr_value = 0
    for lwr_value in asv2taxon[asv_id].keys():
        if lwr_value > highest_lwr_value:
            highest_lwr_value = lwr_value
            asv_taxa = asv2taxon[asv_id][lwr_value]
        elif lwr_value == highest_lwr_value:
            asv_taxa = asv2taxon[asv_id][lwr_value]
        else:
            pass

    asv2definitive_taxon[asv_id] = asv_taxa[len(asv_taxa) - 1]


# Read feature_table.tsv
with open(feature_table_file, 'r') as file:

    for line in file:

        if line.startswith('ASV\t'):
            columns = line.strip().split('\t')
            asv_id = columns.pop(0)
            print(f'{asv_id}\tTaxon\t'+"\t".join(columns))
            continue

        # Split the line into columns
        columns = line.strip().split('\t')
        
        # Assign columns to variables
        asv_id = columns.pop(0)
        if asv_id not in asv2definitive_taxon.keys():
            print(f'{asv_id}\tunplaced\t'+"\t".join(columns))
        else:
            print(f'{asv_id}\t{asv2definitive_taxon[asv_id]}\t'+"\t".join(columns))

