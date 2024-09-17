#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Updates the sample groups in the metadata file")
parser.add_argument('--input', type=str, metavar='sample_metadata_conekt_microbiome.txt',
                    dest='sample_metadata_file',
                    help='The TSV file with sample metadata used in CoNekT Grasses Microbiome',
                    required=True)
parser.add_argument('--association_file', type=str, metavar='association_file.txt',
                    dest='association_file',
                    help='The file with associations between sample identifiers and groups.',
                    required=True)
parser.add_argument('--output', type=str, metavar='sample_metadata_updated.txt',
                    dest='output_file',
                    help='Output file with the updated sample metadata',
                    required=True)

args = parser.parse_args()

sample_metadata_file = args.sample_metadata_file
output_file = args.output_file
association_file = args.association_file

input_group_type = None

with open(association_file, 'r') as fin:

    _, input_group_type = fin.readline().rstrip().split()
    
    new_sample_groups_dict = {}
    
    for line in fin:
        sample_id, sample_group = line.rstrip().split("\t")
        new_sample_groups_dict[sample_id] = sample_group

output_file_obj = open(output_file, 'w')

with open(sample_metadata_file, 'r') as fin:

    colnames = fin.readline().rstrip().split()

    if (colnames[0] == 'SampleID') and\
        (colnames[1] == 'DOI') and\
        (colnames[2] == 'ConditionDescription') and\
        (colnames[3] == 'Replicate') and\
        (colnames[4] == 'PO_anatomy') and\
        (colnames[5] == 'PO_dev_stage') and\
        (colnames[6] == 'PECO') and\
        (colnames[7] == 'ENVO') and\
        (colnames[8] == 'Genotype') and\
        (colnames[9] == 'Sample_groups'):
        output_file_obj.write(f"{colnames[0]}\t{colnames[1]}\t{colnames[2]}\t{colnames[3]}\t{colnames[4]}\t{colnames[5]}\t{colnames[6]}\t{colnames[7]}\t{colnames[8]}\t{colnames[9]}\n")
    else:
        print("The metadata file does not have the expected column names. Exiting...")
        exit(1)

    previous_sample_groups_dict = {}

    for line in fin:
        sample_id, doi, condition_desc, replicate, po_anatomy, po_dev_stage, peco, envo, genotype, sample_groups = line.rstrip().split("\t")
        sample_groups = sample_groups.split(";")

        if sample_id in previous_sample_groups_dict.keys():
            print(f"Sample {sample_id} is duplicated in the metadata file. Exiting...")
            exit(1)
        
        previous_sample_groups_dict[sample_id] = {}

        for group in sample_groups:
            group_type, group_value = group.split(":")

            if input_group_type == group_type:
                if sample_id in new_sample_groups_dict.keys():
                    previous_sample_groups_dict[sample_id][group_type] = new_sample_groups_dict[sample_id]
                else:
                    print(f"Sample {sample_id} does not have an association in the association file. Exiting...")
                    exit(1)
            else:
                previous_sample_groups_dict[sample_id][group_type] = group_value
        
        output_file_obj.write(f"{sample_id}\t{doi}\t{condition_desc}\t{replicate}\t{po_anatomy}\t{po_dev_stage}\t{peco}\t{envo}\t{genotype}"+"\t"+';'.join([f'{key}: {value}' for key, value in previous_sample_groups_dict[sample_id].items()])+"\n")

output_file_obj.close()