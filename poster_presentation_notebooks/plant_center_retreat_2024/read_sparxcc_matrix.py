#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='SparXCC output matrix (\$cor object from R list)')
parser.add_argument('--input', type=str, metavar='SparXCC_output.txt',
                    dest='sparxcc_matrix_file',
                    help='The TSV file with the cross-correlation matrix output from SparXCC',
                    required=True)
parser.add_argument('--output', type=str, metavar='output.txt',
                    dest='output_file',
                    help='Output file with the significant correlations (pairs of OTUs and genes passing the permutation threshold)',
                    required=True)

args = parser.parse_args()

sparxcc_matrix = args.sparxcc_matrix_file
output_file = args.output_file

output_file_obj = open(output_file, 'w')

with open(sparxcc_matrix, 'r') as fin:

    colnames = fin.readline().rstrip().split()
    halfmatrix = int(len(colnames) / 2)
    gene_names = [n.replace('cor.','').replace('"', '') for n in colnames[:halfmatrix]]

    for line in fin:
        otu_name, *line_fields = line.rstrip().split()

        # Retrieve the correlation matrix, the m value and the boolean matrix
        correlations_matrix = line_fields[:halfmatrix]
        m_value = line_fields[halfmatrix]
        m_boolean_matrix = line_fields[halfmatrix + 1:]

        # Get significant correlations
        for i, val in enumerate(m_boolean_matrix):
            if val == 'TRUE':
                output_file_obj.write(otu_name, gene_names[i], correlations_matrix[i], "\n")

output_file_obj.close()