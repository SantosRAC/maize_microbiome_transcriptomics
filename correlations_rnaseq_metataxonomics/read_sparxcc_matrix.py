#!/usr/bin/env python3

# Read file line by line

with open('SparCC_output_day_genus.txt', 'r') as fin:

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
                print(otu_name, gene_names[i], correlations_matrix[i])

