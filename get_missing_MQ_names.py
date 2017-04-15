#Griffin Calme 2017, Python 3

import pandas as pd
import time

print('\n"Fixes blank "Protein names" and "Gene names" entries from MaxQuant output files')
print('Accepts original txt file from MaxQuant or xlsx & csv too')

MQ_file = input('\nEnter protein or PTM (e.g. Phospho (STY)Sites.txt filepath: (default: proteinGroups.txt)\n') or 'proteinGroups.txt'


if MQ_file.endswith('.txt'):
    MQ_df = pd.read_table(MQ_file, dtype=object)
elif MQ_file.endswith('.xlsx'):
    MQ_df = pd.read_excel(MQ_file)
elif MQ_file.endswith('.csv'):
    MQ_df = pd.read_csv(MQ_file)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')


for MQ_index, row in MQ_df.iterrows():

    # Skip reversed and contaminants
    if pd.isnull(row[0]) or row[0].startswith('REV') or row[0].startswith('CON'):
        pass

    else:
        # If 'Protein names' is blank, then get name from fasta
        if pd.isnull(row['Protein names']):
            fasta_header = str(row['Fasta headers'])  # Get fasta header

            protein_name = fasta_header.split(" ", 1)
            protein_name = protein_name[1]  # Slice beginning of protein name
            protein_name = protein_name.split(" OS=", 1)
            protein_name = protein_name[0]  # Slice end of protein name

            row['Protein names'] = protein_name
            print('Found protein name: ' + protein_name)

        # If 'Gene names' is blank, then get name from fasta
        if pd.isnull(row['Gene names']):
            fasta_header = str(row['Fasta headers'])  # Get fasta header
            if "GN=" in fasta_header:
                gene_name = fasta_header.split("GN=", 1)
                gene_name = gene_name[1]  # Slice beginning of gene name
                gene_name = gene_name.split(" ", 1)
                gene_name = gene_name[0]  # Slice end of gene name
                gene_name = gene_name.split(";", 1) # sometimes GN is end of header, followed by semicolon
                gene_name = gene_name[0]

                row['Gene names'] = gene_name
                print('Found gene name: ' + gene_name)

            else:
                pass


current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M_%S", current_time)
output_filename = 'Fixed_MQ_output' + time_string + '.xlsx'
MQ_df.to_excel(output_filename, sheet_name='Sheet1', index=False)

print('\n\nFile saved as ' + output_filename)
