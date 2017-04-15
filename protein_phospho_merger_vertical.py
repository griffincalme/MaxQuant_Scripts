# Griffin Calme, Python 3
# This script merges proteinGroup.txt and Phospho (STY)Sites.txt (or other similar files)

import pandas as pd
import time


print('\n"Protein groups" & "phospho sites" file merger for MaxQuant results')
print('Accepts original txt file from MaxQuant or xlsx & csv too')

protein_file = input('\nEnter protein filepath: (default: proteinGroups.txt)') or 'proteinGroups.txt'
phospho_file = input('Enter phospho filepath: (default: Phospho (STY)Sites.txt)') or 'Phospho (STY)Sites.txt'


if protein_file.endswith('.txt'):
    my_protein_df = pd.read_table(protein_file, dtype=object)
elif protein_file.endswith('.xlsx'):
    my_protein_df = pd.read_excel(protein_file)
elif protein_file.endswith('.csv'):
    my_protein_df = pd.read_csv(protein_file)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')

if phospho_file.endswith('.txt'):
    my_phospho_df = pd.read_table(phospho_file, dtype=object)
elif phospho_file.endswith('.xlsx'):
    my_phospho_df = pd.read_excel(phospho_file)
elif phospho_file.endswith('.csv'):
    my_phospho_df = pd.read_csv(phospho_file)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')

start_secs = time.time()

print('\nAll reversed "REV" peptides will be removed')
print('The following columns will be dropped: \n')


# Could speed up if I could figure out a vectorized way to drop columns
# Would need to figure out how to create a list of the full column names or regex
# Wouldn't be able to slice column index strings like below
def delete_junk_columns(protein_df, phospho_df):
    # Remove columns from Protein Groups excel sheet that contain:
    # "Peptides" + identifier, "identification type", and "sequence coverage"

    # need to vectorize column drops...
    #prot_junk_list = []

    for column_prot in protein_df.columns:
        if column_prot == 'Peptides':
            pass
        elif column_prot[:7] == 'Razor + ':
            pass

        elif column_prot[:9] == 'Peptides ':
            print('Removed: ' + column_prot, 1)
            protein_df = protein_df.drop(column_prot, 1)

        elif column_prot[:20] == 'Identification type ':
            print('Removed: ' + column_prot, 1)
            protein_df = protein_df.drop(column_prot, 1)

        elif column_prot[:18] == 'Sequence coverage ':
            print('Removed: ' + column_prot, 1)
            protein_df = protein_df.drop(column_prot, 1)

    # Remove columns from Phospho Sites excel sheet that contain:
    # All the extra "Localization prob", "score diff", & "PEP", "identification type",
    # "Ratio mod/base", Occupancy, ratio & error
    for column_phospho in phospho_df.columns:

        if column_phospho == 'Localization prob':
            pass
        elif column_phospho == 'Score diff':
            pass
        elif column_phospho == 'PEP':
            pass
        elif column_phospho == 'Score':
            pass
        elif column_phospho == 'Score for localization':
            pass

        elif column_phospho[:18] == 'Localization prob ':
            print('Removed: ' + column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:11] == 'Score diff ':
            print('Removed: ' + column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:4] == 'PEP ':
            print('Removed: ' + column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:6] == 'Score ':
            print('Removed: ' + column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:20] == 'Identification type ':
            print('Removed: ' + column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:15] == 'Ratio mod/base ':
            print('Removed: ' + column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:10] == 'Occupancy ':
            print('Removed: ' + column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

    return protein_df, phospho_df


# Remove any reversed peptide
def remove_REV(protein_df, phospho_df):
    protein_df = protein_df[protein_df['Protein IDs'].str.contains('REV__') == False]
    phospho_df = phospho_df[phospho_df['Protein'].str.contains('REV__') == False]

    protein_df = protein_df.reset_index(drop=True)
    phospho_df = phospho_df.reset_index(drop=True)

    return protein_df, phospho_df


# unpack multiple fasta headers into their own rows
def unpack_IDs(df, header_string, delimiter_string):
    s = df[header_string].str.split(delimiter_string, expand=True).stack()
    i = s.index.get_level_values(0)
    unpacked_fasta_df = df.loc[i].copy()
    unpacked_fasta_df[header_string] = s.values

    return unpacked_fasta_df


# This function flips merged_df from [phos][prot] to [prot][phos]
def swap_phospho_and_protein_sections(merged_df,start_of_protein_string='Protein IDs',
                                      start_of_phospho_string='Proteins'):

    # Get the start location of the protein section
    protIDs_index_prot = merged_df.columns.get_loc(start_of_protein_string)

    # Get the start location of the phospho section
    prot_index_phos = merged_df.columns.get_loc(start_of_phospho_string)

    merged_df_prot_slice = merged_df.iloc[:, protIDs_index_prot:]  # slice from start of prot section to end

    # slice from beginning to end of phos section
    merged_df_phos_slice = merged_df.iloc[:,prot_index_phos:protIDs_index_prot]

    # concat the dfs back together by index
    output_df = pd.concat([merged_df_prot_slice, merged_df_phos_slice], axis=1)

    return output_df


# convert integer index to excel index style ([27] = [AA])
# for assisting the user in finding the start and end of each section
def col_to_excel(col): # col is 0 based
    col += 1
    excel_col = str()
    div = col
    while div:
        (div, mod) = divmod(div-1, 26) # will return (x, 0 .. 25)
        excel_col = chr(mod + 65) + excel_col

    return excel_col


##############


# Call functions
my_protein_df, my_phospho_df = delete_junk_columns(my_protein_df, my_phospho_df)
my_protein_df, my_phospho_df = remove_REV(my_protein_df, my_phospho_df)

protein_unpacked_fasta_df = unpack_IDs(my_protein_df, 'Fasta headers', ';')

# remove extra headers from phospho "Fasta headers" column
my_phospho_df['Fasta headers'] = my_phospho_df['Fasta headers'].apply(lambda x: x.split(';')[0])


# merge phospho file with matching entries from expanded protein file
print('\n\n\nExecuting excel file merge...\n')
print('This could take a minute...')
my_merged_df = my_phospho_df.merge(protein_unpacked_fasta_df, on='Fasta headers', how='inner', suffixes=('_phos', '_prot'))

# swap phospho and protein sections (so that protein is to the right and phospho is to the left
my_output_df = swap_phospho_and_protein_sections(my_merged_df, 'Protein IDs', 'Proteins')


# Save merged dataframe
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M", current_time)

my_output_df.to_excel('Merged_files ' + time_string + '.xlsx', sheet_name='Sheet1', index=False)


################


print('\nExcel file organized as such: \n')
print('[ProteinA] - [PhosphoA1]')
print('[ProteinA] - [PhosphoA2]')
print(' ... ')
print('[ProteinB] - [PhosphoB1]')
print('[ProteinB] - [PhosphoB2]')
print(' ... ')
print('[ProteinC] - [PhosphoC1]')
print('[ProteinC] - [PhosphoC2]')
print(' ... ')


# Get the column locations for the start of each section, for user interface below
protein_location = my_output_df.columns.get_loc('Protein IDs')
protein_location_excel = col_to_excel(protein_location)

phospho_location = my_output_df.columns.get_loc('Proteins')
phospho_location_excel = col_to_excel(phospho_location)

print('\nProtein section begins at Protein IDs, column=' + protein_location_excel)
print('Phospho section begins at Proteins, column=' + phospho_location_excel)

end_secs = time.time()

runsecs = end_secs - start_secs

print('\nSaved as: ' 'Merged_files ' + time_string + '.xlsx')

print('\nTook ' + str(runsecs) + ' seconds')






