import pandas as pd
import time


print('\n"Protein groups" & "phospho sites" file merger for MaxQuant results')
print('Accepts original txt file from MaxQuant or xlsx & csv too')

protein_file = input('\nEnter protein filepath: (default: proteinGroups.txt)') or 'proteinGroups.txt'
phospho_file = input('Enter phospho filepath: (default: Phospho (STY)Sites.txt)') or 'Phospho (STY)Sites.txt'


if protein_file.endswith('.txt'):
    protein_df = pd.read_table(protein_file, dtype=object)
elif protein_file.endswith('.xlsx'):
    protein_df = pd.read_excel(protein_file)
elif protein_file.endswith('.csv'):
    protein_df = pd.read_csv(protein_file)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')

if phospho_file.endswith('.txt'):
    phospho_df = pd.read_table(phospho_file, dtype=object)
elif phospho_file.endswith('.xlsx'):
    phospho_df = pd.read_excel(phospho_file)
elif phospho_file.endswith('.csv'):
    phospho_df = pd.read_csv(phospho_file)
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
            print(column_prot, 1)
            protein_df = protein_df.drop(column_prot, 1)

        elif column_prot[:20] == 'Identification type ':
            print(column_prot, 1)
            protein_df = protein_df.drop(column_prot, 1)

        elif column_prot[:18] == 'Sequence coverage ':
            print(column_prot, 1)
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
            print(column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:11] == 'Score diff ':
            print(column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:4] == 'PEP ':
            print(column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:6] == 'Score ':
            print(column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:20] == 'Identification type ':
            print(column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:15] == 'Ratio mod/base ':
            print(column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

        elif column_phospho[:10] == 'Occupancy ':
            print(column_phospho, 1)
            phospho_df = phospho_df.drop(column_phospho, 1)

    return protein_df, phospho_df


def remove_REV(protein_df, phospho_df):
    # Remove any reversed peptide
    protein_df = protein_df[protein_df['Protein IDs'].str.contains('REV__') == False]
    phospho_df = phospho_df[phospho_df['Protein'].str.contains('REV__') == False]

    protein_df = protein_df.reset_index(drop=True)
    phospho_df = phospho_df.reset_index(drop=True)

    return protein_df, phospho_df


# Could be improved with list comprehension, I couldn't figure out unpacking list of list w/ the comprehension
# http://stackoverflow.com/questions/3899645/list-extend-and-list-comprehension
def get_original_column_order(phospho_df, merged_df, max_phospho_per_protein):
    # These few blocks of code take the alphanumerically sorted columns
    # and attempts to revert back to the original column order

    def append_str_to_all_list_items(list, string):
        numberified_list = [i + string for i in list]
        return numberified_list

    master_phospho_index_list = []

    # Generates phospho columns based on max number of phosphos per protein
    for i in range(1, max_phospho_per_protein):
        phospho_numberified = append_str_to_all_list_items(list=list(phospho_df),
                                                           string=' (phospho ' + str(i) + ')')
        master_phospho_index_list.extend(phospho_numberified)

    organized_column_index_list = list(protein_df) + master_phospho_index_list

    # Reindex columns by original layout (fixes alphanumeric autosorting bug of pandas.concat)
    # https://github.com/pandas-dev/pandas/issues/4588
    de_alphabetized_df = merged_df.reindex_axis(organized_column_index_list, axis=1)
    return de_alphabetized_df


# This function takes the most time, it appends the phosphopeptides to their respective protein row, line by line
# It would greatly benefit from vectorization but such efforts would be difficult
def merge_prot_phospho(protein_df, phospho_df):
    print('\n\n\nExecuting excel file merge...\n')
    max_phospho_per_protein = 0  # Counts the maximum number of phospho sites detected per protein
                                 # Needed for assembling original template to un-alphabetize the output_df

    for protein_index, row in protein_df.iterrows():
        # Create new temporary protein dataframe for every iteration,
        # it will be used to append the phosphos and then be concatenated to the output dataframe
        protein_temp_df = protein_df.iloc[[protein_index]]

        phospho_counter = 1  # Start phospho name count at zero for each new protein

        # If the phospho fasta matches the protein fasta, append phospho row to protein row
        for phospho_index, row in phospho_df.iterrows():

            if row['Fasta headers'] == protein_temp_df.iloc[0]['Fasta headers']:
                if phospho_counter > max_phospho_per_protein:  # keeps track of how many phosphos per protein
                    max_phospho_per_protein = phospho_counter

                # Create a temporary dataframe with the phospho row to be joined to the temporary protein row
                phospho_temp_df = phospho_df.ix[[phospho_index]]

                phospho_temp_df.index = [protein_index]  # Make row index match the protein's index for merging purposes

                #phospho_df = phospho_df.drop([phospho_index], axis=0)  # Remove phospho row (to detect remaining non-matches at end?)

                # Rename columns of phospho row by phospho counter (phospho0, phospho1 ...)
                column_suffix = ' (phospho ' + str(phospho_counter) + ')'
                for column in phospho_temp_df.columns:
                    phospho_temp_df.rename(columns={column:column+column_suffix},inplace=True)
                phospho_counter += 1

                # Append phospho row to the end of protein row
                protein_temp_df = protein_temp_df.join(phospho_temp_df, how='right')

        if protein_index == 0:
            merged_df = protein_temp_df

        else:
            merged_df = pd.concat([merged_df, protein_temp_df])

        print('Protein row done: ' + str(protein_index+1))


    # "max_phospho..." is the count of the most phosphopeptides for any one protein
    return merged_df, max_phospho_per_protein



# Call functions
protein_df, phospho_df = delete_junk_columns(protein_df, phospho_df)
protein_df, phospho_df = remove_REV(protein_df, phospho_df)
merged_df, max_phospho_per_protein = merge_prot_phospho(protein_df, phospho_df)

output_df = get_original_column_order(phospho_df, merged_df, max_phospho_per_protein)


# Save merged dataframe
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M", current_time)

output_df.to_excel('Merged_files ' + time_string + '.xlsx', sheet_name='Sheet1', index=False)

print('\n\nSaved as: ' 'Merged_files ' + time_string + '.xlsx')


print('\nExcel file organized as such: \n')
print('[ProteinA] - [PhosphoA1] - [PhosphoA2] - [PhosphoA3] ... ')
print('[ProteinB] - [PhosphoB1] - [PhosphoB2] - [PhosphoB3] ... ')
print('[ProteinC] - [PhosphoC1] - [PhosphoC2] - [PhosphoC3] ... ')
print(' ... ')
print(' ... ')


end_secs = time.time()

runsecs = end_secs - start_secs

print('\n Took ' + str(runsecs) + ' seconds')





