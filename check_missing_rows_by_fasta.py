# Python 3 Griffin Calme
# This program checks if a result file has dropped any rows from an original file
# It checks for missing rows by fasta header

import pandas as pd
import time


#original = input('\nEnter protein filepath: (default: proteinGroups.txt)') or 'proteinGroups.txt'
original = 'Phospho (STY)Sites.txt'
new = 'Merged_files 04-14-17 23_34.xlsx'
#new = input('Enter phospho filepath: (default: Phospho (STY)Sites.txt)') or 'Phospho (STY)Sites.txt'


if original.endswith('.txt'):
    original_df = pd.read_table(original, dtype=object)
elif original.endswith('.xlsx'):
    original_df = pd.read_excel(original)
elif original.endswith('.csv'):
    original_df = pd.read_csv(original)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')

if new.endswith('.txt'):
    new_df = pd.read_table(new, dtype=object)
elif new.endswith('.xlsx'):
    new_df = pd.read_excel(new)
elif new.endswith('.csv'):
    new_df = pd.read_csv(new)
else:
    raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')

start_secs = time.time()

def remove_REV(phospho_df):
    # Remove any reversed peptide
    phospho_df = phospho_df[phospho_df['Protein'].str.contains('REV__') == False]

    phospho_df = phospho_df.reset_index(drop=True)

    return phospho_df

original_df = remove_REV(original_df)
new_df = remove_REV(new_df)



original_df['Fasta headers'] = original_df['Fasta headers'].apply(lambda x: x.split(';')[0])


original_df = original_df['Fasta headers'].tolist()
new_df = new_df['Fasta headers'].tolist()

list = []
for i in new_df:
    if i not in original_df:
        list.append(i)
        print(i)

#df = pd.merge(new_df, original_df, how='outer', indicator=True)
#rows_in_df1_not_in_df2 = df[df['_merge']=='left_only'][new_df.columns]



# Save merged dataframe
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M", current_time)

#rows_in_df1_not_in_df2.to_excel('rows not in new file ' + time_string + '.xlsx', sheet_name='Sheet1', index=False)

f=open('missing_from_merge.txt','w')
s1='\n'.join(list)
f.write(s1)
f.close()