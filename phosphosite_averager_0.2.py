# Griffin Calme 2017, Python 3
# This program averages the intensity or ratio of phosphosites with the same gene name.

import pandas as pd
import time

print('Make sure to input table with column headers in first row')
# Import the files, user questions
phospho_file = input('\nEnter txt (tab-delimited), .xlsx or .csv, default: Phospho (STY)Sites.txt)\n') or 'Phospho (STY)Sites.txt'


def file_type_checker(file):
    # Check if the right file type
    if file.endswith('.txt'):
        checked_df = pd.read_table(file, dtype=object)
    elif phospho_file.endswith('.xlsx'):
        checked_df = pd.read_excel(file)
    elif phospho_file.endswith('.csv'):
        checked_df = pd.read_csv(file)
    else:
        raise Exception('Please use tab-delimited (.txt), .xlsx or .csv')

    return checked_df


def dataframe_prep(phospho_df, identifier, ratio, column_list):
    # If any column headers are not string, make it so, otherwise program will fail when appending column names
    phospho_df.columns = [str(x) for x in phospho_df.columns]

    # Make identifier column first column
    identifier_column = phospho_df[identifier]
    phospho_df.drop(labels=[identifier], axis=1, inplace=True)
    phospho_df.insert(0, identifier, identifier_column)

    # Zeros become 1e-8
    phospho_df[ratio] = [0.00000001 if x == 0 else x for x in phospho_df[ratio]]

    original_phospho_df = phospho_df.copy()  # Make a copy for later use

    # Remove phosphosite-specific columns for the averaging section
    phospho_df = phospho_df[column_list]

    return phospho_df, original_phospho_df


def is_up_or_down(phospho_df, ratio):
    phospho_df['up_or_down'] = [0 if x < 1 else 1 if x > 1 else .5 for x in phospho_df[ratio]]

    return phospho_df


def phospho_ratio_averager(phospho_df, identifier, ratio):
    # Averages both the column of interest and the up_or_down column
    print('\nAveraging ...')
    phospho_df_counts = phospho_df.groupby([identifier], axis=0, as_index=False).count()
    phospho_df_counts.rename(columns={ratio: 'counts'},
                             inplace=True)  # In addition to average, get number of phosphosites per protein

    # Overwrites transformed folds with means (cannot simply avg two ratios though, must transform first)
    # .67 becomes -1.5 then -.5, avg(-.5 & .5) = 0 {correct} [will add 1 later in this function],
    # whereas avg(.67 & 1.5) = 1.08 {incorrect}, average of .67 (1.5-fold down) and 1.5 (1.5-fold up) should be 1
    phospho_df['transformed folds'] = [x - 1 if x >= 1 else (-1 / x) + 1 for x in phospho_df[ratio]]

    phospho_df = phospho_df.groupby([identifier], axis=0, as_index=False).mean()

    # Use the average of up_or_down to determine up, down, or mixed etc.
    # Mixed does not mean that phosphorylation average is 0, only that one was up and one was down
    phospho_df['up_or_down'] = ['up' if x == 1 else 'down' if x == 0 else 'mixed equal' if x == 0.5
    else 'mixed decrease' if x < 0.5 else 'mixed increase' for x in phospho_df['up_or_down']]

    # Convert mean "transformed" fold-changes back into mean fold-change
    # (-.5 becomes -1.5 then .67)
    phospho_df['transformed folds'] = [x + 1 if x >= 0 else -1 / (x - 1) for x in phospho_df['transformed folds']]

    # Make it explicit that transformed folds column is now the mean
    phospho_df.rename(columns={'transformed folds': 'MEAN fold change ratio'}, inplace=True)

    # Add in the counts column
    phospho_df['Number of phosphosites'] = phospho_df_counts['counts']

    # Drop incorrect averages
    phospho_df = phospho_df.drop(ratio, 1)

    # Add spacer to the right of the averages section
    phospho_df[' '] = ' '

    return phospho_df


def add_individual_phosphos_to_right(avgd_phospho_df, original_df, identifier):
    max_phosphosites = 0  # Counts the maximum number of phosphosites detected
    # Needed for assembling original template to un-alphabetize the output_df

    # Append original columns to the right
    for index, row in avgd_phospho_df.iterrows():
        avgd_phospho_temp_df = avgd_phospho_df.iloc[[index]]
        phospho_counter = 1  # Start count at zero for each unique identifier

        #print(' ')
        for original_index, original_row in original_df.iterrows():
            if original_row[identifier] == row[identifier]:

                if phospho_counter > max_phosphosites:  # Keeps track of the max num phosphos per protein
                    max_phosphosites = phospho_counter  # Dependency for "get_original_column_order"

                print(original_row[identifier])  # User display

                temp_df_original = original_df.iloc[[original_index]].copy()  # Get line for working temp df

                temp_df_original.index = [index]  # reindex to facilitate correct merging

                # Rename columns by counter (0, 1, 2)
                column_suffix = ' (' + str(phospho_counter) + ')'  # Suffix to distinguish identical cols
                for column in temp_df_original.columns:
                    temp_df_original.rename(columns={column: column + column_suffix}, inplace=True)
                phospho_counter += 1

                # Append original row to end of new df row
                avgd_phospho_temp_df = avgd_phospho_temp_df.join(temp_df_original, how='right')

        if index == 0:
            merged_df = avgd_phospho_temp_df.copy()  # if this is the first row, initialize new df

        else:
            merged_df = pd.concat([merged_df, avgd_phospho_temp_df])  # otherwise merge to right of row

    return merged_df, max_phosphosites


# Could be improved with list comprehension, I couldn't figure out unpacking list of list w/ the comprehension
# http://stackoverflow.com/questions/3899645/list-extend-and-list-comprehension
def get_original_column_order(phospho_df, original_df, merged_df, max_phosphosites):
    # These few blocks of code take the alphanumerically sorted columns
    # and attempts to revert back to the original column order

    def append_str_to_all_list_items(list_in, string):
        numberified_list = [x + string for x in list_in]
        return numberified_list

    master_phospho_index_list = []

    # Generates phospho columns based on max number of phosphos per protein
    for i in range(1, max_phosphosites+1):
        phospho_numberified = append_str_to_all_list_items(list_in=list(original_df),
                                                           string=' (' + str(i) + ')')
        master_phospho_index_list.extend(phospho_numberified)

    organized_column_index_list = list(phospho_df) + master_phospho_index_list

    # Reindex columns by original layout (fixes alphanumeric autosorting bug of pandas.concat)
    # https://github.com/pandas-dev/pandas/issues/4588
    de_alphabetized_df = merged_df.reindex_axis(organized_column_index_list, axis=1)
    return de_alphabetized_df


################################


# call filetype_checker function
my_phospho_df = file_type_checker(phospho_file)

my_ratio = 'blank'
my_identifier = 'blank'

while my_ratio not in my_phospho_df.columns:
    my_ratio = input('Input the column name that you would like to average (default: ratio) ') or 'ratio'

print('\nRatios of zero will be converted to 1e-8 !!! \n')

while my_identifier not in my_phospho_df.columns:
    my_identifier = input('Input the identifier column name (default: Gene Name) ') or 'Gene Name'


my_column_list = [my_identifier, my_ratio]

# Call dataframe prep function
my_phospho_df, my_original_phospho_df = dataframe_prep(my_phospho_df, my_identifier, my_ratio, my_column_list)

# Call function to check if phospho is up or down
my_phospho_df = is_up_or_down(my_phospho_df, my_ratio)

# Call averager function
my_phospho_df = phospho_ratio_averager(my_phospho_df, my_identifier, my_ratio)

# Call functions
my_merged_df, my_max_phosphosites = add_individual_phosphos_to_right(my_phospho_df, my_original_phospho_df, my_identifier)
my_merged_df = get_original_column_order(my_phospho_df, my_original_phospho_df, my_merged_df, my_max_phosphosites)


print('\nDone!')
current_time = time.localtime()
time_string = time.strftime("%m-%d-%y %H_%M_%S", current_time)
output_filename = 'averaged_phosphosites' + time_string + '.xlsx'
my_merged_df.to_excel(output_filename, sheet_name='Sheet1', index=False)

print('\n\nFile saved as ' + output_filename)
