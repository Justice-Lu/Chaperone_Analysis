import pandas as pd 
import numpy as np 

def normalize_metric(df, normalize_by, normalize_metric, groupby_cols=['exp_date', 'OR']):
    # Group by 'normalize_by' and 'exp_date' to calculate the mean 'normalize_metric' for the specified group
    normalize_mean = df[df['Sample_name'] == normalize_by].groupby(groupby_cols)[normalize_metric].mean().reset_index()
    normalize_mean.columns = groupby_cols +  ['normalize_mean']
    
    # Merge the mean 'normalize_metric' back into the original dataframe based on 'exp_date'
    df = pd.merge(df, normalize_mean, on=groupby_cols, how='left')
    
    # Divide each row's 'normalize_metric' by the mean value calculated for 'normalize_by'
    normalized_column_name = f'Normalized_{normalize_metric}'
    df[normalized_column_name] = df[normalize_metric] / df['normalize_mean']
    
    # Drop the additional column
    df.drop(columns=['normalize_mean'], inplace=True)
    
    return df

# def filter_expdate_by_samples(df, contain_samples):
#     # Initialize a list to store exp_dates that meet the condition
#     valid_exp_dates = []
    
#     # Iterate through each exp_date
#     for _exp_date in df['exp_date'].unique():
#         # Filter the dataframe for the current exp_date
#         exp_date_df = df[df['exp_date'] == _exp_date]

#         # Check if all contain_samples are present in Sample_name
#         if all([_string in exp_date_df['Sample_name'].values for _string in contain_samples]):
#         # if all([_sample_name for _sample_name in exp_date_df['Sample_name'].values for _string in contain_samples if _sample_name]):
#             # If all are present, add the exp_date to the list of valid exp_dates
#             valid_exp_dates.append(_exp_date)

#     # Filter the original dataframe based on the valid exp_dates
#     filtered_df = df[df['exp_date'].isin(valid_exp_dates)]
    
#     return filtered_df

def filter_expdate_by_samples(df, contain_samples, match_mode='exact'):
    """
    Filters the DataFrame to only include rows from exp_dates that contain all specified samples,
    either by exact match or substring presence.

    Parameters:
    - df: pandas.DataFrame with 'exp_date' and 'Sample_name' columns.
    - contain_samples: list of sample names or substrings to look for.
    - match_mode: 'exact' (default) or 'contains'. 
        'exact' requires Sample_name values to exactly match.
        'contains' requires the substrings to appear somewhere in the Sample_name values.

    Returns:
    - Filtered DataFrame with only rows from exp_dates that meet the condition.
    """
    valid_exp_dates = []

    for _exp_date in df['exp_date'].unique():
        exp_date_df = df[df['exp_date'] == _exp_date]
        sample_names = exp_date_df['Sample_name'].astype(str)

        if match_mode == 'exact':
            if all(_string in sample_names.values for _string in contain_samples):
                valid_exp_dates.append(_exp_date)
        elif match_mode == 'contains':
            if all(sample_names.str.contains(_string, regex=False).any() for _string in contain_samples):
                valid_exp_dates.append(_exp_date)
        else:
            raise ValueError("match_mode must be either 'exact' or 'contains'")

    return df[df['exp_date'].isin(valid_exp_dates)]