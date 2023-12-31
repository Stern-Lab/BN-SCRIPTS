from datetime import datetime
import pandas as pd
import os

def join_time_samples(f_name1, f_name2):
    """
    Creating input format required for BN script
    :param:
    :return: A txt file formatted as an input to the BN script
    """
    df1 = pd.read_csv(f"./{f_name1}.csv")
    df2 = pd.read_csv(f"./{f_name2}.csv")

    merged_df = pd.merge(df1, df2, how='outer', on="mutation")

    # Fill NaN with zeros
    merged_df = merged_df.fillna(0)
    
    txt = ""
    for ind, row in merged_df.iterrows():
        txt = txt + "{:.6f}".format(row['frequency_x']) + "\t" + "{:.6f}".format(row['frequency_y']) + "\n"

    with open("./results.txt", 'w') as file:
        file.write(txt)

    return