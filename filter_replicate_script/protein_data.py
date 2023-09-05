import pandas as pd
import time
import re

PATH = r"filter_replicate_script/mutation_type_full_mutation.csv"

def create_protein_dict():
    print("Creating protein dictionary...")
    start = time.time()
    protein_df = pd.read_csv(PATH)
    protein_dict = dict()
    pattern = r"\['(.*?)'\]"

    for ind, row in protein_df.iterrows():
        match = re.search(pattern, row['mutation_type'])
        mut_type = ""
        if match: 
            mut_type = match.group(1)
            
        protein_dict[row['full_mutation']] = (mut_type, row['protein'])


    print("Dictionary created")
    print(f"Elapsed time: {(time.time() - start)} sec")
    return protein_dict
    