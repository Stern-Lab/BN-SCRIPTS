from datetime import datetime
import pandas as pd
import numpy as np
import Filter_Usecase.filter_replicate_script.list_PROBLEMATIC_positions as prob
import os

PROBLEMATIC = prob.list_PROBLEMATIC_positions()

def _assign_transition_type(mutation_type):
    """
    Provides a transition state associated with a given mutation type.
    :param mutation_type: The mutation type for which an associated transition state will be provided.
    :return: A transition state associated with the given mutation type.
    """
    transitions = ['AG', 'GA', 'TC', 'CT']
    oxidation = ['CA', 'GT']
    transversions = ['AC', 'TG', 'TA', 'AT', 'GC', 'CG']

    if not isinstance(mutation_type, str):
        return 'err'
    if mutation_type in transitions:
        return 'ts'
    elif mutation_type in oxidation:
        return 'ox'
    elif mutation_type[0] == '-':
        return 'ins'
    elif mutation_type[1] == '-':
        return 'del'
    elif mutation_type in transversions:
        return 'tv'
    else:
        return 'ref'

def enrich_mutation(df):
    """
    Assigns basic information on mutation using three columns of the given df: ref_base (the reference of this
    position), read_base (what was actually read by the machine), ref_pos (the position in the genome).
    For example mutation A1234G is a transition of A (ref_base) in position 1234 (ref_pos) to a G (read_base).
    :param df: A Dataframe to which 'transition', 'type', and 'mutation' columns will be added.
    :return: A dataframe identical to the received one but with the aforementioned columns added.
    """
    df['transition'] = df['ref_base'] + df['read_base']
    df['type'] = df.transition.map(_assign_transition_type)
    df['mutation'] = df['ref_base'] + df['ref_pos'].astype(int).map(str) + df['read_base']
    return df

def filter_ref(df):
    """
    Filters a dataframe from non-mutation positions(positions where the ref_base is the same as read_base).
    :param df: A Dataframe from which ref mutations will be filtered out.
    :return: A new Dataframe identical to the received one, but without reference mutations.
    """
    return df.loc[df['type'] != 'ref']

def filter_non_mutations(df):
    ret_df = df.loc[(df['ref_base'] != df['read_base'])]
    return ret_df[~ret_df['ref_pos'].astype(int).isin(PROBLEMATIC)]

def new_freq_insert(rep_df, coverage, freq, base_count):
    for ind, row in rep_df.iterrows():
        if row["coverage"] <= coverage:
            if row["frequency"] <= freq:
                rep_df.loc[ind, "new_freq"] = "X"

            elif row["frequency"] > freq:
                rep_df.loc[ind, "new_freq"] = -1
        
        elif row["coverage"] > coverage:
            if row["frequency"] <= freq:
                rep_df.loc[ind, "new_freq"] = 0

            elif row["frequency"] > freq:
                if row["base_count"] < base_count:
                    rep_df.loc[ind, "new_freq"] = -1
                elif row["base_count"] > base_count:
                    rep_df.loc[ind, "new_freq"] = row["frequency"]
    return rep_df

def final_freq_calc(freq1, freq2, coverage, freq, base_count):
    


def filter(tsv1, tsv2, freq, coverage, base_count, protein_dict, result_dir):
    # Read file & add columns & filter mutations to each pair of replica's freqs file
    
    # All mutation WITHOUT FILTERING
    rep1_df_all = pd.read_csv(tsv1, sep='\t')
    rep1_df_all = filter_non_mutations(rep1_df_all)
    rep1_df_all = enrich_mutation(rep1_df_all)
    rep1_df_all = filter_ref(rep1_df_all)

    rep2_df_all = pd.read_csv(tsv2, sep='\t')
    rep2_df_all = filter_non_mutations(rep2_df_all)
    rep2_df_all = enrich_mutation(rep2_df_all)
    rep2_df_all = filter_ref(rep2_df_all)
    
    num_of_mut_rep1 = rep1_df_all.shape[0]
    num_of_mut_rep2 = rep2_df_all.shape[0]
    
    # Save all mutation except ref and PROBLEMATIC
    rep1_df = rep1_df_all[["ref_pos", "mutation", "base_count", "coverage", "frequency"]].copy()
    rep1_df.to_csv(f"{result_dir}/freq1_all_mutations.csv", index=False)
    rep2_df = rep2_df_all[["ref_pos", "mutation", "base_count", "coverage", "frequency"]].copy()
    rep2_df.to_csv(f"{result_dir}/freq2_all_mutations.csv", index=False)

    # Insert new frequency to each mutation according to each file independently
    rep1_df["new_freq"] = "new"
    rep2_df["new_freq"] = "new"

    new_rep1_df = new_freq_insert(rep1_df, coverage, freq, base_count)
    new_rep2_df = new_freq_insert(rep2_df, coverage, freq, base_count)
    
    new_rep1_df.to_csv(f"{result_dir}/freq1_independent.csv", index=False)
    new_rep2_df.to_csv(f"{result_dir}/freq2_independent.csv", index=False)
    
    # Insert new frequency to each mutation according to both freqs


    # Create DF for next phase (BN algorithem)
    filtered_rep1 = rep1_df[["mutation", "frequency"]]
    filtered_rep1.to_csv(f"{result_dir}/mut_freq_3.csv", index=False)
    
    filtered_rep2 = rep2_df[["mutation", "frequency"]]
    filtered_rep2.to_csv(f"{result_dir}/mut_freq_4.csv", index=False)
    indepandent = 
    # Merge data framse based on common mutation after filtering
    merged_df = pd.merge(rep1_df, rep2_df, how="inner", on="mutation")
    merged_all_df = pd.merge(rep1_df_all, rep2_df_all, how="inner", on="mutation")
    num_of_mut_merged = merged_df.shape[0]
    num_of_mut_merged_all = merged_all_df.shape[0]
    
    # Add protein (mutation type) column
    merged_df["mutation_type"] = ""

    # Usecase submitting
    for ind, row in merged_df.iterrows():
        rep1_freq = row['frequency_x']
        rep2_freq = row['frequency_y']
        if row['mutation'] in protein_dict:
            merged_df.loc[ind, 'mutation_type'] = protein_dict[row['mutation']][0]

        big_freq = max(rep1_freq, rep2_freq)
        small_freq = min(rep1_freq, rep2_freq)
    
    return usecase_df, 1, 1, 1, 1