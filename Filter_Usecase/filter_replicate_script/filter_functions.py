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

def filter_deletion_insertion(df):
    ret_df = df.loc[df['ref_base'] != "-"]
    ret_df = ret_df.loc[df['read_base'] != "-"]
    return ret_df

def filter_non_mutations(df):
    ret_df = df.loc[(df['ref_base'] != df['read_base'])]
    return ret_df[~ret_df['ref_pos'].astype(int).isin(PROBLEMATIC)]

def calc_weighted_avg(bs1, bs2, cvg1, cvg2):
    return (bs1 + bs2) / (cvg1 + cvg2)

def phase0_filter(tsv, filter_indels):
    rep_df_all = pd.read_csv(tsv, sep='\t')
    if filter_indels:
        rep1_df_all = filter_deletion_insertion(rep_df_all)
    rep_df_all = filter_non_mutations(rep_df_all)
    rep_df_all = enrich_mutation(rep_df_all)
    rep_df_all = filter_ref(rep_df_all)

    return rep_df_all

def phase1_filter(rep_df, coverage, freq, base_count):
    for ind, row in rep_df.iterrows():
        if row["coverage"] < coverage:
            rep_df.loc[ind, "new_freq"] = -1
        
        else:
            if row["frequency"] < freq:
                rep_df.loc[ind, "new_freq"] = 0

            else:
                if row["base_count"] < base_count:
                    rep_df.loc[ind, "new_freq"] = -1
                else:
                    rep_df.loc[ind, "new_freq"] = row["frequency"]
    return rep_df

def phase2_filter(freq1, freq2, protein_dict):
    # Merge df inner method
    merged_df = pd.merge(freq1, freq2, on="mutation", how="inner")
    
    # Add protein (mutation type) column
    merged_df["mutation_type"] = ""
    
    # Add final frequency column
    merged_df["final_freq"] = "new"
    merged_df["tot_cov"] = -1
    merged_df["tot_base_count"] = -1

    for ind, row in merged_df.iterrows():
        if row["mutation"] in protein_dict:
            merged_df.loc[ind, 'mutation_type'] = protein_dict[row['mutation']][0]
        
        freq1 = float(row["new_freq_x"])
        cov1 = int(row["coverage_x"])
        bc1 = int(row["base_count_x"])
        freq2 = float(row["new_freq_y"])
        cov2 = int(row["coverage_y"])
        bc2 = int(row["base_count_y"])
        diff = abs(freq1 - freq2)
        
        if (freq1 == -1) or (freq2 == -1): # If one of the replicate's frequnecies NA
            merged_df.loc[ind, "final_freq"] = -1
        else:
            if (freq1 == 0) and (freq2 == 0): # If both replicate's frequnecies are zero
                merged_df.loc[ind, "final_freq"] = 0
                merged_df.loc[ind, "tot_cov"] = cov1 + cov2
                merged_df.loc[ind, "tot_base_count"] = bc1 + bc2
            elif (freq1 != 0) and (freq2 != 0): # If both replicate's frequnecies are not zero
                if (freq1 < 0.5) and (freq2 < 0.5): # If both replicate's frequnecies are low difference limit is 0.1, else 0.3
                    diff_limit = 0.1
                else:
                    diff_limit = 0.3
                
                if diff <= diff_limit: # Check if frequnecies diff in limit and calculate weighted avg, else NA
                    merged_df.loc[ind, "final_freq"] = calc_weighted_avg(row["base_count_x"], row["base_count_y"], cov1, cov2)
                    merged_df.loc[ind, "tot_cov"] = cov1 + cov2
                    merged_df.loc[ind, "tot_base_count"] = bc1 + bc2
                else:
                    merged_df.loc[ind, "final_freq"] = -1
            else: # If one replicate frequnecy is zero and the other one isn't
                merged_df.loc[ind, "final_freq"] = -1

    return merged_df

def filter(tsv1, tsv2, freq, coverage, base_count, protein_dict, result_dir, filter_indels):
    # Read file & add columns & filter mutations to each pair of replica's freqs file
    
    # Phase 0 filtering
    rep1_df_all = phase0_filter(tsv1, filter_indels)
    rep2_df_all = phase0_filter(tsv2, filter_indels)
    
    # Save all mutation except ref and PROBLEMATIC
    rep1_df = rep1_df_all[["ref_pos", "mutation", "base_count", "coverage", "frequency"]].copy()
    rep1_df.to_csv(f"{result_dir}/freq1_all_mutations.csv", index=False)
    rep2_df = rep2_df_all[["ref_pos", "mutation", "base_count", "coverage", "frequency"]].copy()
    rep2_df.to_csv(f"{result_dir}/freq2_all_mutations.csv", index=False)

    # Insert new frequency to each mutation according to each file independently
    rep1_df["new_freq"] = "new"
    rep2_df["new_freq"] = "new"

    new_rep1_df = phase1_filter(rep1_df, coverage, freq, base_count)
    new_rep2_df = phase1_filter(rep2_df, coverage, freq, base_count)
    
    new_rep1_df.to_csv(f"{result_dir}/freq1_independent_filtering.csv", index=False)
    new_rep2_df.to_csv(f"{result_dir}/freq2_independent_filtering.csv", index=False)
    
    # Insert new frequency to each mutation according to both freqs
    merged_df = phase2_filter(new_rep1_df, new_rep2_df, protein_dict)
    if "new" in merged_df["final_freq"].unique():
        print("***Error***\nIn decision tree phase 1, there is a case that's not covered!!!")
    merged_df.to_csv(f"{result_dir}/merged.csv", index=False)

    # Create DF for next phase (BN algorithem)
    final_df = merged_df[["mutation", "final_freq", "tot_cov", "tot_base_count"]]
    final_df.to_csv(f"{result_dir}/frequnecies.csv", index=False)
    
    # Return number of each frequncy type (totatl, NA, 0, W_avg)
    return final_df.shape[0], (final_df["final_freq"] == -1).sum(), (final_df["final_freq"] == 0).sum(), (final_df["final_freq"] > 0).sum()