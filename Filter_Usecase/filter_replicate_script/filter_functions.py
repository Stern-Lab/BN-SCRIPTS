from datetime import datetime
import pandas as pd
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

def filter_df(df, freq, coverage, base_count):
    """
    Filters a dataframe of mutations according to the threshold arguments provided at the top and according to
    the list of known problematic position (provided by the PROBLEMATIC assignment at the top).
    :param df: Dataframe to be filtered by received frequency coverage and base_count thresholds.
    :param freq: A frequency threshold according to which the received dataframe is to be filtered.
    :param coverage: A coverage threshold according to which the received dataframe is to be filtered.
    :param base_count: A base_count threshold according to which the received dataframe is to be filtered.
    :return: A new Dataframe identical to the received one, filtered according to the received value thresholds.
    """
    ret_df = df.loc[(df['frequency'] >= freq) & (df['ref_base'] != df['read_base']) &
                (df['coverage'] >= coverage) & (df['base_count'] >= base_count)]
    return ret_df[~ret_df['ref_pos'].astype(int).isin(PROBLEMATIC)]

def filter(tsv1, tsv2, patient, timepoint, freq, coverage, base_count, protein_dict):
    """
    Filters a dataframe of mutations according to the threshold arguments provided and according to
    the list of known problematic position (provided by the PROBLEMATIC assignment at the top).
    :param df: Dataframe to be filtered by received frequency coverage and base_count thresholds.
    :param freq: A frequency threshold according to which the received dataframe is to be filtered.
    :param coverage: A coverage threshold according to which the received dataframe is to be filtered.
    :param base_count: A base_count threshold according to which the received dataframe is to be filtered.
    :return: A new Dataframe identical to the received one, filtered according to the received value thresholds.
    """
    # Creates results directory
    date_time_str = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    patient_dir = r"./Filter_Usecase/results/" + patient
    if not os.path.exists(patient_dir):
        os.makedirs(patient_dir)
    run_dir = f"results_({freq}_{coverage}_{base_count})_{date_time_str}"
    res_dir = f"{patient_dir}/{timepoint}/{run_dir}"
    os.makedirs(res_dir)

    # Read file & add columns & filter mutations to each pair of replica's freqs file
    
    # All mutation WITHOUT FILTERING
    rep1_df_all = pd.read_csv(tsv1, sep='\t')
    rep1_df_all = filter_df(rep1_df_all, 0, 0, 1)
    rep1_df_all = enrich_mutation(rep1_df_all)
    rep1_df_all = filter_ref(rep1_df_all)

    rep2_df_all = pd.read_csv(tsv2, sep='\t')
    rep2_df_all = filter_df(rep2_df_all, 0, 0, 1)
    rep2_df_all = enrich_mutation(rep2_df_all)
    rep2_df_all = filter_ref(rep2_df_all)
    
    num_of_mut_rep1 = rep1_df_all.shape[0]
    num_of_mut_rep2 = rep2_df_all.shape[0]

    
    # All mutation WITH FILTERING
    rep1_df = filter_df(rep1_df_all, freq, coverage, base_count)
    rep1_df = enrich_mutation(rep1_df)
    rep1_df = filter_ref(rep1_df)
    rep1_df.to_csv(f"{res_dir}/replicate1.csv", index=False)

    rep2_df = filter_df(rep2_df_all, freq, coverage, base_count)
    rep2_df = enrich_mutation(rep2_df)
    rep2_df = filter_ref(rep2_df)
    rep2_df.to_csv(f"{res_dir}/replicate2.csv", index=False)
    
    # Create DF for next phase (BN algorithem)
    filtered_rep1 = rep1_df[["mutation", "frequency"]]
    filtered_rep1.to_csv(f"{res_dir}/mut_freq_1.csv", index=False)
    
    filtered_rep2 = rep2_df[["mutation", "frequency"]]
    filtered_rep2.to_csv(f"{res_dir}/mut_freq_2.csv", index=False)

    # Merge data framse based on common mutation after filtering
    merged_df = pd.merge(rep1_df, rep2_df, how="inner", on="mutation")
    merged_all_df = pd.merge(rep1_df_all, rep2_df_all, how="inner", on="mutation")
    num_of_mut_merged = merged_df.shape[0]
    num_of_mut_merged_all = merged_all_df.shape[0]
    
    # Add protein (mutation type) column
    merged_df["mutation_type"] = ""

    # More filtering based on both files
    merged_df["CriticalDelta"] = "No"
    merged_df["UseCaseGroup"] = 0

    # Usecase submitting
    for ind, row in merged_df.iterrows():
        rep1_freq = row['frequency_x']
        rep2_freq = row['frequency_y']
        if row['mutation'] in protein_dict:
            merged_df.loc[ind, 'mutation_type'] = protein_dict[row['mutation']][0]

        big_freq = max(rep1_freq, rep2_freq)
        small_freq = min(rep1_freq, rep2_freq)

        # Usecase 1
        if (0.8 <= big_freq) and (0.5 <= small_freq):
            merged_df.loc[ind, 'UseCaseGroup'] = 1
            continue
        
        # Usecase 2
        elif (0.8 <= big_freq) and (small_freq < 0.5):
            merged_df.loc[ind, 'CriticalDelta'] = "Yes"
            merged_df.loc[ind, 'UseCaseGroup'] = 2
        
        # Usecase 3
        elif (0.5 <= big_freq < 0.8) and (0.5 <= small_freq < 0.8):
            merged_df.loc[ind, 'UseCaseGroup'] = 3
            continue
        
        elif (0.5 <= big_freq < 0.8) and (small_freq < 0.5):
            # Usecase 4
            if abs(big_freq - small_freq) < 0.3:
                merged_df.loc[ind, 'UseCaseGroup'] = 4
            # Usecase 7
            else:
                merged_df.loc[ind, 'UseCaseGroup'] = 7
            continue

        elif (big_freq < 0.5) and (small_freq < 0.5):
            # Usecase 5
            if abs(big_freq - small_freq) >= 0.1:
                merged_df.loc[ind, 'UseCaseGroup'] = 5
            # Usecase 6
            else:
                merged_df.loc[ind, 'UseCaseGroup'] = 6
            continue
        
        # No appropriate usecase
        else:
            continue
    
    merged_df.to_csv(f"{res_dir}/merged.csv", index=False)
    usecase_df = merged_df[["ref_pos_x", "mutation", "mutation_type", "base_count_x", "base_count_y", "coverage_x", "coverage_y", "frequency_x", "frequency_y", "CriticalDelta", "UseCaseGroup"]]
    usecase_df.to_csv(f"{res_dir}/usecase.csv", index=False)
    output_df = merged_df[["mutation", "frequency_x", "frequency_y"]]
    output_df.to_csv(f"{res_dir}/output3.csv", index=False)
    
    return usecase_df, num_of_mut_rep1, num_of_mut_rep2, num_of_mut_merged, num_of_mut_merged_all