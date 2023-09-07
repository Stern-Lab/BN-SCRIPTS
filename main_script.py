import filter_replicate_script.filter_functions as ff
import filter_replicate_script.protein_data as prot
import time
from datetime import datetime
import pandas as pd
import numpy as np
import os

FREQ = 0.01
COVERAGE = 100
BASECOUNT = 50

def get_freq_file(sample_id, v3_dirs, v4_dirs):
    if sample_id in v3_dirs:
        freq_file = PATH_V3 + "/" + sample_id + "/freqs.tsv"
        return freq_file, True

    elif sample_id in v4_dirs:
        freq_file = PATH_V4 + "/" + sample_id + "/freqs.tsv"
        return freq_file, True
    
    else:
        print(f"Directory {sample_id} wasn't found in V3 or V4. Skipping...")
        return "", False

def usecase_calc(uc_df, res_df, cols_list, ind, mut_cnt):
    for col in cols_list:
        res_df.loc[ind, col] = 0

    for key in uc_df["UseCaseGroup"].value_counts().index.tolist():
        res_df.loc[ind, cols_list[key-1]] = (uc_df["UseCaseGroup"].value_counts()[key])/mut_cnt

    res_df.loc[ind, "criticalDelta_cnt"] = len(uc_df[uc_df["CriticalDelta"] != "No"])
    return res_df

def get_bn_prep_files(patient_id, timepoint):
    newest_modified_time = None
    newest_modified_directory = ""
    for dir, _, files in os.walk(f"./result/{patient_id}/{timepoint}"):
            modified_time = os.path.getmtime(dir)
            if (newest_modified_time is None) or (modified_time > newest_modified_time):
                    newest_modified_directory = dir
                    newest_modified_time = modified_time
    return f"./result/{patient_id}/{timepoint}/" + newest_modified_directory + f"/output3_{FREQ}_{COVERAGE}_{BASECOUNT}.csv"

if __name__ == "__main__":
    print("Main script is starting...")
    print(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    start = time.time()
    try:
        while True:
            user_input = input("Enter 1 for cluster and 2 for local run: ")
            # Prepartion actions
            if (user_input == "1"):
                PATH_V3 = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V3"
                PATH_V4 = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V4"
                PATIENTS = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/all_patients_global_content_initials.csv"
                break
            elif (user_input == "2"):
                PATH_V3 = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V3"
                PATH_V4 = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V4"
                PATIENTS = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/all_patients_global_content_initials.csv"
                break
            else:
                user_input = input(
                    "Wrong input please try again\nEnter 1 for cluster and 2 for local run: ")
        
        protein_dict = prot.create_protein_dict() # Protein dictionary
        all_patients_df = pd.read_csv(PATIENTS) # Load patients data
        sample_size = all_patients_df.shape[0]
        res_cols = ["Patient", "Timepoint", "sample_ID", "Date","Ct", "merged_total_mutations", "criticalDelta_cnt",
                    "UC1_freq", "UC2_freq", "UC3_freq", "UC4_freq", "UC5_freq", "UC6_freq"]
        results_df = pd.DataFrame(columns=res_cols) # Create results data frame
        
        # Get list of all directories
        v3_dirs = os.listdir(PATH_V3)
        v4_dirs = os.listdir(PATH_V4)

        # Update results data frame with all inforamtion needed
        ind = 0
        for i, curr_row in all_patients_df.iterrows():

            print(f"Progress: {(i/sample_size*100):.2f}%")

            # Check if sample is not nan
            if (np.isnan(curr_row["sample"])):
                print("Sample id is nan. Skipping iteration...")
                continue

            # Get patient info
            curr_sample_id = str(int(curr_row["sample"]))
            curr_patient_id = curr_row["patient_ID"][:2]
            curr_timepoint = int(curr_row["time_since_first_sampling"])

            print(f"Filtering patient {curr_patient_id} timepoint {curr_timepoint}.")
            
            # Find replicate1 file
            s1_rep1, found = get_freq_file(curr_sample_id, v3_dirs, v4_dirs)
            if not (found):
                 continue
            # Find replicate2 file
            s1_rep2, found = get_freq_file(curr_sample_id + "_L001", v3_dirs, v4_dirs)
            if not (found):
                 continue
            
            ind += 1
            # Update patients data to results_df
            results_df.loc[ind, "Timepoint"] = curr_timepoint
            results_df.loc[ind, "Patient"] = curr_patient_id
            results_df.loc[ind, "sample_ID"] = curr_sample_id
            results_df.loc[ind, "Date"] = curr_row["sampling date"]
            results_df.loc[ind, "Ct"] = curr_row["mean ct"]
            
            # filter timepoint and update data to results
            usecase_df = ff.filter(s1_rep1, s1_rep2, curr_patient_id, curr_timepoint, FREQ, COVERAGE, BASECOUNT, protein_dict)
            mutation_cnt = usecase_df.shape[0]
            results_df.loc[ind, "merged_total_mutations"] = mutation_cnt
            
            # Update usecases statistics
            us_col_list = ["UC1_freq", "UC2_freq", "UC3_freq", "UC4_freq", "UC5_freq", "UC6_freq"]
            results_df = usecase_calc(usecase_df, results_df, us_col_list, ind, mutation_cnt)

            # # Check if previous sample has a pair to compare
            # if (ind == 1): # Skip first index
            #     continue
            
            # # Get previus sample information
            # prev_row = results_df.loc[ind-1]
            # prev_sample_id = str(int(prev_row["sample_ID"]))
            # prev_patient_id = prev_row["Patient"]
            # prev_timepoint = int(prev_row["Timepoint"])

            # # Check if two samples are from same patient
            # if curr_patient_id != prev_patient_id:
            #     continue
            
            # # Create BN input files
            # print(f"Creating bottlneck input files for {curr_patient_id} at time point {prev_timepoint} and time point {curr_timepoint}")
            
            # # Find input files (filter script output files) for prev row
            # prev_bn_file = get_bn_prep_files(prev_patient_id, prev_timepoint)

            
        # Count non-zero parameters
        results_size = results_df.shape[0] + 1
        results_df.loc[results_size, "Patient"] = "COUNT NON-ZERO"
        us_col_list.append("criticalDelta_cnt")
        for col in (us_col_list):
            results_df.loc[results_size, col] = len(results_df[results_df[col] != 0])
        
        # Save data frame as a file
        results_df.to_csv("./Results.csv", index=False)
    
    except Exception as e:
        print("An error has occured!\nTerminating script...")
        print(f"Main script elapsed time: {(time.time() - start)} sec")
        print("Log:")
        print(e)
        exit(1)

    print("***Script finished!***")
    print(f"Script elapsed time: {(time.time() - start)} sec")
