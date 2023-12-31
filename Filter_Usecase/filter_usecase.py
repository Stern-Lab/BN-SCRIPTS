import Filter_Usecase.filter_replicate_script.filter_functions as ff
import Filter_Usecase.filter_replicate_script.protein_data as prot
import time
from datetime import datetime
import pandas as pd
import numpy as np
import os

def get_freq_file(sample_id, v3_dirs, v4_dirs, PATH_V3, PATH_V4):
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

def main():
    print("Data filtering and Usecase table creation script is starting...")
    print(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    start = time.time()
    try:
        while True:
            user_input = input("Enter 1 for cluster and 2 for local run: ")
            
            # Prepartion actions
            if (user_input == "1"):
                PATH_V3 = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V3"
                PATH_V4 = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V4"
                PATIENTS = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_analysis/replicates_2023/all_patients_global_content_initials_V4.csv"
                break
            elif (user_input == "2"):
                PATH_V3 = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V3"
                PATH_V4 = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V4"
                PATIENTS = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_analysis/replicates_2023/all_patients_global_content_initials_V4.csv"
                break
            else:
                user_input = input("Wrong input please try again\nEnter 1 for cluster and 2 for local run: ")
        
        protein_dict = prot.create_protein_dict() # Protein dictionary
        all_patients_df = pd.read_csv(PATIENTS) # Load patients data
        sample_size = all_patients_df.shape[0]
        res_cols = ["Patient", "Timepoint", "sample_ID", "Date","Ct", "total_mutations_Before_Filtering_rep1", "total_mutations_Before_Filtering_rep2",
                    "total_mutations_Before_Filtering_merged", "total_mutations_After_Filtering_merged","criticalDelta_cnt", "UC1_freq", "UC2_freq",
                    "UC3_freq", "UC4_freq", "UC5_freq", "UC6_freq"]
        results_df = pd.DataFrame(columns=res_cols) # Create results data frame
        
        # Get list of all directories
        v3_dirs = os.listdir(PATH_V3)
        v4_dirs = os.listdir(PATH_V4)

        # Update results data frame with all inforamtion needed
        change_filter = input("Defaults filtering paramters are FREQ = 0.01, COVERAGE = 100, BASECOUNT = 50.\nDo you want to enter filtering parameters (y/n)? ")
        while True:
            if (change_filter == "n"):
                FREQ = 0.01
                COVERAGE = 100
                BASECOUNT = 50
                break
            elif (change_filter == "y"):
                FREQ = float(input("Minimal Frequnecy = "))
                COVERAGE = int(input("Minimal Coverage = "))
                BASECOUNT = int(input("Minimal Basecount = "))
                break
            else:
               change_filter = input("Wrong input!\nPlease enter (y/n): ") 

        ind = 0
        for i, curr_row in all_patients_df.iterrows():

            print(f"Progress: {(i/sample_size*100):.2f}%")

            # Check if sample is not nan
            if (np.isnan(curr_row["sample"])):
                print("Sample id is empty. Skipping iteration...")
                continue

            # Get patient info
            curr_sample_id = str(int(curr_row["sample"]))
            if (pd.isna(curr_row["patient_ID"])):
                print("Patient id is empty. Skipping iteration...")
                continue
            curr_patient_id = curr_row["patient_ID"][:2]

            if (np.isnan(curr_row["time_since_first_sampling"])):
                print("Sample time is empty. Skipping iteration...")
                continue
            curr_timepoint = int(curr_row["time_since_first_sampling"])

            print(f"Filtering patient {curr_patient_id} timepoint {curr_timepoint}.")
            
            # Find replicate1 file
            s1_rep1, found = get_freq_file(curr_sample_id, v3_dirs, v4_dirs, PATH_V3, PATH_V4)
            if not (found):
                continue
            
            # Find replicate2 file
            s1_rep2, found = get_freq_file(curr_sample_id + "_L001", v3_dirs, v4_dirs, PATH_V3, PATH_V4)
            if not (found):
                continue
            
            # Update index
            ind += 1
            
            # Update patients data to results_df
            results_df.loc[ind, "Timepoint"] = curr_timepoint
            results_df.loc[ind, "Patient"] = curr_patient_id
            results_df.loc[ind, "sample_ID"] = curr_sample_id
            results_df.loc[ind, "Date"] = curr_row["sampling date"]
            results_df.loc[ind, "Ct"] = curr_row["mean_ct"]
            
            # filter timepoint and update data to results
            usecase_df, num_of_mut_rep1, num_of_mut_rep2, num_of_mut_merged, num_of_mut_merged_all = ff.filter(s1_rep1, s1_rep2, curr_patient_id, curr_timepoint, FREQ, COVERAGE, BASECOUNT, protein_dict)
            results_df.loc[ind, "total_mutations_Before_Filtering_rep1"] = num_of_mut_rep1
            results_df.loc[ind, "total_mutations_Before_Filtering_rep2"] = num_of_mut_rep2
            results_df.loc[ind, "total_mutations_Before_Filtering_merged"] = num_of_mut_merged_all
            results_df.loc[ind, "total_mutations_After_Filtering_merged"] = num_of_mut_merged
            
            # Update usecases statistics
            us_col_list = ["UC1_freq", "UC2_freq", "UC3_freq", "UC4_freq", "UC5_freq", "UC6_freq"]
            results_df = usecase_calc(usecase_df, results_df, us_col_list, ind, num_of_mut_merged)
        
        # Save data frame as a file
        results_df.to_csv(r"./Filter_Usecase/Results.csv", index=False)
    
    except Exception as e:
        print("An error has occured!\nTerminating script...")
        print(f"Main script elapsed time: {(time.time() - start)} sec")
        print("Log:")
        print(e)
        exit(1)

    print("***Filter Script finished successfully!***")
    print(f"Script elapsed time: {(time.time() - start)} sec")
