import filter_replicate_script.filter_functions as ff
import filter_replicate_script.protein_data as prot
import time
import pandas as pd
import numpy as np
import os

PATH_V3 = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V3"
PATH_V4 = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V4"
PATIENTS = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/all_patients_global_content_initials.csv"

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
            
            


if __name__ == "__main__":
    print("Main script is starting...")
    start = time.time()
    try:
        # Prepartion actions
        # protein_dict = prot.create_protein_dict() # Protein dictionary
        all_patients_df = pd.read_csv(PATIENTS) # Load patients data
        sample_size = all_patients_df.shape[0]
        res_cols = ["Timepoints", "Patient", "sample_ID_1", "sample_ID_2", "Date_1", "Date_2", "merged_Total_mutations",
                    "UC1_freq", "UC2_freq", "UC3_freq", "UC4_freq", "UC5_freq", "UC6_freq"]
        results_df = pd.DataFrame(columns=res_cols) # Create results data frame
        # Get list of all directories
        v3_dirs = os.listdir(PATH_V3)
        v4_dirs = os.listdir(PATH_V4)

        # Update results data frame with all inforamtion needed
        for ind, curr_row in all_patients_df.iterrows():
            if (ind == sample_size-1):
                continue

            next_row = all_patients_df.loc[ind+1]
            print(f"Progress: {(ind/sample_size*100):.2f}%")

            # Check if two samples are related to same patient
            if (np.isnan(curr_row["sample"]) or np.isnan(next_row["sample"])):
                print("Sample id is nan. Skipping iteration...")
                continue

            # Get patient info
            curr_sample_id = str(int(curr_row["sample"]))
            next_sample_id = str(int(next_row["sample"]))
            curr_patient_id = curr_row["patient_ID"][:2]
            next_patient_id = next_row["patient_ID"][:2]
            curr_timepoint = int(curr_row["time_since_first_sampling"])
            next_timepoint = int(next_row["time_since_first_sampling"])


            # Check if two samples are from same patient
            if curr_patient_id != next_patient_id:
                continue

            print(f"Processing patient {curr_patient_id} timepoint {curr_timepoint} and timepoint {next_timepoint}")
            
            # Update patients data to results_df
            results_df.loc[ind, "Timepoints"] = f"TP1_{curr_timepoint}_TP2_{next_timepoint}"
            results_df.loc[ind, "Patient"] = curr_patient_id
            results_df.loc[ind, "sample_ID_1"] = curr_sample_id
            results_df.loc[ind, "sample_ID_2"] = next_sample_id
            results_df.loc[ind, "Date_1"] = curr_row["sampling date"]
            results_df.loc[ind, "Date_2"] = next_row["sampling date"]

            # Find sample1/replicate1 file
            s1_rep1, found = get_freq_file(curr_sample_id, v3_dirs, v4_dirs)
            if not(found):
                 continue
            # Find sample1/replicate2 file
            s1_rep2, found = get_freq_file(curr_sample_id + "_L001", v3_dirs, v4_dirs)
            if not (found):
                 continue
            
            # filter time point 1

            # Find sample2/replicate1 file
            s2_rep1 = get_freq_file(next_sample_id, v3_dirs, v4_dirs)
            if not (found):
                 continue
            # Find sample2/replicate2 file
            s2_rep2 = get_freq_file(next_sample_id + "_L001", v3_dirs, v4_dirs)
            if not (found):
                 continue
        
        results_df.to_csv("./Results.csv", index=False)
    
    except Exception as e:
        print("An error has occured!\nTerminating script...")
        print(f"Main script elapsed time: {(time.time() - start)} sec")
        print("Log:")
        print(e)
        exit(1)

    print("Script finished!")
    print(f"Script elapsed time: {(time.time() - start)} sec")
