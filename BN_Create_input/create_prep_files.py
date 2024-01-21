import time
from datetime import datetime
import pandas as pd
import os

def get_res_dir():
    while True:
        # user_input = input("Enter 1 for cluster and 2 for local run: ")
        user_input = "1"
        if (user_input == "1"):
            RESULTS = r"/sternadi/home/volume1/ido/BN-SCRIPTS/Filter_Usecase/results"
            break
        elif (user_input == "2"):
            RESULTS = r"Z:\home\volume1\ido\BN-SCRIPTS\Filter_Usecase\results"
            break
        else:
            user_input = input("Wrong input please try again\nEnter 1 for cluster and 2 for local run: ")
    return RESULTS 

def get_bn_prep_files(root, patient_id, timepoint):
    """
    input: root folder
    output: latest modified subdirectory from level 1
    """
    newest_modified_time = None
    newest_modified_directory = ""
    # Get a list of all items (files and directories) in the specified directory
    all_items = os.listdir(root)
    # Filter out only directories
    directories = [item for item in all_items if os.path.isdir(os.path.join(root, item))]
    for subdir in directories:
        full_path = os.path.join(root, subdir)
        modified_time = os.path.getmtime(full_path)
        if (newest_modified_time is None) or (modified_time > newest_modified_time):
                newest_modified_directory = full_path
                newest_modified_time = modified_time
    return f"{newest_modified_directory}/{patient_id}/{timepoint}/frequnecies.csv"

def main():
    print("BN prepartion files script is starting...")
    run_start = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    print(run_start)
    start = time.time()
    try:
        RESULTS = get_res_dir()
        
        all_samples = pd.read_csv(RESULTS + "/Results.csv")
        sample_size = all_samples.shape[0]

        # Update results data frame with all inforamtion needed
        ind = 0
        for i, curr_row in all_samples.iterrows():
            print(f"Progress: {(i/sample_size*100):.2f}%")            

            if(i==0):
                # Get previous timepoint patient info
                prev_patient_id = curr_row["Patient"]
                prev_timepoint = str(int(curr_row["Timepoint"]))
                continue

            # Get current patient info
            curr_patient_id = curr_row["Patient"]
            curr_timepoint = str(int(curr_row["Timepoint"]))

            # Check if samples's patient name are the same
            if (prev_patient_id != curr_patient_id):
                prev_patient_id = curr_patient_id
                prev_timepoint = curr_timepoint
                continue

            print(f"Creating prepration files for {curr_patient_id} for timepoints {prev_timepoint} and {curr_timepoint}.")
            
            # Create results directory
            res_dir = f"./BN_Create_input/results/{run_start}/{curr_patient_id}"
            if not os.path.exists(res_dir):
                os.makedirs(res_dir)
            
            # Find filtering results
            t1_results = get_bn_prep_files(RESULTS, curr_patient_id, prev_timepoint)            
            t1_df = pd.read_csv(t1_results)
            t2_results = get_bn_prep_files(RESULTS, curr_patient_id, curr_timepoint)            
            t2_df = pd.read_csv(t2_results)

            # Keep all mutations that appears in one of the time points
            merged_df = pd.merge(t1_df, t2_df, how='outer', on= 'mutation') # Ask Natalie
            
            # Re-organize df (change column's name)
            merged_df.rename(columns={"final_freq_x":f"frequency_{prev_timepoint}", "final_freq_y":f"frequency_{curr_timepoint}"}, inplace=True)

            # Save df
            merged_df.to_csv(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/merged_{prev_timepoint}_{curr_timepoint}.csv", index=False)
    
            # Drop mutations that's not in user's usecase choice
            merged_df = merged_df[(merged_df[f'frequency_{prev_timepoint}'] != -1) & (merged_df[f'frequency_{curr_timepoint}'] != -1)] # Drop NA mutations
            merged_df = merged_df[(merged_df[f'frequency_{prev_timepoint}'] != 0) | (merged_df[f'frequency_{curr_timepoint}'] != 0)] # Drop both zero mutations
            merged_df.to_csv(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/res_{prev_timepoint}_{curr_timepoint}.csv", index=False)
    
            # Create text file for Bottleneck algorithm
            txt = ""
            for i, row in res_df.iterrows():
                tot_coverage = str.format('{0:.6f}', row['tot_cov_x'])
                tot_base_count = str.format('{0:.6f}', row['tot_cov_y'])

                txt += (str.format('{0:.6f}', row[f'frequency_{prev_timepoint}']) + "\t" + str.format('{0:.6f}',f'frequency_{curr_timepoint}') + "\t" + tot_coverage + "\t" + tot_base_count + "\n")

            with open(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/frequencies_{prev_timepoint}_{curr_timepoint}.txt", 'w') as file:
                file.write(txt)

            # Next time point
            prev_patient_id = curr_patient_id
            prev_timepoint = curr_timepoint

    except Exception as e:
        print("An error has occured!\nTerminating script...")
        print(f"Main script elapsed time: {(time.time() - start)} sec")
        print("Log:")
        print(e)
        exit(1)

    print("***BN prepration files Script finished!***")
    print(f"Script elapsed time: {(time.time() - start)} sec")
