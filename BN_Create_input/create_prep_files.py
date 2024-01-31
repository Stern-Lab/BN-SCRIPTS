import time
from datetime import datetime
import pandas as pd
import os

def get_res_dir():
    while True:
        user_input = input("Enter 1 for cluster and 2 for local run: ")
        # user_input = "1"
        if (user_input == "1"):
            RESULTS = r"/sternadi/home/volume1/ido/BN-SCRIPTS/Filter_Usecase/results"
            break
        elif (user_input == "2"):
            RESULTS = r"Z:\home\volume1\ido\BN-SCRIPTS\Filter_Usecase\results"
            break
        else:
            user_input = input("Wrong input please try again\nEnter 1 for cluster and 2 for local run: ")
    return RESULTS 

def results_dir_choice(res_root_dir):
    # Get a list of all items (files and directories) in the specified directory
    all_items = os.listdir(res_root_dir)
    # Filter out only directories
    directories = [item for item in all_items if os.path.isdir(os.path.join(res_root_dir, item))]
    user_txt = "\nResults directories:\n"
    for i, res_dir in enumerate(directories):
        user_txt += f"{i+1}) {res_dir}\n"
    
    print(user_txt)
    res_dir_ind = input("Choose results directory index: ")

    return directories[int(res_dir_ind)-1]

def get_bn_prep_files(dir, patient_id, timepoint):
    """
    input: root folder
    output: latest modified subdirectory from level 1
    """
    return f"{dir}/{patient_id}/{timepoint}/{patient_id}_T{timepoint}_frequnecies.csv"

def main():
    print("BN prepartion files script is starting...")
    run_start = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    log_txt = f"{run_start}\n"
    print(run_start)
    start = time.time()
    try:
        # Get results directory
        RESULTS = get_res_dir()
        user_res_dir = results_dir_choice(RESULTS)

        log_txt += f"Filtering results input path: {RESULTS}/{user_res_dir}\n"
        log_txt += f"Results output path: ./BN_Create_input/results/{run_start}\n"

        
        all_samples = pd.read_csv(fr"{RESULTS}/{user_res_dir}/Results.csv")
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
            log_txt += f"Creating prepration files for {curr_patient_id} for timepoints {prev_timepoint} and {curr_timepoint}.\n"
            
            # Create results directory
            res_dir = f"./BN_Create_input/results/{run_start}/{curr_patient_id}"
            if not os.path.exists(res_dir):
                os.makedirs(res_dir)
            
            # Find filtering results
            t1_results = get_bn_prep_files(fr"{RESULTS}/{user_res_dir}", curr_patient_id, prev_timepoint)  
            t1_df = pd.read_csv(t1_results)
            t2_results = get_bn_prep_files(fr"{RESULTS}/{user_res_dir}", curr_patient_id, curr_timepoint)            
            t2_df = pd.read_csv(t2_results)

            # Keep all mutations that appears in both of the time points
            merged_df = pd.merge(t1_df, t2_df, how='inner', on= 'mutation')
            
            # Re-organize df (change column's name)
            merged_df.rename(columns={"final_freq_x":f"frequency_{prev_timepoint}", "final_freq_y":f"frequency_{curr_timepoint}"}, inplace=True)
            merged_df.rename(columns={"tot_cov_x":f"tot_cov_{prev_timepoint}", "tot_cov_y":f"tot_cov_{curr_timepoint}"}, inplace=True)
            merged_df.rename(columns={"tot_base_count_x":f"tot_base_count_{prev_timepoint}", "tot_base_count_y":f"tot_base_count_{curr_timepoint}"}, inplace=True)

            # Save df
            merged_df.to_csv(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/merged_{prev_timepoint}_{curr_timepoint}.csv", index=False)
    
            # Drop mutations that's not in user's usecase choice
            merged_df = merged_df[(merged_df[f'frequency_{prev_timepoint}'] != -1) & (merged_df[f'frequency_{curr_timepoint}'] != -1)] # Keep mutations without NA
            merged_df = merged_df[(merged_df[f'frequency_{prev_timepoint}'] != 0) | (merged_df[f'frequency_{curr_timepoint}'] != 0)] # Keep mutations if one of the freqs different than 0
            merged_df.to_csv(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/final_{prev_timepoint}_{curr_timepoint}.csv", index=False)
    
            # Create text file for Bottleneck algorithm
            txt = ""
            for i, row in merged_df.iterrows():
                t2_tot_coverage = str.format('{0:.6f}', row[f'tot_cov_{curr_timepoint}'])
                t2_tot_base_count = str.format('{0:.6f}', row[f'tot_base_count_{curr_timepoint}'])

                txt += (str.format('{0:.6f}', row[f'frequency_{prev_timepoint}']) + "\t" + str.format('{0:.6f}',row[f'frequency_{curr_timepoint}']) + "\t" + t2_tot_coverage + "\t" + t2_tot_base_count + "\n")

            with open(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/frequencies_{prev_timepoint}_{curr_timepoint}.txt", 'w') as file:
                file.write(txt)

            # Next time point
            prev_patient_id = curr_patient_id
            prev_timepoint = curr_timepoint
        
        tot_time = time.time() - start
        # Save log data as a file
        print("Saving log.txt file...")
        with open(f"./BN_Create_input/results/{run_start}/log.txt", 'w') as log_file:
            log_txt += f"Script elapsed time: {tot_time} sec"
            log_file.write(log_txt)
        
        print("***BN prepration files Script finished!***")
        print(f"Script elapsed time: {(tot_time)} sec")

    except Exception as e:
        print("An error has occured!\nTerminating script...")
        print(f"Main script elapsed time: {(time.time() - start)} sec")
        print("Log:")
        print(e)
        # Save log data as a file
        with open(f"./BN_Create_input/results/{run_start}/log.txt", 'w') as log_file:
            log_txt += f"{e}\n"
            log_txt += f"Script elapsed time: {tot_time} sec"
            log_file.write(log_txt)
        exit(1)

