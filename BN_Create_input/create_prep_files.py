import time
from datetime import datetime
import pandas as pd
import os

def get_res_dir():
    while True:
        user_input = input("Enter 1 for cluster and 2 for local run: ")
        if (user_input == "1"):
            RESULTS = r"/sternadi/home/volume1/ido/BN-SCRIPTS/Filter_Usecase"
            break
        elif (user_input == "2"):
            RESULTS = r"Z:/home/volume1/ido/BN-SCRIPTS/Filter_Usecase"
            break
        else:
            user_input = input("Wrong input please try again\nEnter 1 for cluster and 2 for local run: ")
    return RESULTS 

def get_bn_prep_files(directory, patient_id, timepoint):
    newest_modified_time = None
    newest_modified_directory = ""
    for dir, _, files in os.walk(f"{directory}/results/{patient_id}/{timepoint}"):
            modified_time = os.path.getmtime(dir)
            if (newest_modified_time is None) or (modified_time > newest_modified_time):
                    newest_modified_directory = dir
                    newest_modified_time = modified_time
    return f"{newest_modified_directory}/frequnecies.csv"

def calc_weighted_avg(bs1, bs2, cvg1, cvg2):
    return (bs1 + bs2) / (cvg1 + cvg2)

def main():
    print("BN prepartion files script is starting...")
    run_start = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    print(run_start)
    start = time.time()
    try:
        RESULTS = get_res_dir()
        
        all_samples = pd.read_csv(f"{RESULTS}/Results.csv")
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
            
            # Find filtering and usecase results
            t1_results = get_bn_prep_files(RESULTS, curr_patient_id, prev_timepoint)            
            t1_df = pd.read_csv(t1_results)
            t2_results = get_bn_prep_files(RESULTS, curr_patient_id, curr_timepoint)            
            t2_df = pd.read_csv(t2_results)

            # Keep mutations that appears in both time points
            merged_df = pd.merge(t1_df, t2_df, how= 'inner', on= 'mutation')

            # Re-organize DF (change column's name)
            merged_df.drop(['ref_pos_x_y', 'mutation_type_y'], inplace=True, axis=1)
            merged_df.rename(columns={"ref_pos_x_x": "ref_pos", "mutation_type_x": "mutation_type", "base_count_x_x":f"base_count_{prev_timepoint}_1",
                                    "base_count_y_x":f"base_count_{prev_timepoint}_2", "coverage_x_x":f"coverage_{prev_timepoint}_1", "coverage_y_x":f"coverage_{prev_timepoint}_2",
                                    "frequency_x_x":f"frequency_{prev_timepoint}_1", "frequency_y_x":f"frequency_{prev_timepoint}_2"}, inplace=True)
            merged_df.rename(columns={"base_count_x_y":f"base_count_{curr_timepoint}_1", "base_count_y_y":f"base_count_{curr_timepoint}_2",
                                    "coverage_x_y":f"coverage_{curr_timepoint}_1", "coverage_y_y":f"coverage_{curr_timepoint}_2",
                                    "frequency_x_y":f"frequency_{curr_timepoint}_1", "frequency_y_y":f"frequency_{curr_timepoint}_2"}, inplace=True)

            # Save DF
            res_dir = f"./BN_Create_input/results/{run_start}/{curr_patient_id}"
            if not os.path.exists(res_dir):
                os.makedirs(res_dir)
            merged_df.to_csv(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/merged_{prev_timepoint}_{curr_timepoint}.csv", index=False)

            res_df = merged_df[['ref_pos', 'mutation', 'mutation_type', f'base_count_{prev_timepoint}_1', f'base_count_{prev_timepoint}_2',
                                f'base_count_{curr_timepoint}_1', f'base_count_{curr_timepoint}_2', f'frequency_{prev_timepoint}_1',
                                f'frequency_{prev_timepoint}_2', f'frequency_{curr_timepoint}_1', f'frequency_{curr_timepoint}_2',
                                f'coverage_{prev_timepoint}_1', f'coverage_{prev_timepoint}_2', f'coverage_{curr_timepoint}_1',
                                f'coverage_{curr_timepoint}_2', f'UseCaseGroup_x', 'UseCaseGroup_y']]
    
            # Drop mutations that's not in user's usecase choice
            res_df = res_df[(res_df['UseCaseGroup_x'].isin(uc_list)) & (res_df['UseCaseGroup_y'].isin(uc_list))]
            res_df.to_csv(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/res_{prev_timepoint}_{curr_timepoint}.csv", index=False)
            
            # Add columns for final file
            new_col1 = f"final_freq_{prev_timepoint}"
            new_col2 = f"final_freq_{curr_timepoint}"
            res_df[new_col1] = ""
            res_df[new_col2] = ""
            res_df['timepoint_2_coverage'] = ""
            res_df['timepoint_2_basecount'] = ""
            
            # Calculate frequencies by Usecase
            for i, row in res_df.iterrows():
                # Previous time point
                # Take weighted average for all usecases
                final_freq1 = calc_weighted_avg(row[f'base_count_{prev_timepoint}_1'], row[f'base_count_{prev_timepoint}_2'], row[f'coverage_{prev_timepoint}_1'], row[f'coverage_{prev_timepoint}_2'])
                
                # Current time point
                # Take weighted average for all usecases
                final_freq2 = calc_weighted_avg(row[f'base_count_{curr_timepoint}_1'], row[f'base_count_{curr_timepoint}_2'], row[f'coverage_{curr_timepoint}_1'], row[f'coverage_{curr_timepoint}_2'])

                coverage_sum2 = (row[f'coverage_{curr_timepoint}_1'] + row[f'coverage_{curr_timepoint}_2'])
                base_count_sum2 = (row[f'base_count_{curr_timepoint}_1'] + row[f'base_count_{curr_timepoint}_2'])
            
                res_df.loc[i, new_col1] = round(final_freq1, 6)
                res_df.loc[i, new_col2] = round(final_freq2, 6)
                res_df.loc[i, 'timepoint_2_coverage'] = coverage_sum2
                res_df.loc[i, 'timepoint_2_basecount'] = base_count_sum2

            # Save DF
            prep_df = res_df[['mutation', new_col1, new_col2, 'timepoint_2_coverage', 'timepoint_2_basecount']]
            prep_df.to_csv(f"./BN_Create_input/results/{run_start}/{curr_patient_id}/prep_{prev_timepoint}_{curr_timepoint}.csv", index=False)

            # Create text file for Bottleneck Algorithem
            txt = ""
            for i, row in prep_df.iterrows():
                tot_coverage = str.format('{0:.6f}', row['timepoint_2_coverage'])
                tot_base_count = str.format('{0:.6f}', row['timepoint_2_basecount'])

                txt += (str.format('{0:.6f}', row[new_col1]) + "\t" + str.format('{0:.6f}',row[new_col2]) + "\t" + tot_coverage + "\t" + tot_base_count + "\n")

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
