import os
import subprocess
import time
from datetime import datetime

def get_res_dir():
    while True:
        user_input = input("Enter 1 for cluster and 2 for local run: ")
        if (user_input == "1"):
            RESULTS = r"/sternadi/home/volume1/ido/BN-SCRIPTS/BN_Create_input/results"
            break
        elif (user_input == "2"):
            RESULTS = r"Z:/home/volume1/ido/BN-SCRIPTS/BN_Create_input/results"
            break
        else:
            user_input = input("Wrong input please try again\nEnter 1 for cluster and 2 for local run: ")
    return RESULTS

def get_latest_res_dir(root):
    """
    input: root folder
    output: latest modified subdirectory from level 1
    """
    newest_modified_time = None
    newest_modified_directory = ""
    for subdir in os.listdir(root):
        full_path = os.path.join(root, subdir)
        modified_time = os.path.getmtime(full_path)
        if (newest_modified_time is None) or (modified_time > newest_modified_time):
                newest_modified_directory = full_path
                newest_modified_time = modified_time
    return newest_modified_directory

def choose_method(path):
    approx_command = f"Rscript {path}/BN-SCRIPTS/Run_BN_script/BB_bottleneck-master/Bottleneck_size_estimation_approx.r --file"
    exact_command = f"Rscript {path}/BN-SCRIPTS/Run_BN_script/BB_bottleneck-master/Bottleneck_size_estimation_exact.r --file"
    while True:
        method = input("1- For the approximate code run\n2- For the exact code run\nchoose method: ")
        if method == "1":
            return approx_command
        elif method == "2":
            return exact_command
        else:
            continue

def main():
    print("\nRunning R Bottleneck estimation script on preparation files for latest results directory in ./BN_Create_input/results")
    run_start = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    print(run_start)
    start = time.time()

    try:
        # Get latest results directory
        RESULTS = get_res_dir()
        results_dir = get_latest_res_dir(RESULTS)
        path = RESULTS[:-35]
        
        # Choose script method (approximate or excat)
        method = choose_method(path)
        
        # Change parameters if needed
        params = "--plot_bool TRUE --var_calling_threshold 0.03 --Nb_min 1 --Nb_max 200 --Nb_increment 1 --confidence_level .95"
        print(f"\nDefault command parameters:")
        print(params)
        
        while True:
            change = input("Do you want to change any input parameters (y/n): ")
            if change == "y":
                Nb_max = input("Enter Nb_max value: ")
                params = f"--plot_bool TRUE --var_calling_threshold 0.03 --Nb_min 1 --Nb_max {Nb_max} --Nb_increment 1 --confidence_level .95"
                break
            elif change == "n":
                break
            else:
                continue
        
        print(f"\nRESULTS DIR: {results_dir}")
                
        # Loop through patients
        for patient in os.listdir(results_dir):
            print(f"-----------PATIENT {patient[-2:]}-----------")
            txt_files = [f'{results_dir}/{patient}/{file}' for file in os.listdir(f'{results_dir}/{patient}') if file.endswith('.txt')]
            
            # Run R for file per patient
            for file in txt_files:
                command_line = f'{method} "{file}" {params}'
                process = subprocess.Popen(command_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                # Get the output and errors if any
                output, error = process.communicate()

                # Check for errors
                if error:
                    print(f"Error occurred:\n {error.decode('utf-8')}")
                else:
                    print(f"R script executed successfully:\n {output.decode('utf-8')}")
    
    except Exception as e:
        print("An error has occured!\nTerminating script...")
        print(f"Main script elapsed time: {(time.time() - start)} sec")
        print("Log:")
        print(e)
        exit(1)

    print("***BN estimation Script finished!***")
    print(f"BN estimation Script elapsed time: {(time.time() - start)} sec")

    return