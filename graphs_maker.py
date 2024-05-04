import os
import pandas as pd
# import matplotlib.pyplot as plt
# from sklearn.linear_model import LinearRegression

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

def data_to_graph(x_axis, y_axis, p, tp, filtered):
    
    # Create scatter plot
    plt.scatter(x_axis, y_axis, color='blue')
                
    # Create regression line
    reg = LinearRegression()
    reg.fit(x_axis, y_axis)
    reg_line = reg.predict(x_axis)
    plt.plot(x_axis, reg_line, color='black', linestyle='--')
    plt.title(f'{p}_T{tp}')

    if filtered:
        plt.savefig(fr"{tp_path}/{p}_T{tp}_filtered.jpg")
    else:
        plt.savefig(fr"{tp_path}/{p}_T{tp}_not_filtered.jpg")

    return

if __name__ == "__main__":
    try:
        RESULTS = r'/sternadi/home/volume1/ido/BN-SCRIPTS/Filter_Usecase/results'
        GRAPHS_DIR = r"./graphs"
        res_dir = RESULTS + "/" + results_dir_choice(r'/sternadi/home/volume1/ido/BN-SCRIPTS/Filter_Usecase/results')
        patiens_dir = [d for d in os.listdir(res_dir) if (os.path.isdir(os.path.join(res_dir, d)))]
        for p in patiens_dir:
            print(f"\n***Patients {p}***")
            full_p_dir = res_dir + f"/{p}"
            time_points = [d for d in os.listdir(full_p_dir) if (os.path.isdir(os.path.join(full_p_dir, d)))]
            for tp in time_points:
                print(f"--TP {tp}--")
                tp_path = full_p_dir + f"/{tp}"

                df_file = fr"{tp_path}/{p}_T{tp}_merged.csv"

                # Not filtered graph
                merged_df = pd.read_csv(df_file)
                # data_to_graph(merged_df["frequency_x"], merged_df["frequency_y"], p, tp, filtered=False)
                
                # Filtered graph
                merged_df_filtered = merged_df[merged_df['final_freq'] > 0]
                # data_to_graph(merged_df_filtered["frequency_x"], merged_df_filtered["frequency_y"], p, tp, filtered=True)

                print(merged_df.shape[0], merged_df_filtered.shape[0])

        print("All patients Done!!")

    except Exception as e:
        print("ERROR!!!")
        print(e)
        exit(1)
