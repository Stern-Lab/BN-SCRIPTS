import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import seaborn as sns
import scipy

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
    # res_dir_ind = 12
    return directories[int(res_dir_ind)-1]

def data_to_graph(df, ax, p, tp, cls):
    sns.scatterplot(data=df, x='frequency_x', y='frequency_y', ax=ax)
    x = df['frequency_x']
    y = df['frequency_y']
    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    r_squared = r_value ** 2
    r_x = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    r_y = [slope * xi + intercept for xi in r_x]
    ax.plot(r_x, r_y, color='red', label='True regression line')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xlim(0, 1)
    ax.set_xticks([0.2, 0.4, 0.6, 0.8])
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8])
    p_val_str = "{:.5e}".format(p_value)
    ax.set_title(f'Patient: {p} Timepoint: {tp}', fontsize=14)
    
    # ax.text(0.05, 0.95, f"R\u00b2={r_squared:.2f}\nEF={cls}\nP_val={p_val_str}", transform=ax.transAxes, fontsize=12,
    #             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax.text(0.05, 0.95, f"R\u00b2={r_squared:.2f}\nP_val={p_val_str}", transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

def create_plot(EF_df, filterd):
    fig, axes = plt.subplots(8, 6, figsize=(25, 5*8), sharex="none", sharey="none")
    axes_f = axes.flat
    
    i = 0
    for p in PAT_DIRS:
        if p == "N9" and filterd:
            continue
        print(f"\n***Patients {p} ***")
        full_p_dir = res_dir + f"/{p}"
        time_points = [d for d in os.listdir(full_p_dir) if (os.path.isdir(os.path.join(full_p_dir, d)))]
        for tp in time_points:
            print(f"--TP {tp}--")
            print("Creating graph...", end=" ")
            tp_path = full_p_dir + f"/{tp}"
            df_file = fr"{tp_path}/{p}_T{tp}_merged.csv"

            # Read data frame and filter if True
            df = pd.read_csv(df_file)
            if filterd:
                df = df[df['final_freq'] > 0]
            
            # Extract EF value for p/tp
            EF_filtered_df = EF_df[(EF_df['Patient'] == p) & (EF_df['Timepoint'] == int(tp))]
            if not EF_filtered_df.empty:
                cls = EF_filtered_df['EF_Value'].iloc[0]
            else:
                cls = 0

            # Add graph to p/tp
            data_to_graph(df, axes_f[i], p, tp, cls)
            print("Done!")
            i += 1
    # Save fig after all patients and timepoints
    if filterd:
        plt.savefig(fr"Z:/home/volume1/ido/BN-SCRIPTS/graphs/filtered.png")
    else:
        plt.savefig(fr"Z:/home/volume1/ido/BN-SCRIPTS/graphs/not_filtered.png")
    print("Graphs saved, All patients Done!!\n")

    return
    
if __name__ == "__main__":
    try:
        # RESULTS = r'/sternadi/home/volume1/ido/BN-SCRIPTS/Filter_Usecase/results'
        RESULTS = r'Z:/home/volume1/ido/BN-SCRIPTS/Filter_Usecase/results'
        GRAPHS_DIR = r"./graphs"
        res_dir = RESULTS + "/" + results_dir_choice(r'Z:/home/volume1/ido/BN-SCRIPTS/Filter_Usecase/results')
        PAT_DIRS = [d for d in os.listdir(res_dir) if (os.path.isdir(os.path.join(res_dir, d)))]
        EF_df = pd.read_csv(rf"Z:\home\volume1\ido\CalcDepthEF\Data\results.csv")

        print("Creating graph with all data")
        create_plot(EF_df, False)
        print("Creating graph with filtered data")
        create_plot(EF_df, True)        

    except Exception as e:
        print("ERROR!!!")
        print(e)
        exit(1)
