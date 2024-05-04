import os

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

if __name__ == "__main__":
    GRAPHS_DIR = r"./graphs"
    
    return