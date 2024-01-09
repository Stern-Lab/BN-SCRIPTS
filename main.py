from datetime import datetime
import time
import Filter_Usecase.filter_usecase as filterpy
import BN_Create_input.create_prep_files as createprepfile
import Run_BN_script.run_bn_script as bn

if __name__ == "__main__":
    try:
        start = time.time()
        print("Main script script is starting...")
        print(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
        print("What do you want to do?")
        print("1-Run filter and create prepartion files scripts")
        print("2-Run only filter script")
        print("3-Run only prepartion files scripts")
        print("4-Run Bottleneck estimation R script")
        while True:
            user_input = input("Choose: ")
            if user_input == "1":
                filterpy.main()
                createprepfile.main()
                break
            elif user_input == "2":
                filterpy.main()
                break
            elif user_input == "3":
                createprepfile.main()
                break
            elif user_input == "4":
                bn.user_main()
                break
            else:
                print("Wrong input. Please enter 1, 2, 3 or 4!")
            
    except Exception as e:
        print("An error has occured!\nTerminating script...")
        print(f"Main script elapsed time: {(time.time() - start)} sec")
        print("Log:")
        print(e)
        exit(1)

    print("====Main Script finished successfully!====")
    print(f"Main Script elapsed time: {(time.time() - start)} sec")
