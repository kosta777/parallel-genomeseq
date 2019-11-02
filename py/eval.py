import pandas as pd
import sys

if __name__ == '__main__':
    if sys.argv[1] == "hello":
        print("hello")

    elif sys.argv[1] == "sw_solve_small":
        if len(sys.argv)<3:
            align_output_file = "../data/align_output.csv"
        else:
            align_output_file = sys.argv[2]

        df_ao = pd.read_csv(align_output_file)
        df_ao = df_ao.set_index('index')
        df_ao['delta_pos'] = df_ao['pos_pred'] - df_ao['POS']
        df_report = df_ao[df_ao['delta_pos']!=0]
        diff_size = df_report.shape[0]
        total_size = df_ao.shape[0]
        if diff_size > 0:
            print("%d/%d alignments different from ground truth"%(diff_size,total_size))
            print("May be caused by cost function. There is often no unique correct solution.")
            print(df_report)
        else:
            print("No diffs")


    
        

