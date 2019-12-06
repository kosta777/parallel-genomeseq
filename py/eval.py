import pandas as pd
import sys
import os
from os.path import dirname, abspath, join as pjoin
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

cwd = dirname(abspath(__file__))
project_dir = pjoin(cwd,"..")

timing1_path = pjoin(project_dir,"data/timings/timing_20191122_1650_ompfg_oxtest.csv")

if __name__ == '__main__':
    if sys.argv[1] == "hello":
        print("hello")

    elif sys.argv[1] == "timing1":

        def poly_fit(x,w0,w1,w2):
            return(w0+w1*x+w2*(x**2))

        df = pd.read_csv(timing1_path)
        df_colnames = list(df.keys())
#        t1_key = df_colnames[3]
        t2_key = df_colnames[4]

#        x_data = np.log10(df['n_threads'].values)
        x_data = df['n_threads'].values

#        t1 = (df[t1_key].values)/(1E6)
#        t2 = (df[t2_key].values)/(1E6)
#        t1 = (df[t1_key].values)/( (df[df['n_threads']==1][t1_key]).mean())
        t2 = (df[t2_key].values)/( (df[df['n_threads']==1][t2_key]).mean())

        w_opt,w_cov = curve_fit(poly_fit,np.log10(x_data),t2)
        x_fit = np.linspace(np.log10(np.min(x_data)),np.log10(np.max(x_data)),1000)
        x_fit_plot = 10**x_fit
        t_fit = poly_fit(x_fit,*w_opt)

        f1 = plt.figure()
#        ax1 = f1.add_subplot(111)
        ax1 = f1.add_subplot(111,xscale="log")

#        ax1.scatter(x_data,t1,s=10.0)
        ax1.scatter(x_data,t2,s=10.0)

        ax1.plot(x_fit_plot,t_fit,linewidth=1.0,color="red")

        ax1.minorticks_on()
        ax1.grid(which='major',linestyle='-',linewidth=0.5)
        ax1.grid(which='minor',linestyle=':',linewidth=0.5)

        ax1.set_title("OMP parallelization of anti-diagonal DP matrix construction (fine grain)")
        ax1.set_xlabel(r"$N_{threads}$")
        ax1.set_xticks(2**(np.arange(0,7)))
        ax1.set_xticklabels(2**(np.arange(0,7)))
        ax1.set_ylabel("Normalized Construction Time") #10k* 30k D

        plt.show()





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


    
        

