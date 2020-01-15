import pandas as pd
import sys
import os
from os.path import dirname, abspath, join as pjoin
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import argparse

thisfile_dir = dirname(abspath(__file__))
project_dir = pjoin(thisfile_dir,"..")

timing1_path = pjoin(project_dir,"data/timings/ompfg_timing_results.csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--option',default='hello') #ompfg
    parser.add_argument('-y','--yaxis',default='abs_time') #abs_time, normed_time
    parser.add_argument('-p','--plot_type',default='box_plot') #box_plot, scatter
    parser.add_argument('-f','--fit',default='false') #true
    args = parser.parse_args()
    argsdict = vars(args)
    print(argsdict)

    if argsdict['option'] == "hello":
        print("hello")

    elif argsdict['option'] == "ompfg":

        def poly_fit(x,w0,w1,w2):
            return(w0+w1*x+w2*(x**2))

        df = pd.read_csv(timing1_path)
        df_colnames = list(df.keys())
        t2_key = df_colnames[4]

        x_data = df['n_threads'].values

        f1 = plt.figure()
#        ax1 = f1.add_subplot(111,xscale="log")
        ax1 = f1.add_subplot(111)

        if argsdict['yaxis'] == "abs_time":
            df['y'] = (df[t2_key].values)/(1E6)
            ylabel_txt = "Abs Construction Time (s)"
        elif argsdict['yaxis'] == "normed_time":
            seq_mean = (df[df['n_threads']==1][t2_key]).mean()
            df['y'] = (df[t2_key].values)/seq_mean
            ylabel_txt = "Normalized Construction Time"
        elif argsdict['yaxis'] == "speedup":
            seq_mean = (df[df['n_threads']==1][t2_key]).mean()
            df['y'] = seq_mean/(df[t2_key].values)
        elif argsdict['yaxis'] == "gcups":
            gc = (1E4)*(3E4)/(1E9) #fixed
            df['y'] = gc/(df[t2_key].values/1E6)
            df['y_reci'] = 1/df['y']
#            ylabel_txt = "GCUPS (nb: gc=0.3 fixed)"
            ylabel_txt = "GCUPS"


        if argsdict['plot_type'] == "scatter":
            x_data = np.log2(x_data)
            ax1.scatter(x_data,df['y'].values,s=10.0)
        elif argsdict['plot_type'] == "box_plot":
            x_data_unique = np.unique(df['n_threads'].values)
            bp_data = [df[df['n_threads']==x]['y'].values for x in x_data_unique]

#            flier_settings = dict(markersize=3.5, markerfacecolor='g', marker='D' )
#            flier_settings = dict(markersize=3.5, markerfacecolor='black', marker='o' )
#            ax1.boxplot(x=bp_data,positions=np.log2(x_data_unique),widths=0.15,flierprops=flier_settings)
            ax1.boxplot(boxprops=dict(linewidth=1.5),x=bp_data,positions=np.log2(x_data_unique),widths=0.15,showfliers=False)

        if argsdict['fit'] == 'poly':
            w_opt,w_cov = curve_fit(poly_fit,np.log2(x_data),df['y'].values)
            x_fit = np.linspace(np.log2(np.min(x_data)),np.log2(np.max(x_data)),1000)
            t_fit = poly_fit(x_fit,*w_opt)
            ax1.plot(x_fit,t_fit,linewidth=1.0,color="red")
        elif argsdict['fit'] == 'hmean':
            y_hmean = np.zeros(len(x_data_unique)) #harmonic mean
            for i in range(len(x_data_unique)):
                y_hmean[i] = 1/(df[df['n_threads']==x_data_unique[i]]['y_reci'].mean())
            ax1.plot(np.log2(x_data_unique),y_hmean,linewidth=1.0,color="red",label="Harmonic mean")
            ax1.legend(loc='upper left',fontsize=12)
            ax1.scatter(np.log2(df['n_threads'].values),df['y'],s=5.0,color='black',marker='o')

        ax1.minorticks_on()
        ax1.grid(which='major',linestyle='-',linewidth=0.5)
        ax1.grid(which='minor',linestyle=':',linewidth=0.5)

#        ax1.set_title("OMP parallelization of anti-diagonal DP matrix construction (fine grain)")
        ax1.set_xlabel(r"$\mathregular{N_{threads}}$", fontsize=16)
#        ax1.set_xticks(2**(np.arange(0,7)))
        ax1.set_xticks(np.arange(0,7))
        ax1.set_xticklabels(2**(np.arange(0,7)))

        ax1.set_ylabel(ylabel_txt, fontsize=16) #10k* 30k D

        ax1.set_aspect(aspect=(3/4)/ax1.get_data_ratio())
        plt.show()


    elif argsdict['option'] == "sw_solve_small":
        parser.add_argument('-aln','--align_file',default="data/align_output.csv")
        args = parser.parse_args()
        argsdict = vars(args)
        print(argsdict)

        align_output_file = pjoin(project_dir,argsdict['align_file'])

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


    
        

