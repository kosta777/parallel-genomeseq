import numpy as np
import pandas as pd
import os
from os.path import dirname, abspath, join as pjoin
import sys
import datetime

cwd = dirname(abspath(__file__))

small_fa = pjoin(cwd,"data_small/genome.chr22.5K.fa") #reference genome
small_sam = pjoin(cwd,"data_small/output_tiny_30xCov.mod.sam") #sam result
small_fq_list = [pjoin(cwd,"data_small/output_tiny_30xCov1.fq"), pjoin(cwd,"data_small/output_tiny_30xCov2.fq")] #reads

long_fa = pjoin(cwd,"data/genome.chr22.fa") #reference genome
long_fq_list_5 = [pjoin(cwd,"data/output_5xCov1.fq"), pjoin(cwd,"data/output_5xCov2.fq")] #reads
long_fq_list_10 = [pjoin(cwd,"data/output_10xCov1.fq"), pjoin(cwd,"data/output_10xCov2.fq")] #reads
long_fq_list_30 = [pjoin(cwd,"data/output_30xCov1.fq"), pjoin(cwd,"data/output_30xCov2.fq")] #reads

def general_readtxt(open_file_path):
    f = open(open_file_path,"r")
    fq = f.read()
    f.close()
    return(fq)

class SAM:
    def __init__(self,sam_str):
        self.sam_fields = ['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL']
        self.sam_str = sam_str
        sam = np.array(sam_str.split("\n")[:-1])
        meta_indexes = [i for i in range(len(sam)) if sam[i][0]=="@"]
        self.meta = [sam[i] for i in meta_indexes]
        
        sam = np.delete(sam,meta_indexes)
        self.df = pd.DataFrame(columns=self.sam_fields)
#        print(self.df)
        for i in range(len(sam)):
            row_list = sam[i].split("\t")
            row_dict = {self.sam_fields[j]:[row_list[j]] for j in range(len(row_list))}
#            print(row_dict)
            self.df = self.df.append(pd.DataFrame(row_dict))
        self.df = self.df.reset_index()


def single_fq_2_np(open_file_path):
    f = open(open_file_path,"r")
    fq = f.read()
    f.close()
    fq = fq.split("\n")[:-1]
    
    fq = np.reshape(fq,(int(len(fq)/4),4))
    return(fq) 


def read_fa(open_file_path):
    f = open(open_file_path,"r")
    fa = f.read()
    f.close()
    fa = fa.split("\n")[1:-1]
    fa = "".join(fa)
    return(fa)

def read_multi_fq(fq_file_list):
    fq = []
    for path in fq_file_list:
#        print(path)
        print(np.shape(single_fq_2_np(path)))
        fq.append(single_fq_2_np(path))
    fq = np.array(fq)
    return(fq)

if __name__ == '__main__':
#    import align
#    from cost_and_score_fns import *

    if sys.argv[1] == 'parse_fq':
        fq = single_fq_2_np(small_fq_list[1])
        np.savetxt('small_fq_1.txt',fq[:,1],fmt="%s")

    if sys.argv[1] == 'sam_df':
        sam_raw = general_readtxt(small_sam)
        sam1 = SAM(sam_raw)
        print(sam1.df)

    elif sys.argv[1] == 'data_test':
        print('starting seqdata_test')
        sam_raw = general_readtxt(small_sam)
        sam1 = SAM(sam_raw)

        fa = read_fa(small_fa)
        T = fa

        align_mode = 'smith_waterman'
        print(align_mode)
        result_df = pd.DataFrame()
        for i in range(12,sam1.df.shape[0]):
            timer = []
            Q = sam1.df.iloc[i]['SEQ']
#            print('timer0')
            timer.append(datetime.datetime.now())

            align1 = align.alignment(Q,T,align_mode)
#            print('timer1')
            timer.append(datetime.datetime.now())

#            align1.compute_D(score_smith_waterman)
#            align1.compute_D(score_sw_smallrun)
            align1.compute_D_sw(score_sw_smallrun)
#            print('timer2')
            timer.append(datetime.datetime.now())

            align1.D_get_align()
#            print('timer3')
            timer.append(datetime.datetime.now())
             
            timer = [(timer[i+1]-timer[i]).total_seconds() for i in range(len(timer)-1)]

            result_dict = {
                    "SW_recovered": align1.pos_result['B'],
                    "POS": int(sam1.df.iloc[i]['POS']),
                    "PNEXT": int(sam1.df.iloc[i]['PNEXT']) }
            result_dict = {k:[result_dict[k]] for k in list(result_dict.keys())}
            result_df = result_df.append(pd.DataFrame(result_dict))

            result_str = "timedelta:"+str(timer)+", SW_recovered=%d, POS=%d, PNEXT=%d"%(
                    align1.pos_result['B'],
                    int(sam1.df.iloc[i]['POS']),
                    int(sam1.df.iloc[i]['PNEXT']))
            print("%d, "%(i)+result_str)


    elif sys.argv[1] == 'fq_fq_compr':
        fq = read_multi_fq(small_fq_list)
        print(np.shape(fq))
        sample_seq = fq[0][0][1]
        print(sample_seq)
    
        for i in range(np.shape(fq)[0]):
            for j in range(np.shape(fq)[1]):
                print(j)
                if fq[i][j][1] in [sample_seq,sample_seq[::-1]]:
                    print("%d,%d"%(i,j))
    
    elif sys.argv[1] == 'fq_fa_compr':
        result_arr = [[],[]]
        fq = read_multi_fq(small_fq_list)
        fa = read_fa(small_fa)
        print(np.shape(fq))
    #    sample_seq = fq[0][0][1]
    #    print(sample_seq)
    
        for i in range(np.shape(fq)[0]):
            for j in range(np.shape(fq)[1]):
    #            print(j)
                fq_subseq = fq[i][j][1]
                if fq_subseq in fa:
                    print("%d,%d"%(i,j))
                    loc = (fa.find(fq_subseq))
                    result_arr[i].append([j,loc])
    
        df = [pd.DataFrame(ra) for ra in result_arr]
    
    
    elif sys.argv[1] == 'sam_fa_compr':
        result_arr = []
        sam = single_sam_2_np(small_sam)
        fa = read_fa(small_fa)
        print(np.shape(sam))
    
        for i in range(len(sam)):
            
            if sam[i] in fa:
                print("%d"%(i))
                loc = (fa.find(sam[i]))
                result_arr.append(loc)
    
    #    df = [pd.DataFrame(ra) for ra in result_arr]
        print(result_arr)
