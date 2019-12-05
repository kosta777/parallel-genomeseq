import numpy as np
import pandas as pd
import os
from os.path import dirname, abspath, join as pjoin
import sys
import datetime
import pickle

cwd = dirname(abspath(__file__))
project_dir = pjoin(cwd,"..")

small_fa = pjoin(cwd,"data_small/genome.chr22.5K.fa") #reference genome
small_sam = pjoin(cwd,"../data/data_small/output_tiny_30xCov.mod.sam") #sam result
small_fq_list = [pjoin(cwd,"data_small/output_tiny_30xCov1.fq"), pjoin(cwd,"data_small/output_tiny_30xCov2.fq")] #reads

long_fa = pjoin(project_dir,"data/genome.chr22.fa") #reference genome
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

def mpi_prepare(text_str):
    text = text_str.split("\n")
    f = open("../data/data_small/mpi_test_tiny.fq", "w")
    for i in range(1, len(text), 4):
        f.write(text[i]+'\n')
    f.close()

def uniprot_prepare(file_path): 
    token = '>sp'
    chunks = []
    current_chunk = []

    for line in open("../data/uniprot/uniprot_sprot.fasta"):
       if line.startswith(token) and current_chunk: 
          # if line starts with token and the current chunk is not empty
          chunks.append(current_chunk[:]) #  add not empty chunk to chunks
          current_chunk = [] #  make current chunk blank
       # just append a line to the current chunk on each iteration
       current_chunk.append(line)

    chunks.append(current_chunk)  #  append the last chunk outside the loop
    i=0
    for chunk in chunks:
        with open("../data/uniprot/"+str(i)+".fasta", 'w') as f:
            for el in chunk:
                f.write("%s" % el)
        i=i+1
    with open("../data/uniprot/stats.txt", 'w+') as f:
        f.write("%d" % i)

def uniprot_prepare_single(file_path): 
    token = '>sp'
    chunks = []
    current_chunk = ""

    for line in open("../data/uniprot/uniprot_sprot.fasta"):
       if line.startswith(token) and current_chunk != "": 
          # if line starts with token and the current chunk is not empty
          chunks.append(current_chunk[:]) #  add not empty chunk to chunks
          current_chunk = "" #  make current chunk blank
       # just append a line to the current chunk on each iteration
       if not line.startswith(token):
          current_chunk = current_chunk + line[:-1]

    chunks.append(current_chunk)  #  append the last chunk outside the loop
    i=0
    with open("../data/uniprot/database.fasta", 'a+') as f:
        for chunk in chunks:
            f.write("%s\n" % chunk)
        i=i+1
    with open("../data/uniprot/stats.txt", 'w+') as f:
        f.write("%d" % i)

 

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

    elif sys.argv[1] == 'sam_df':
        sam_raw = general_readtxt(small_sam)
        sam1 = SAM(sam_raw)
        print(sam1.df)

    elif sys.argv[1] == 'create_single_string_fa':
        parsed_fa_string = read_fa(small_fa)
        new_fa_file = open("data/test_fa.fa","w")
#        new_fa_file = open("data/test_fa.fa","a")
        new_fa_file.write(parsed_fa_string)
        new_fa_file.close()

    elif sys.argv[1] == 'create_single_string_fa_long':
        open("data/test_fa_long.fa", "w").close()
        new_fa_file = open("data/test_fa_long.fa","a")

        i=0
#        with open(small_fa,"r") as original_fa:
        with open(long_fa,"r") as original_fa:
            for fa_line in original_fa:
                if i>0:
                    new_fa_file.write(fa_line.split("\n")[0])
                    if i%1000==0 or i<10:
                        print("i=%d: "%(i)+fa_line.split("\n")[0])
    #                    print(fa_line.split("\n")[1:-1])
                i+=1
        new_fa_file.close()


    elif sys.argv[1] == 'gen_ref_custom':
        custom_ref_len = 30 * 1000
        test_start_pos = 300000 * 60
        alphabet = ['A','C','G','T','N']

        open("data/test_fa_custom1.fa", "w").close()
        new_fa_file = open("data/test_fa_custom1.fa","a")
        line_idx=0
        nt_idx = 0

        custom_ref_str = ""

        with open(long_fa,"r") as original_fa:
            for fa_line in original_fa:
                nts_in_line = 0
                if line_idx>0:
                    line_str=fa_line.split("\n")[0]
                    nts_in_line = len(line_str)
                    if nt_idx >= test_start_pos and nt_idx < test_start_pos+custom_ref_len:
                        custom_ref_str+=line_str
                        new_fa_file.write(line_str)
                    if line_idx%1000==0 or line_idx<10:
                        print("line_idx=%d: "%(line_idx)+line_str)
    #                    print(fa_line.split("\n")[1:-1])
                line_idx+=1
                nt_idx += nts_in_line
        new_fa_file.close()

        print("results: ")
        for c in alphabet:
            print("%s:%d"%(c,custom_ref_str.count(c)))


    elif sys.argv[1] == 'gen_input_custom':
        import random

        custom_read_len = 10000
        n_reads_to_gen = 100

        f = open(  pjoin(project_dir,"data/test_fa_custom1.fa"),"r")
        custom_ref = f.read()
        f.close()

        print(custom_ref)
        print(len(custom_ref))
        custom_ref_len = len(custom_ref)

        gt_df = pd.DataFrame({
            "QNAME":[],
            "SEQ":[],
            "POS":[],
            })


        for i in range(n_reads_to_gen):
            start_idx = int(random.random() * (custom_ref_len-custom_read_len))
            gt_df = gt_df.append(
                        pd.DataFrame({
                        "QNAME":["custom_read_%d"%(i)],
                        "SEQ":[custom_ref[start_idx:start_idx+custom_read_len]],
                        "POS":[str(start_idx)],
                        })
                    )

        gt_df = gt_df[["QNAME","SEQ","POS"]]
        gt_df = gt_df.reset_index(drop=True)
        gt_df.index.name = "index"
        print(gt_df.head())
        print(gt_df.keys())
#        gt_df = gt_df.reset_index()
        open(pjoin(project_dir,"data/ground_truth_custom10k.csv"),"w").close()
        gt_df.to_csv( pjoin(project_dir, "data/ground_truth_custom10k.csv"))

        open(pjoin(project_dir,"data/ground_truth_custom10k_readsonly.txt"),"w").close()
        reads_only_file = open( pjoin(project_dir,"data/ground_truth_custom10k_readsonly.txt"),"a")
        for i in range(gt_df.shape[0]):
            reads_only_file.write(gt_df.iloc[i]['SEQ']+'\n')
        reads_only_file.close()


    elif sys.argv[1] == 'gen_input_125':
        sam_raw = general_readtxt(small_sam)
        sam1 = SAM(sam_raw)

        gt_df = pd.DataFrame({
            "QNAME":sam1.df["QNAME"],
            "SEQ":sam1.df["SEQ"],
            "POS":sam1.df["POS"],
            })

        gt_df.index.name = "index"

        gt_df.to_csv("../data/ground_truth.csv")

    elif sys.argv[1] == 'read_fa':
        fa = read_fa(small_fa)



    elif sys.argv[1] == 'create_lsv_reads_from_sam':  #lsv="line separated values"
        sam_raw = general_readtxt(small_sam)
        sam1 = SAM(sam_raw)
        print(sam1.df)

        new_reads_file = open("data/test_reads.txt","a")
        for i in range(sam1.df.shape[0]):
            new_reads_file.write(sam1.df.iloc[i]['SEQ']+'\n')
        new_reads_file.close()

    elif sys.argv[1] == 'mpi_prepare_input':
        mpi_prepare(general_readtxt("../data/data_small/output_tiny_30xCov1.fq"))
    elif sys.argv[1] == 'uniprot_prepare':
        uniprot_prepare("")


###############



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
