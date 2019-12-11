import argparse
from os.path import dirname, abspath, join as pjoin

import pandas as pd
from reader import *

thisfile_dir = dirname(abspath(__file__))
project_dir = pjoin(thisfile_dir, "..")

long_fa = pjoin(project_dir, "data/chr22.fa")  # reference genome

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--option', default='hello')

    parser.add_argument('-sp', '--start_pos', default=300000 * 60)
    parser.add_argument('-rfl', '--ref_len', default=30 * 1000)
    parser.add_argument('-o1', '--output_ref_path', default="data/custom_ref_1.fa")
    parser.add_argument('-rmN', '--remove_N', default='true')

    parser.add_argument('-rdl', '--read_len', default=10000)
    parser.add_argument('-nr', '--n_reads', default=100)
    parser.add_argument('-in', '--input_ref_path', default="data/custom_ref_1.fa")
    parser.add_argument('-o2', '--output_reads_path', default="data/custom_reads_1.csv")

    args = parser.parse_args()
    argsdict = vars(args)
    print(argsdict)

    if argsdict['option'] == 'gen_ref_custom':
        #        test_start_pos = 300000 * 60
        #        custom_ref_len = 30 * 1000
        #        output_ref_path = "data/test_fa_custom1.fa"

        test_start_pos = int(argsdict['start_pos'])
        custom_ref_len = int(argsdict['ref_len'])
        output_ref_path = pjoin(project_dir, argsdict['output_ref_path'])

        alphabet = ['A', 'C', 'G', 'T', 'N']

        new_fa_file = open(output_ref_path, "w")
        line_idx = 0
        nt_idx = 0

        custom_ref_str = ""

        with open(long_fa, "r") as original_fa:
            for fa_line in original_fa:
                nts_in_line = 0
                if line_idx > 0:
                    line_str = (fa_line.split("\n")[0]).upper()
                    nts_in_line = len(line_str)
                    if nt_idx >= test_start_pos and nt_idx < test_start_pos + custom_ref_len:
                        custom_ref_str += line_str
                    if line_idx % 10000 == 0 or line_idx < 10:
                        print("line_idx=%d: " % (line_idx) + line_str)
                #                    print(fa_line.split("\n")[1:-1])
                line_idx += 1
                nt_idx += nts_in_line
        if argsdict['remove_N'] == 'true':
            custom_ref_str = custom_ref_str.replace('N', '')
        new_fa_file.write(custom_ref_str)
        new_fa_file.close()

        print("custom ref generated, ACGT counts: ")
        for c in alphabet:
            print("%s:%d" % (c, custom_ref_str.count(c)))

    elif argsdict['option'] == 'gen_reads_custom':
        import random

        #        custom_read_len = 10000
        #        n_reads_to_gen = 100
        #        input_ref_path = "data/test_fa_custom1.fa"
        #        output_reads_path = "data/ground_truth_custom10k.csv"

        custom_read_len = int(argsdict['read_len'])
        n_reads_to_gen = int(argsdict['n_reads'])
        input_ref_path = argsdict['input_ref_path']
        output_reads_path = argsdict['output_reads_path']
        output_readsonly_path = output_reads_path.split(".")[0] + "_readsonly.txt"

        f = open(pjoin(project_dir, input_ref_path), "r")
        custom_ref = f.read()
        f.close()

        print(custom_ref)
        print(len(custom_ref))
        custom_ref_len = len(custom_ref)

        gt_df = pd.DataFrame({
            "QNAME": [],
            "SEQ": [],
            "POS": [],
        })

        for i in range(n_reads_to_gen):
            start_idx = int(random.random() * (custom_ref_len - custom_read_len))
            gt_df = gt_df.append(
                pd.DataFrame({
                    "QNAME": ["custom_read_%d" % (i)],
                    "SEQ": [custom_ref[start_idx:start_idx + custom_read_len]],
                    "POS": [str(start_idx)],
                })
            )

        gt_df = gt_df[["QNAME", "SEQ", "POS"]]
        gt_df = gt_df.reset_index(drop=True)
        gt_df.index.name = "index"
        print(gt_df.head())
        print(gt_df.shape)
        print(gt_df.keys())
        #        gt_df = gt_df.reset_index()
        open(pjoin(project_dir, output_reads_path), "w").close()
        gt_df.to_csv(pjoin(project_dir, output_reads_path))

        open(pjoin(project_dir, output_readsonly_path), "w").close()
        reads_only_file = open(pjoin(project_dir, output_readsonly_path), "a")
        for i in range(gt_df.shape[0]):
            reads_only_file.write(gt_df.iloc[i]['SEQ'] + '\n')
        reads_only_file.close()
