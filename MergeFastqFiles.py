import time
import os
# from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform

import Util
import Logic
import LogicPrep
#################### st env ####################
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
SYSTEM_NM = platform.system()

if SYSTEM_NM == 'Linux':
    # REAL
    pass
else:
    # DEV
    REF_DIR = "G:/"

FASTQ = 'FASTQ/'
IN = 'input/'
OU = 'output/'

MIS_MATCH_CUT_OFF = 1

MERGE_FILE_LIST = 'merge_all.txt'

os.makedirs(WORK_DIR + IN, exist_ok=True)
os.makedirs(WORK_DIR + OU, exist_ok=True)

# do not use multi process for editing files
# TOTAL_CPU = mp.cpu_count()
# MULTI_CNT = int(TOTAL_CPU*0.8)
#################### en env ####################


def merge_fastq_files(file_list):
    for file_path_arr in file_list:
        result_file_name = 'all_' + '_'.join(file_path_arr[0].split('/')[-1].split('_')[:-1]) + '.extendedFrags.fastq'
        print('merge :', result_file_name)
        # result_path = '/'.join(file_path_arr[0].split('/')[:-1]) + '/'
        # with open(REF_DIR + result_path + result_file_name, 'w') as f_out:
        # with open('./' + OU + result_file_name, 'w') as f_out:

        # # output path must be different drive not input drive
        with open('C:/Users/terry/Desktop/' + OU + result_file_name, 'w') as f_out:
            for file_path in file_path_arr:
                f_in = open(REF_DIR + file_path)
                print('read :', file_path)
                for i, line in enumerate(f_in):
                    f_out.write(line)
                f_in.close()
        print('done :', result_file_name)


def non_multi_merge_files():
    util = Util.Utils()

    file_list = util.read_tsv_ignore_N_line('./' + IN + MERGE_FILE_LIST, 0)
    merge_fastq_files(file_list)

# do not use multi process for editing files
"""
def mutil_merge_files(multi_cnt):
    util = Util.Utils()

    file_list = util.read_tsv_ignore_N_line('./' + IN + MERGE_FILE_LIST, 0)
    len_file_list = len(file_list)

    if multi_cnt > len_file_list:
        multi_cnt = len_file_list

    # divide data_list by multi_cnt
    splited_file_list = np.array_split(file_list, multi_cnt)
    file_list.clear()
    print("platform.system() : ", SYSTEM_NM)
    print("total cpu_count : ", str(TOTAL_CPU))
    print("will use : ", str(multi_cnt))
    pool = mp.Pool(processes=multi_cnt)

    pool.map(merge_fastq_files, splited_file_list)
    pool.close()
"""


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # mutil_merge_files(MULTI_CNT)
    non_multi_merge_files()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))