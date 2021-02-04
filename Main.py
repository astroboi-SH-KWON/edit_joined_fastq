import time
import os
# from Bio import SeqIO
import multiprocessing as mp
# import numpy as np
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
    WORK_DIR = "D:/000_WORK/JangHyeWon_LeeMinYung/20200703/WORK_DIR/"
    REF_DIR = "D:/000_WORK/000_reference_path/"

FASTQ = 'FASTQ/'
IN = 'input/'
OU = 'output/'

MIS_MATCH_CUT_OFF = 1

BRCD_FILE = "barcode_list.txt"
BRCD1_POS =[0, 9]
LEN_BRCD = 9
GAP_ARR = [0, 1, 2, 3]
CONST1 = "AGTACGTACGAGTC"  # 14 bp
CONST2 = "GTACTCGCAGTAGTC"  # 15 bp
# CONST2 = 'TTTTT'  # test.fastq
CONST2_POS = [30, 37]
# CONST2_POS = [70, 90]  # test.fastq

os.makedirs(WORK_DIR + IN, exist_ok=True)
os.makedirs(WORK_DIR + OU, exist_ok=True)

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
#################### en env ####################


def split_FASTQ_by_1500x1500_cell_id():
    # trgt_fastq = '253_2_S7_L001_R1_001.fastq'
    trgt_fastq = 'test_253_2_S7_L001_R1_001.fastq'

    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    brcd_list = util.read_tsv_ignore_N_line(WORK_DIR + IN + BRCD_FILE)
    brcd_arr = [tmp_arr[0] for tmp_arr in brcd_list]

    result_dict = {}
    with open(WORK_DIR + FASTQ + trgt_fastq) as f_in:
        while True:
            tmp_id = f_in.readline()
            if tmp_id.replace(' ', '') == '':
                break
            ngs_read = f_in.readline()
            thrd_lin = f_in.readline()
            qual_lin = f_in.readline()

            brcd1_fr_ngs = ngs_read[BRCD1_POS[0]: BRCD1_POS[1]]
            if brcd1_fr_ngs not in brcd_arr:
                continue

            cons2_pos = ngs_read.find(CONST2, BRCD1_POS[1])
            len_cons2 = len(CONST2)

            if CONST2_POS[0] <= cons2_pos <= CONST2_POS[1]:
                brcd2_fr_ngs = ngs_read[cons2_pos - LEN_BRCD: cons2_pos]
                if brcd2_fr_ngs not in brcd_arr:
                    continue

                brcd_paired_key = brcd1_fr_ngs + "_" + brcd2_fr_ngs
                if brcd_paired_key in result_dict:
                    result_dict[brcd_paired_key] += (tmp_id + ngs_read + thrd_lin + qual_lin)
                else:
                    result_dict.update({brcd_paired_key: tmp_id + ngs_read + thrd_lin + qual_lin})

            else:  # get n_mismatch const2 by MIS_MATCH_CUT_OFF
                for i in range(CONST2_POS[0], CONST2_POS[1]):
                    const_candi = ngs_read[i: i + len_cons2].upper()

                    # check length
                    if len(const_candi) != len_cons2:
                        continue

                    cnt = 0
                    for idx in range(len_cons2):
                        if const_candi[idx] != CONST2[idx]:
                            cnt += 1
                        if cnt > MIS_MATCH_CUT_OFF:
                            break
                    if cnt <= MIS_MATCH_CUT_OFF:
                        brcd2_fr_ngs = ngs_read[i - LEN_BRCD: i]
                        if brcd2_fr_ngs not in brcd_arr:
                            continue

                        brcd_paired_key = brcd1_fr_ngs + "_" + brcd2_fr_ngs
                        if brcd_paired_key in result_dict:
                            result_dict[brcd_paired_key] += (tmp_id + ngs_read + thrd_lin + qual_lin)
                        else:
                            result_dict.update({brcd_paired_key: tmp_id + ngs_read + thrd_lin + qual_lin})

    fastq_dir = trgt_fastq.replace(".fastq", "") + "/"
    os.makedirs(WORK_DIR + FASTQ + fastq_dir, exist_ok=True)
    for brcd_key, tmp_str in result_dict.items():
        with open(WORK_DIR + FASTQ + fastq_dir + brcd_key + '_' + trgt_fastq, 'w') as f_ou:
            f_ou.write(tmp_str)


def copy_n_paste_const2_bfr_trgt_fastq():
    trgt_fastq = '253_2_S7_L001_R1_001.fastq'
    with open(WORK_DIR + FASTQ + trgt_fastq) as f_in:
        with open(WORK_DIR + FASTQ + 'const2_bfr_' + trgt_fastq, 'w') as f_ou:
            while True:
                tmp_id = f_in.readline()
                if tmp_id.replace(' ', '') == '':
                    break
                ngs_read = f_in.readline()
                thrd_lin = f_in.readline()
                qual_lin = f_in.readline()

                cons2_pos = ngs_read.find(CONST2, BRCD1_POS[1])
                len_cons2 = len(CONST2)

                if CONST2_POS[0] <= cons2_pos <= CONST2_POS[1]:
                    new_ngs_read = ngs_read[cons2_pos: cons2_pos + len_cons2] + ngs_read
                    new_qual_lin = qual_lin[cons2_pos: cons2_pos + len_cons2] + qual_lin

                    f_ou.write(tmp_id)
                    f_ou.write(new_ngs_read)
                    f_ou.write(thrd_lin)
                    f_ou.write(new_qual_lin)

                else:  # get n_mismatch const2 by MIS_MATCH_CUT_OFF
                    for i in range(CONST2_POS[0], CONST2_POS[1]):
                        const_candi = ngs_read[i: i + len_cons2].upper()

                        # check length
                        if len(const_candi) != len_cons2:
                            continue

                        cnt = 0
                        for idx in range(len_cons2):
                            if const_candi[idx] != CONST2[idx]:
                                cnt += 1
                            if cnt > MIS_MATCH_CUT_OFF:
                                break
                        if cnt <= MIS_MATCH_CUT_OFF:
                            quali_candi = qual_lin[i: i + len_cons2]
                            f_ou.write(tmp_id)
                            f_ou.write(const_candi + ngs_read)
                            f_ou.write(thrd_lin)
                            f_ou.write(quali_candi + qual_lin)


def switch_brcd1_aftr_const2_in_trgt_fastq():
    trgt_fastq = '253_2_S7_L001_R1_001.fastq'
    # trgt_fastq = 'test.fastq'  # test.fastq
    with open(WORK_DIR + FASTQ + trgt_fastq) as f_in:
        with open(WORK_DIR + FASTQ + 'brcd1_aftr_const2_' + trgt_fastq, 'w') as f_ou:
            while True:
                tmp_id = f_in.readline()
                if tmp_id.replace(' ', '') == '':
                    break
                ngs_read = f_in.readline()
                thrd_lin = f_in.readline()
                qual_lin = f_in.readline()

                brcd1 = ngs_read[BRCD1_POS[0]: BRCD1_POS[1]]
                cons2_pos = ngs_read.find(CONST2, BRCD1_POS[1])
                len_cons2 = len(CONST2)

                if CONST2_POS[0] <= cons2_pos <= CONST2_POS[1]:
                    new_ngs_read = ngs_read[BRCD1_POS[1]: cons2_pos + len_cons2] + brcd1 + ngs_read[
                                                                                           cons2_pos + len_cons2:]
                    new_qual_lin = qual_lin[BRCD1_POS[1]: cons2_pos + len_cons2] + qual_lin[BRCD1_POS[0]: BRCD1_POS[
                        1]] + qual_lin[cons2_pos + len_cons2:]

                    f_ou.write(tmp_id)
                    f_ou.write(new_ngs_read)
                    f_ou.write(thrd_lin)
                    # f_ou.write(str(cons2_pos) + '\n')  # test.fastq
                    f_ou.write(new_qual_lin)


def add_given_seq_befr_trgt_fastq():
    given_seq = "A" * 15
    given_qual = "F" * 15

    trgt_fastq = '253_2_S7_L001_R1_001.fastq'
    with open(WORK_DIR + FASTQ + trgt_fastq) as f_in:
        with open(WORK_DIR + FASTQ + given_seq + '_bfr_' + trgt_fastq, 'w') as f_ou:
            while True:
                tmp_id = f_in.readline()
                if tmp_id.replace(' ', '') == '':
                    break
                ngs_read = f_in.readline()
                thrd_lin = f_in.readline()
                qual_lin = f_in.readline()

                f_ou.write(tmp_id)
                f_ou.write(given_seq + ngs_read)
                f_ou.write(thrd_lin)
                f_ou.write(given_qual + qual_lin)


tmp_list = [
    '@A00547:139:HLCM3DSXY:1:1101:20274:1282 1:N:0:CGCATGAT+AAGCCTGA'
    , ''
    , ''
    , ''
    , ''
    , ''
]


def test():
    # trgt_fastq = 'const2_bfr_253_2_S7_L001_R1_001.fastq'
    trgt_fastq = '253_2_S7_L001_R1_001.fastq'
    with open(WORK_DIR + FASTQ + trgt_fastq) as f_in:
        with open(WORK_DIR + FASTQ + 'test_' + trgt_fastq, 'w') as f_ou:
            cnt = 0
            while True:
                tmp_id = f_in.readline()
                if tmp_id.replace(' ', '') == '':
                    break

                if cnt >= 10000:
                    break
                ngs_read = f_in.readline()
                thrd_lin = f_in.readline()
                qual_lin = f_in.readline()

                f_ou.write(tmp_id)
                f_ou.write(ngs_read)
                f_ou.write(thrd_lin)
                f_ou.write(qual_lin)
                cnt += 1


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # test()
    # split_FASTQ_by_1500x1500_cell_id()
    # add_given_seq_befr_trgt_fastq()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))