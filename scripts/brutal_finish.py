#!/usr/bin/env python3
"""
Created on Mon Oct  1 14:29:26 2018

@author: michal
"""

import glob
from os.path import join
from os import remove
import time
from supramolecular_search.config import configure
from multiprocessing import Pool


def merge_logs(log_list):
    log_final = log_list[0]
    logs = log_list[1]

    log_files = glob.glob(logs)
    for log_file in log_files:
        if log_file == log_final:
            continue

        new_log = open(log_file, "r")
        final_log = open(log_final, "a+")
        line = new_log.readline()
        while line:
            final_log.write(line)
            line = new_log.readline()

        new_log.close()
        final_log.close()
        remove(log_file)


def merge_progress_files(scratch):
    cif_processed = 0
    log_files = glob.glob(join(scratch, "partialProgress*"))
    for log_file in log_files:
        log = open(log_file, "r")
        cif_processed += int(log.readline())
        log.close()
        remove(log_file)

    cif_no_file = open("logs/cif2process.log", "r")
    cif_no = int(cif_no_file.readline())
    cif_no_file.close()

    progress_summary = open("logs/progressSummary.log", "w")
    progress_summary.write(
        "Przetworzono: " + str(cif_processed) + "/" + str(cif_no) + "\n"
    )
    progress_summary.close()


config = configure()
number_of_processes = config["N"]
scratch = config["scratch"]

#
pool = Pool(number_of_processes)

arguments_list = []
arguments_list.append(["logs/anionPi.log", join(scratch, "anionPi*.log")])
arguments_list.append(["logs/cationPi.log", join(scratch, "cationPi*.log")])
arguments_list.append(["logs/piPi.log", join(scratch, "piPi*.log")])
arguments_list.append(["logs/anionCation.log", join(scratch, "anionCation*.log")])
arguments_list.append(["logs/hBonds.log", join(scratch, "hBonds*.log")])
arguments_list.append(["logs/metalLigand.log", join(scratch, "metalLigand*.log")])
arguments_list.append(["logs/planarAnionPi.log", join(scratch, "planarAnionPi*.log")])
arguments_list.append(["logs/linearAnionPi.log", join(scratch, "linearAnionPi*.log")])
arguments_list.append(["logs/methylPi.log", join(scratch, "methylPi*.log")])
arguments_list.append(["logs/additionalInfo.log", join(scratch, "additionalInfo*.log")])

pool.map(merge_logs, arguments_list)

merge_progress_files(scratch)

time_stop = time.time()
time_file = open("logs/timeStop.log", "w")
time_file.write(str(time_stop))
time_file.close()
