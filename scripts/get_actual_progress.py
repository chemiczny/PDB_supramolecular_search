#!/usr/bin/env python3
"""
Created on Sat May 12 14:46:55 2018

@author: michal
"""

import glob
import time
import datetime
from os.path import isfile

cif_processed = 0
log_files = glob.glob("logs/partialProgress*")
for log_file in log_files:
    log = open(log_file, "r")
    cif_processed += int(log.readline())
    log.close()

cif_no_file = open("logs/cif2process.log", "r")
cif_no = int(cif_no_file.readline())
cif_no_file.close()

time_file = open("logs/timeStart.log", "r")
time_start = float(time_file.readline())
time_file.close()

time_stop = -1
if isfile("logs/timeStop.log"):
    time_file = open("logs/timeStop.log", "r")
    time_stop = float(time_file.readline())
    time_file.close()

time_actual = time.time()

progress = float(cif_processed) / cif_no * 100
if abs(progress - 100) < 0.00001 or time_stop > time_start:
    time_actual = time_stop

time_taken = time_actual - time_start

time_estimated = time_taken / cif_processed * (cif_no - cif_processed)
pretty_time_taken = str(datetime.timedelta(seconds=time_taken))
pretty_time_estimated = str(datetime.timedelta(seconds=time_estimated))

print("##########################################")
print("################PROGRESS##################")
print("##########################################")
print("Processed: ", cif_processed, "/", cif_no)
print(progress, "%")
print("Time taken: ", pretty_time_taken)
if time_stop < time_start:
    print("Estimated time remaining: ", pretty_time_estimated)
