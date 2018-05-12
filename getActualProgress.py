#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 12 14:46:55 2018

@author: michal
"""

import glob
import time, datetime
    
cifProcessed = 0
log_files = glob.glob("logs/partialProgress*")
for log_file in log_files:
    log = open(log_file, 'r')
    cifProcessed += int(log.readline())
    log.close()

cifNoFile = open("logs/cif2process.log", 'r')
cifNo = int(cifNoFile.readline())
cifNoFile.close()

timeFile = open("logs/timeStart.log", 'r')
timeStart = float(timeFile.readline())
timeFile.close()

timeActual = time.time()
timeTaken = timeActual - timeStart

timeEstimated = timeTaken/cifProcessed * (cifNo  - cifProcessed)
prettyTimeTaken = str(datetime.timedelta(seconds = timeTaken))
prettyTimeEstimated = str( datetime.timedelta(seconds = timeEstimated) )

progress = float(cifProcessed)/cifNo * 100
print("##########################################")
print("#################POSTEP###################")
print("##########################################")
print("Przetworzono: ", cifProcessed, "/", cifNo)
print(progress, "%")
print("W czasie: ", prettyTimeTaken)
print("Szacowany pozostały czas: ", prettyTimeEstimated)