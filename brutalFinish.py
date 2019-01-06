#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:29:26 2018

@author: michal
"""
import glob
from os import remove
import time

def mergeLogs( logFinal, logs  ):
    final_log = open(logFinal, "a+")
    log_files = glob.glob(logs)
    for log_file in log_files:
        if log_file == logFinal:
            continue
        
        new_log = open(log_file, 'r')
        
        line = new_log.readline()
        while line:
            final_log.write(line)
            line = new_log.readline()
        
        new_log.close()
        remove(log_file)
    
    final_log.close()
    
def mergeProgressFiles():
    cifProcessed = 0
    log_files = glob.glob("logs/partialProgress*")
    for log_file in log_files:
        log = open(log_file, 'r')
        cifProcessed += int(log.readline())
        log.close()
        remove(log_file)
    
    cifNoFile = open("logs/cif2process.log", 'r')
    cifNo = int(cifNoFile.readline())
    cifNoFile.close()
    
    progressSummary = open("logs/progressSummary.log", 'w')
    progressSummary.write("Przetworzono: " +str( cifProcessed ) + "/"+  str(cifNo)+"\n")
    progressSummary.close()
#
mergeLogs("logs/anionPi.log", "logs/anionPi*.log"  )
mergeLogs("logs/cationPi.log", "logs/cationPi*.log"  )
mergeLogs("logs/piPi.log", "logs/piPi*.log"  )
mergeLogs("logs/anionCation.log", "logs/anionCation*.log"  )
mergeLogs("logs/hBonds.log", "logs/hBonds*.log"  )
mergeLogs("logs/additionalInfo.log", "logs/additionalInfo*.log"  )

mergeProgressFiles()

timeStop = time.time()
timeFile = open("logs/timeStop.log", 'w')
timeFile.write(str(timeStop))
timeFile.close()