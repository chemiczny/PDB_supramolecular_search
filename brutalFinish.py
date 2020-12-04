#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:29:26 2018

@author: michal
"""
import glob
from os.path import join
from os import remove
import time
from configure import configure
from multiprocessing import Pool

def mergeLogs( logList  ):
    logFinal = logList[0]
    logs = logList[1]
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
    
def mergeProgressFiles(scratch):
    cifProcessed = 0
    log_files = glob.glob(join(scratch, "partialProgress*"))
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
    
    
config = configure()
numberOfProcesses = config["N"]
scratch = config["scratch"]
        
#
pool = Pool(numberOfProcesses)

argumentsList = []
argumentsList.append( ["logs/anionPi.log", join(scratch, "anionPi*.log")]  )
argumentsList.append(["logs/cationPi.log", join(scratch, "cationPi*.log") ]  )
argumentsList.append(["logs/piPi.log", join(scratch, "piPi*.log")])  
argumentsList.append(["logs/anionCation.log", join(scratch, "anionCation*.log")] )
argumentsList.append(["logs/hBonds.log", join(scratch, "hBonds*.log")]  )
argumentsList.append(["logs/metalLigand.log", join(scratch, "metalLigand*.log")]  )
argumentsList.append(["logs/planarAnionPi.log", join(scratch, "planarAnionPi*.log" )] )
argumentsList.append(["logs/linearAnionPi.log", join(scratch, "linearAnionPi*.log") ] )
argumentsList.append(["logs/methylPi.log", join(scratch, "methylPi*.log") ] )
argumentsList.append(["logs/additionalInfo.log", join(scratch, "additionalInfo*.log")] )

pool.map(mergeLogs, argumentsList)

mergeProgressFiles(scratch)

timeStop = time.time()
timeFile = open("logs/timeStop.log", 'w')
timeFile.write(str(timeStop))
timeFile.close()
