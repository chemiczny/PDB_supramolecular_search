#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 14:29:26 2018

@author: michal
"""
import glob
from os.path import isfile
from os import remove
import time
import json
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
    
    
numberOfProcesses = 6

configurationFileName = "config.json"

if isfile(configurationFileName):
    configFile = open(configurationFileName)
    config = json.load(configFile)
    configFile.close()
    
    if "N" in config:
        numberOfProcesses = config["N"]
        
#
pool = Pool(numberOfProcesses)

argumentsList = []
argumentsList.append( ["logs/anionPi.log", "logs/anionPi*.log"]  )
argumentsList.append(["logs/cationPi.log", "logs/cationPi*.log" ] )
argumentsList.append(["logs/piPi.log", "logs/piPi*.log"]  )
argumentsList.append(["logs/anionCation.log", "logs/anionCation*.log"]  )
argumentsList.append(["logs/hBonds.log", "logs/hBonds*.log"]  )
argumentsList.append(["logs/metalLigand.log", "logs/metalLigand*.log"]  )
argumentsList.append(["logs/additionalInfo.log", "logs/additionalInfo*.log"]  )
argumentsList.append(["logs/metalLigand.log", "logs/metalLigand*.log"]  )

pool.map(mergeLogs, argumentsList)

mergeProgressFiles()

timeStop = time.time()
timeFile = open("logs/timeStop.log", 'w')
timeFile.write(str(timeStop))
timeFile.close()