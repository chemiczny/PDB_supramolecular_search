#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 14:22:30 2018

@author: michal
"""
import shlex 

class primitiveCif2Dict:
    def __init__(self, cif, interestingKeys):
        self.cifFile = open(cif, 'r')
        self.result = {}
        self.interestingKeys = interestingKeys
        
        line = self.cifFile.readline()
        self.loop = False
        self.loopKeysSize = 0
        self.interestingKeyIndex = {}
        
        while line:
            if self.loop:
                self.loopCase(line)
            else:
                self.outOfTheLoopCase(line)
            
            line = self.cifFile.readline()
        
        self.cifFile.close()
        
    def loopCase(self, line):
        if "loop_" in line.lower():
            self.loop = True
            self.loopKeysSize = 0
            self.interestingKeyIndex = {}
            return
        
        if line.startswith("#"):
            self.loop = False
            return 
        
        lineSpl = line.split()
        if len(lineSpl) == 1 and lineSpl[0][0] == "_":
            self.loopKeysSize +=1
            key = self.interestingKeyInLine(line)
            if key:
                self.interestingKeyIndex[key] = self.loopKeysSize-1
        else:
            if self.interestingKeyIndex:
                loopData = []
                while not "_" in line and not "#" in line:
                                    
                    if not line.startswith(";"):
                        lineSpl = shlex.split(line)
                        loopData += lineSpl
                    else:
                        colonCounter = line.count(";")
                        newData = ""
                        while colonCounter < 2:
                            newData += line.strip()
                            line = self.cifFile.readline()
                            colonCounter += line.count(";")
                            
                        newData += line.strip()
                        loopData.append(newData)
    
                    line = self.cifFile.readline()
                    
                loopSize = len(loopData)
                for key in self.interestingKeyIndex:
                    index2look = self.interestingKeyIndex[key] 
                    while index2look < loopSize:
                        self.appendValue2Key(key, loopData[ index2look ])
                        index2look += self.loopKeysSize
                        
                
                self.loop = False
                self.loopKeysSize = 0
                self.interestingKeyIndex = {}
    
    def outOfTheLoopCase(self, line):
        if "loop_" in line.lower():
            self.loop = True
            self.loopKeysSize = 0
            self.interestingKeyIndex = {}
            return
        elif line.startswith(";"):
            return
        else:
            if not self.fastInterestingKeyInLine(line):
                return
            
            lineSpl = shlex.split(line)
            key = self.interestingKeyInLine(lineSpl[0])
            if key:
                self.appendValue2Key(key, lineSpl[-1])
        
    def fastInterestingKeyInLine(self, line):
        for key in self.interestingKeys:
            if key in line:
                return True
            
        return False
    
    def interestingKeyInLine(self, line):
        lineStrip =line.strip()
        for key in self.interestingKeys:
            if key == lineStrip:
                return key
            
        return False
    
    def appendValue2Key(self, key, value):
        if not key in self.result:
            self.result[key] = [ value ]
        else:
            self.result[key].append( value )
    
    
                    
if __name__ == "__main__":
    
    cif = "cif/2n5t.cif"
    print(cif)
    test = primitiveCif2Dict(cif, ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" , "_exptl.method" ])
    print(test.result)
        
    