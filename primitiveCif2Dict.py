#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 14:22:30 2018

@author: michal
"""
import shlex 

class primitiveCif2Dict:
    def __init__(self, cif, interestingKeys):
        cifFile = open(cif, 'r')
        self.result = {}
        self.interestingKeys = interestingKeys
        
        line = cifFile.readline()
        self.loop = False
        self.loopKeysSize = 0
        self.interestingKeyIndex = {}
        
        while line:
            if self.loop:
                self.loopCase(line)
            else:
                self.outOfTheLoopCase(line)
            
            line = cifFile.readline()
        
        cifFile.close()
        
    def loopCase(self, line):
        lineSpl = line.split()
        if len(lineSpl) == 1 and lineSpl[0][0] == "_":
            self.loopKeysSize +=1
            key = self.interestingKeyInLine(line)
            if key:
                self.interestingKeyIndex[key] = self.loopKeysSize-1
        else:
            if self.interestingKeyIndex:
                lineSpl = shlex.split(line)
                if len(lineSpl) == self.loopKeysSize:
                    if self.interestingKeyIndex:
                        for key in self.interestingKeyIndex:
                            self.appendValue2Key(key, lineSpl[ self.interestingKeyIndex[key] ])
            else:
                self.loop = False
                self.loopKeysSize = 0
                self.interestingKeyIndex = {}
    
    def outOfTheLoopCase(self, line):
        if "loop_" in line.lower():
            self.loop = True
            return
        elif ";" == line.strip()[0]:
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
            self.result[key] = value
        else:
            self.result[key] = [ self.result[key] , value ]
    
    
                    
#if __name__ == "__main__":
#    cifs = glob("cif/*cif")
#    timeStart = time.time()
#    
#    for cif in cifs:
#        print(cif)
#        print(readResolutionAndMethod(cif))
#        
#    timeStop = time.time()
#    
#    print(timeStop-timeStart)
#    
#    timeStart = time.time()
#    
#    for cif in cifs:
#        print(cif)
#        test = primitiveCif2Dict(cif, ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" , "_exptl.method" ])
#        print(test.result)
#        
#    timeStop = time.time()
#    
#    print(timeStop-timeStart)
    