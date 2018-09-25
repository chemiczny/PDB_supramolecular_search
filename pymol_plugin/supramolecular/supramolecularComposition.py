#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:31:51 2018

@author: michal
"""
from anionPi import AnionPiGUI
from piPi import PiPiGUI
from cationPi import CationPiGUI
from anionCation import AnionCationGUI
from os import path
import pandas as pd
from cgo_arrow import cgo_arrow

import sys
if sys.version_info[0] < 3:
    from pymol import cmd
    
class SupramolecularComposition:
    def __init__(self, pageAnionPi, pagePiPi, pageCationPi, pageAnionCation):
        self.guiAnionPi = AnionPiGUI(pageAnionPi)
        self.guiPiPi = PiPiGUI(pagePiPi)
        self.guiCationPi = CationPiGUI( pageCationPi)
        self.guiAnionCation = AnionCationGUI( pageAnionCation)
            
        self.guis = [ self.guiAnionPi, self.guiPiPi, self.guiCationPi, self.guiAnionCation]
        self.actionLabels = [ "AnionPi", "PiPi", "CationPi", "AnionCation" ]
        self.actionLabels2Objects = { "AnionPi" : self.guiAnionPi, 
                                     "PiPi" : self.guiPiPi, 
                                     "CationPi" : self.guiCationPi,
                                     "AnionCation" : self.guiAnionCation }
        
        self.actualKeys = []
        
    def merge(self, actionMenu):
        headersId ={ 
                    "Pi" : [ "Pi acid Code" , "Pi acid chain" , "Piacid id", "CentroidId"] , 
                    "Anion" : ["Anion code", "Anion chain" , "Anion id"] , 
                    "Cation" : ["Cation code", "Cation chain", "Cation id"]  }
        
        selectedData = []
        
        
        for label in self.actionLabels:
            if actionMenu[label]["checkValue"].get() > 0:
                if not "filtered" in  self.actionLabels2Objects[label].logData:
                    continue
                else:
                    selectedData.append(label)
            else:
                self.actionLabels2Objects[label].dataIsMerged = False
                    
        if len(selectedData) < 2:
            return
        
        uniqueData = []
        self.actualKeys = []
        for key in selectedData:
            headers = [ "PDB Code", "Model No" ] 
            
            for headerKey in headersId:
                if headerKey in key:
                    headers += headersId[headerKey]
                    
            if "Anion" in key and "Cation" in key:
                headers += headersId["Pi"]
                    
            if len(uniqueData) == 0:
                uniqueData = self.actionLabels2Objects[ key ].logData["filtered"][ headers ].drop_duplicates()
            else:
                uniqueData = pd.merge( uniqueData,  self.actionLabels2Objects[ key ].logData["filtered"][ headers ], on = list( set(self.actualKeys) & set(headers) ))
                uniqueData = uniqueData.drop_duplicates()
                
            self.actualKeys = list(set( self.actualKeys + headers ))
            
            
        for key in selectedData:
            headers = [ "PDB Code", "Model No" ] 
            
            for headerKey in headersId:
                if headerKey in key:
                    headers += headersId[headerKey]
                    
            if "Anion" in key and "Cation" in key:
                headers += headersId["Pi"]
                    
            mergingKeys = list(set(self.actualKeys) & set(headers) )
            tempDataFrame = uniqueData[ mergingKeys   ].drop_duplicates()
            mergedData = pd.merge( self.actionLabels2Objects[ key ].logData["filtered"], tempDataFrame, on = mergingKeys )
            self.actionLabels2Objects[ key ].printFilterResults( mergedData )
            self.actionLabels2Objects[ key ].dataIsMerged = True
            
            
    
    def showAll(self, showMenu):
        if not self.actualKeys:
            return
        
        selectedData = []
        lastSelectionMenu = ""
        lastSelectionTime = -1
        
        for label in self.actionLabels:
            newTime =  self.actionLabels2Objects[label].actualDisplaying["selectionTime"]
            if newTime > lastSelectionTime:
                lastSelectionTime = newTime
                lastSelectionMenu = label
                
            if showMenu[label]["checkValue"].get() > 0:
                if not "filtered" in  self.actionLabels2Objects[label].logData:
                    continue
                else:
                    selectedData.append(label)
                    
        if len(selectedData) < 1:
            return
        
        if lastSelectionTime < 0:
            return
        
        if not lastSelectionMenu in selectedData:
            return
        
        self.deleteMergedArrows()
        lastData = self.actionLabels2Objects[lastSelectionMenu].actualDisplaying["rowData"]
        
        lastDataHeaders = [ "PDB Code", "Model No" ] 
        
        headersId ={ 
                    "Pi" : [ "Pi acid Code" , "Pi acid chain" , "Piacid id", "CentroidId"] , 
                    "Anion" : ["Anion code", "Anion chain" , "Anion id"] , 
                    "Cation" : ["Cation code", "Cation chain", "Cation id"]  }
        
        for headerKey in headersId:
            if headerKey in lastSelectionMenu:
                lastDataHeaders += headersId[headerKey]
        
        for key in selectedData:
            if key == lastSelectionMenu:
                continue
            
            
            headers = [ "PDB Code", "Model No" ] 
            
            for headerKey in headersId:
                if headerKey in key:
                    headers += headersId[headerKey]
                    
            if "Anion" in key and "Cation" in key:
                headers += headersId["Pi"]
                    
            mergingKeys = list(set(lastDataHeaders) & set(headers) )
            tempDataFrame = lastData[ mergingKeys   ]
            mergedData = pd.merge( self.actionLabels2Objects[ key ].logData["filtered"], tempDataFrame, on = mergingKeys )
            
            for index, row in mergedData.iterrows():
                arrowBegin, arrowEnd = self.actionLabels2Objects[ key ].getArrowFromRow(row)
                uniqueArrowName = self.actionLabels2Objects[ key ].arrowName+"A"+str(index)
                cgo_arrow(arrowBegin, arrowEnd, 0.1, name = uniqueArrowName, color = self.actionLabels2Objects[ key ].arrowColor )
                
    def deleteMergedArrows(self):
        for arrow in cmd.get_names_of_type("object:cgo"):
            if "rrowA" in arrow:
                cmd.delete(arrow)
        
    
    def readAllLogsFromDir(self, logDir):
        for guiKey in self.actionLabels2Objects:
            basename = guiKey[0].lower()+guiKey[1:]+".log"
            logFileName = path.join( logDir, basename )
            if path.isfile(logFileName):
                self.actionLabels2Objects[guiKey].logData["logFile"] = logFileName
                self.actionLabels2Objects[guiKey].openLogFile()
                
    def selectCifDir(self, cifDir):
        for gui in self.guis:
            gui.logData["cifDir"] = cifDir
            
    def grid(self):
        for gui in self.guis:
            gui.grid()
            