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
from hBonds import HBondsGUI
from metalLigand import MetalLigandGUI

from linearAnionPi import LinearAnionPiGUI
from planarAnionPi import PlanarAnionPiGUI

from QMGUI import QMGUI
from jobStatusGUI import JobStatusGUI
from os import path
import pandas as pd
from cgo_arrow import cgo_arrow
import sys
import json

try:
    from pymol import cmd
except:
    pass

if sys.version_info[0] < 3:
    import tkMessageBox
    import tkFileDialog
else:
    from tkinter import messagebox as tkMessageBox
    from tkinter import filedialog as tkFileDialog
    
class SupramolecularComposition:
    def __init__(self, pageAnionPi, pagePiPi, pageCationPi, pageAnionCation, pageHBonds, pageMetalLigand, pageLinearAnionPi, pagePlanarAnionPi, pageQM, pageJobStatus):
        self.guiAnionPi = AnionPiGUI(pageAnionPi, self.parallelSelect, "AnionPi")
        self.guiPiPi = PiPiGUI(pagePiPi, self.parallelSelect, "PiPi")
        self.guiCationPi = CationPiGUI( pageCationPi, self.parallelSelect, "CationPi")
        self.guiAnionCation = AnionCationGUI( pageAnionCation, self.parallelSelect, "AnionCation")
        self.guiHBonds = HBondsGUI( pageHBonds , self.parallelSelect, "HBonds")
        self.guiMetalLigand = MetalLigandGUI(pageMetalLigand, self.parallelSelect, "MetalLigand")
        self.guiLinearAnionPi = LinearAnionPiGUI(pageLinearAnionPi, self.parallelSelect, "LinearAnionPi")
        self.guiPlanarAnionPi = PlanarAnionPiGUI(pagePlanarAnionPi, self.parallelSelect, "PlanarAnionPi")
        self.guiQM = QMGUI(pageQM)
        self.guiJobStatus = JobStatusGUI(pageJobStatus)
            
        self.guis = [ self.guiAnionPi, self.guiPiPi, self.guiCationPi, self.guiAnionCation, self.guiHBonds, 
                     self.guiMetalLigand, self.guiLinearAnionPi, self.guiPlanarAnionPi ]
        
        self.actionLabels = [ "AnionPi", "PiPi", "CationPi", "AnionCation", "HBonds", "MetalLigand", "LinearAnionPi", "PlanarAnionPi" ]
        self.actionLabels2Objects = { "AnionPi" : self.guiAnionPi, 
                                     "PiPi" : self.guiPiPi, 
                                     "CationPi" : self.guiCationPi,
                                     "AnionCation" : self.guiAnionCation,
                                     "HBonds" : self.guiHBonds,
                                     "MetalLigand" : self.guiMetalLigand,
                                     "LinearAnionPi" : self.guiLinearAnionPi,
                                     "PlanarAnionPi" : self.guiPlanarAnionPi}
        
        self.actualKeys = []
        
    def merge(self, actionMenu):        
        selectedData = []
        data2use = {}
        
        for label in self.actionLabels:
            use = actionMenu[label]["checkValue"].get() > 0 
            exclude = actionMenu[label]["checkValueExclude"].get() > 0 
            if use or exclude  :
                if not "filtered" in  self.actionLabels2Objects[label].logData:
                    self.actionLabels2Objects[label].dataIsMerged = False
                    continue
                else:
                    if use and exclude:
                        self.actionLabels2Objects[label].dataIsMerged = False
                        continue
                    
                    data2use[label] = use
                    selectedData.append(label)
            else:
                self.actionLabels2Objects[label].dataIsMerged = False
                    
        if len(selectedData) < 2:
            return
        
        uniqueData = []
        dataExcluded =[]
        self.actualKeys = []
        excludedKeys = []
        
        for key in selectedData:                
            headers = self.getHeadersForGUI(key)
            
            
            if data2use[key]:
                if len(uniqueData) == 0:
                    uniqueData = self.actionLabels2Objects[ key ].logData["filtered"][ headers ].drop_duplicates()
                elif len(self.actionLabels2Objects[ key ].logData["filtered"]) > 0:
                    uniqueData = pd.merge( uniqueData,  self.actionLabels2Objects[ key ].logData["filtered"][ headers ], on = list( set(self.actualKeys) & set(headers) ))
                    uniqueData = uniqueData.drop_duplicates()
                self.actualKeys = list(set( self.actualKeys + headers ))
                
            else:
                if len(dataExcluded) == 0:
                    dataExcluded = self.actionLabels2Objects[ key ].logData["filtered"][ headers ].drop_duplicates()
                elif len(self.actionLabels2Objects[ key ].logData["filtered"]) > 0:
                    dataExcluded = pd.merge( dataExcluded,  self.actionLabels2Objects[ key ].logData["filtered"][ headers ], on = list( set(self.actualKeys) & set(headers) ))
                    dataExcluded = dataExcluded.drop_duplicates()
                excludedKeys = list(set( excludedKeys + headers ))
                
                
        if len(dataExcluded) > 0:
            mergingKeys = list(set(self.actualKeys) & set(excludedKeys) )
            subMerged = pd.merge( uniqueData,  dataExcluded , on = mergingKeys, how='left', indicator=True )
            uniqueData = subMerged[ subMerged['_merge'] == 'left_only' ]
            
        
        for key in selectedData:
            headers = self.getHeadersForGUI(key)
                    
            if len(uniqueData) == 0 :
                tkMessageBox.showwarning(title = "Merging error!", message = "No data left after merge!")
                break
            
            mergingKeys = list(set(self.actualKeys) & set(headers) )
            tempDataFrame = uniqueData[ mergingKeys   ].drop_duplicates()
            mergedData = pd.merge( self.actionLabels2Objects[ key ].logData["filtered"], tempDataFrame, on = mergingKeys )
            self.actionLabels2Objects[ key ].logData["filtered"] = mergedData
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
        uniqueData = self.actionLabels2Objects[lastSelectionMenu].actualDisplaying["rowData"]
                
        lastDataHeaders = self.getHeadersForGUI(lastSelectionMenu)
        
        selection = ""
        keySum = set(lastDataHeaders)
        for key in selectedData:
            headers = self.getHeadersForGUI(key)
            
            mergingKeys = list(keySum & set(headers) )
            uniqueData = pd.merge( self.actionLabels2Objects[ key ].logData["filtered"], uniqueData, on = mergingKeys )
            keySum |= set(headers)
            
        for key in selectedData:
            headers = self.getHeadersForGUI(key)
            mergingKeys = list(set(headers) & keySum)
            tempDataFrame = uniqueData[ mergingKeys   ].drop_duplicates()
            mergedData = pd.merge( self.actionLabels2Objects[ key ].logData["filtered"], tempDataFrame, on = mergingKeys )
            for index, row in mergedData.iterrows():
                arrowBegin, arrowEnd = self.actionLabels2Objects[ key ].getArrowFromRow(row)
                uniqueArrowName = self.actionLabels2Objects[ key ].arrowName+"A"+str(index)
                cgo_arrow(arrowBegin, arrowEnd, 0.1, name = uniqueArrowName, color = self.actionLabels2Objects[ key ].arrowColor )
                
                if selection == "":
                    selection = "(" + self.actionLabels2Objects[key].getSelectionFromRow(row) + " ) "
                else:
                    selection += "or (" + self.actionLabels2Objects[key].getSelectionFromRow(row) + " ) "
                    
        selectionName = "suprSelectionFull"
        cmd.select(selectionName, selection)
        cmd.show( "sticks" , selectionName  )
        cmd.center(selectionName)
        cmd.zoom(selectionName)
                
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
            
        self.guiQM.grid()
        self.guiJobStatus.grid()
        
    def loadState(self):
        jsonFile = tkFileDialog.askopenfilename(title = "Select file", filetypes = (("Json files","*.json"), ("all files","*.*")) )
        
        if not jsonFile:
            return
            
        with open(jsonFile, 'r') as fp:
            state = json.load(fp)
        
        for label in self.actionLabels2Objects:
            if label in state:
                self.actionLabels2Objects[label].loadState( state[label] )
                
        if "QMGUI" in state:
            self.guiQM.loadState(state["QMGUI"])
            
        if "JobStatus" in state:
            self.guiJobStatus.loadState(state["JobStatus"])
    
    def saveState(self):
        state = {}
        for label in self.actionLabels2Objects:
            state[label] = self.actionLabels2Objects[label].getState()
            
        state["QMGUI"] = self.guiQM.getState()
        state["JobStatus"] = self.guiJobStatus.getState()
            
        file2save = tkFileDialog.asksaveasfilename(defaultextension = ".json", filetypes = (("Json files","*.json"),  ("all files","*.*")) )
        if file2save:
            with open(file2save, 'w') as fp:
                json.dump(state, fp)
                
    def setParallelSelection(self, state):
        if state == 1:
            for gui in self.guis:
                gui.parallelSelection = True
        else:
            for gui in self.guis:
                gui.parallelSelection = False
    
    def getHeadersForGUI(self, guiKey):
        lastDataHeaders = [ "PDB Code", "Model No" ] 
        
        headersId ={ 
                    "Pi" : [ "Pi acid Code" , "Pi acid chain" , "Piacid id", "CentroidId"] , 
                    "Anion" : ["Anion code", "Anion chain" , "Anion id"] , 
                    "Cation" : ["Cation code", "Cation chain", "Cation id"] ,
                    "HBonds" : ["Anion code", "Anion chain" , "Anion id"],
                    "Metal" : ["Cation code", "Cation chain", "Cation id"],
                    "Ligand" : ["Anion code", "Anion chain" , "Anion id"] }
        
        for headerKey in headersId:
            if headerKey in guiKey:
                lastDataHeaders += headersId[headerKey]
                
        if "Anion" in guiKey and "Cation" in guiKey:
            lastDataHeaders += headersId["Pi"]
                
        return lastDataHeaders
                
    def parallelSelect(self, name, selectedRow):
        selectedHeaders = self.getHeadersForGUI(name)
        for guiName in self.actionLabels2Objects:
            if guiName == name:
                continue
            
            if not self.actionLabels2Objects[guiName].dataIsMerged:
                continue
            
            headers = self.getHeadersForGUI(guiName)
            commonHeaders = set(selectedHeaders) & set(headers)
#            selectedHeaders = list(set(selectedHeaders) | set(headers))
            header2value = {}
            
            for head in commonHeaders:
                header2value[head] = selectedRow[head].values[0]
                
            self.actionLabels2Objects[guiName].selectRowInTree(header2value)
            