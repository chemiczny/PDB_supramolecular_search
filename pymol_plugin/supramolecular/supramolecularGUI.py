#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 19:39:10 2018

@author: michal
"""

import pandas as pd
import sys
from os import path

if sys.version_info[0] < 3:
    import Tkinter
    from Tkinter import LEFT, RIGHT
    import tkMessageBox, tkFileDialog
    from pymol import cmd, plugins
    import ttk
else:
    import tkinter as Tkinter
    from tkinter import LEFT, RIGHT
    from tkinter import filedialog as tkFileDialog
    from tkinter import messagebox as tkMessageBox
    import tkinter.ttk as ttk

from cgo_arrow import cgo_arrow

def isStrInt( str2check):
    try:
        int(str2check)
        return True
    except:
        return False
    
class SupramolecularGUI:
    def __init__(self, page):
        self.page = page
        self.logData= {"logFile": False, "data" : None, 
               "cifDir" : None, "arrowExists" : False, "displaying" : False, "displayingAround" : False}
        self.numericalParameters = {}
        self.checkboxVars = {}
        self.listParameters = {}
        self.sorting_keys2header = {}
        self.headers = []
        self.currentMolecule = { "PdbCode" : None  }
        self.additionalCheckboxes = []
        
    def grid(self):
        self.gridLogFilePart()
        self.gridNumericalParameters()
        self.gridAroundPart()
        self.gridAdditionalCheckboxes()
        self.gridListParameters()
        self.gridSortingParameters()
        self.gridSearchWidgets()
        self.gridTree()
        self.gridSaveFiltered()
    
    def gridLogFilePart(self):
        self.but_log = Tkinter.Button(self.page, text = "Load log file", command = self.getLogFile, width = 10)
        self.but_log.grid(row = 0, column = 0)
        
        self.ent_log = Tkinter.Entry(self.page, width =20)
        self.ent_log.configure(state = "readonly")
        self.ent_log.grid(row = 0, column = 1, columnspan = 2)
        
        self.lab_use = Tkinter.Label(self.page, text = "Use filter")
        self.lab_use.grid(row = 0, column = 3)
        
    def getLogFile(self):
        self.logData["logFile"] = tkFileDialog.askopenfilename(title = "Select file", filetypes = (("Log files","*.log"), ("Txt files", "*.txt") ,("CSV files","*.csv"), ("Dat files", "*.dat"), ("all files","*.*")) )
        self.ent_log.configure(state = "normal")
        self.ent_log.delete(0,"end")
        self.ent_log.insert(0, self.logData["logFile"])
        self.ent_log.configure(state = "readonly")
        
        try:
            self.logData["data"] = pd.read_csv(self.logData["logFile"], sep = "\t").fillna("NA")
            self.updateMenu()
        except:
            tkMessageBox.showwarning(title = "Error!", message = "Pandas cannot parse this file")
            
    def updateMenu(self):
        for parameter in self.listParameters:
            parametersData = self.logData["data"][ self.listParameters[parameter]["header"]  ].unique()
            parametersData = sorted(parametersData)
            self.logData[parameter] = parametersData
            
            self.listParameters[parameter]["listbox"].delete(0, "end")
            for row in parametersData:
                self.listParameters[parameter]["listbox"].insert("end", row)
    
    def setNumericalParameters(self, numericalParameters):
        self.numericalParameters = numericalParameters
    
    def gridNumericalParameters(self):
        self.actualRow = 1
    
        for parameter in sorted(self.numericalParameters.keys()):
            self.numericalParameters[parameter]["entry_low"] = Tkinter.Entry(self.page, width = 10)
            self.numericalParameters[parameter]["entry_low"].grid(row = self.actualRow, column = 0)
            
            self.numericalParameters[parameter]["label"] = Tkinter.Label(self.page, width = 10, text = "< "+parameter+" <")
            self.numericalParameters[parameter]["label"].grid(row=self.actualRow, column =1)
            
            self.numericalParameters[parameter]["entry_high"] = Tkinter.Entry(self.page, width = 10)
            self.numericalParameters[parameter]["entry_high"].grid(row = self.actualRow, column = 2)
            
            self.checkboxVars[parameter] = Tkinter.IntVar()
        
            self.numericalParameters[parameter]["checkbox"] = Tkinter.Checkbutton(self.page, variable = self.checkboxVars[parameter] )
            self.numericalParameters[parameter]["checkbox"].grid(row=self.actualRow, column =  3)
            
            self.actualRow+=1
            
    def listFilter(self, key):
        template = self.listParameters[key]["entry"].get()
        template = template.upper()
        
        templateLen = len(template)
        self.listParameters[key]["listbox"].delete(0, "end")
        
        for row in self.logData[key]:
            if template == row[:templateLen]:
                self.listParameters[key]["listbox"].insert("end", row)
            
    def gridAroundPart(self):
        self.lab_around = Tkinter.Label(self.page, text = "Around")
        self.lab_around.grid(row = self.actualRow, column = 0)
        
        self.chkvar_around = Tkinter.IntVar()
        
        self.chk_around = Tkinter.Checkbutton(self.page, variable = self.chkvar_around, command = self.showAround)
        self.chk_around.grid(row = self.actualRow, column = 1)
        
        self.ent_around = Tkinter.Entry(self.page, width = 10)
        self.ent_around.grid(row = self.actualRow, column = 2)
        self.ent_around.insert("end", "5.0")
        
        self.actualRow += 1
        
    def showAround(self):
        if not self.logData["displaying"]:
            return
        
        radius = self.ent_around.get()
        try:
            radius = float(radius)
        except:
            return
        
        if self.chkvar_around.get() > 0 and not self.logData["displayingAround"]:
            selectionAroundName = "suprAround"
            cmd.select( selectionAroundName,  "byres ( suprSelection around "+str(radius)+" ) " )
            cmd.show("lines", selectionAroundName)
            self.logData["displayingAround"] = True
        elif self.chkvar_around.get() == 0 and self.logData["displayingAround"]:
            cmd.hide("lines" , "suprAround")
            cmd.delete("suprAround")
            self.logData["displayingAround"] = False
            
        cmd.deselect()
    
    def setListParameters(self, listParameters):
        self.listParameters = listParameters
        
    def gridListParameters(self):
        self.actualColumn = 4
        for parameter in sorted(self.listParameters.keys()):
            self.listParameters[parameter]["label"] = Tkinter.Label(self.page, text = parameter)
            self.listParameters[parameter]["label"].grid(row = 0, column = self.actualColumn)
            
            self.checkboxVars[parameter] = Tkinter.IntVar()
            
            self.listParameters[parameter]["checkbox"] = Tkinter.Checkbutton(self.page, variable = self.checkboxVars[parameter] )
            self.listParameters[parameter]["checkbox"].grid(row = 0, column = self.actualColumn+1)
            
            self.listParameters[parameter]["listbox"] = Tkinter.Listbox(self.page, width =10, height = 15, exportselection = False)
            self.listParameters[parameter]["listbox"].grid(row = 1, column = self.actualColumn, rowspan = 15, columnspan = 2)
            
            self.listParameters[parameter]["entry"] = Tkinter.Entry(self.page, width = 5)
            self.listParameters[parameter]["entry"].grid(row =20, column = self.actualColumn)
             
            self.listParameters[parameter]["button"] = Tkinter.Button(self.page, width = 1, text = "*", command = lambda arg = parameter : self.listFilter(arg))
            self.listParameters[parameter]["button"].grid(row = 20, column = self.actualColumn + 1)
            
            self.actualColumn += 2
            
    def setSortingParameters(self, keys2header, keysCol):
        self.sorting_keys2header = keys2header
        self.sorting_keys_col1 = keysCol
        self.sorting_keys_col2 = [ "Ascd", "Desc" ]
        
        self.sortingMenu = []
        
    def gridSortingParameters(self):
        if not self.sorting_keys2header:
            return
        
        for i in range(4):
            self.sortingMenu.append( { }  )
            
            self.sortingMenu[i]["label"] = Tkinter.Label(self.page, text = "Sorting"+str(i))
            self.sortingMenu[i]["label"].grid(row = 0, column =  self.actualColumn)
            
            self.sortingMenu[i]["chk_value"] = Tkinter.IntVar()
            
            self.sortingMenu[i]["chk_butt"] = Tkinter.Checkbutton(self.page, variable = self.sortingMenu[i]["chk_value"])
            self.sortingMenu[i]["chk_butt"] .grid(row = 0, column = self.actualColumn+1)
            
            self.sortingMenu[i]["sorting_key"] = Tkinter.IntVar()
            
            actual_row = 1
            for value, key in enumerate(self.sorting_keys_col1):
                self.sortingMenu[i][key] = Tkinter.Radiobutton(self.page, text = key, variable = self.sortingMenu[i]["sorting_key"], value = value, indicatoron=0, width = 8 )
                self.sortingMenu[i][key].grid(row = actual_row, column = self.actualColumn, columnspan = 2)
                actual_row += 1
                
            self.sortingMenu[i]["sortingTypeValue"] = Tkinter.IntVar()
            
            for value, key in enumerate(self.sorting_keys_col2):
                self.sortingMenu[i][key]= Tkinter.Radiobutton(self.page, text = key, variable = self.sortingMenu[i]["sortingTypeValue"], value = value, indicatoron=0, width = 8 )
                self.sortingMenu[i][key].grid(row = actual_row, column = self.actualColumn, columnspan = 2)
                actual_row += 1
                
            self.actualColumn += 2
            
    def gridSaveFiltered(self):
        self.but_saveFiltered = Tkinter.Button(self.page, width = 10, command = self.saveFiltered, text = "Save filtered")
        self.but_saveFiltered.grid(row = 50, column = 1)
        
    def setAdditionalCheckboxes(self, additionalCheckboxes):
        self.additionalCheckboxes = additionalCheckboxes
        
    def gridAdditionalCheckboxes(self):
        if self.additionalCheckboxes:
            for i in range(len(self.additionalCheckboxes)):
                self.additionalCheckboxes[i]["labTk"] = Tkinter.Label(self.page, text = self.additionalCheckboxes[i]["label"])
                self.additionalCheckboxes[i]["labTk"].grid(row = self.actualRow, column =0)
                
                self.additionalCheckboxes[i]["chkVar"] = Tkinter.IntVar()
                
                self.additionalCheckboxes[i]["chk"] = Tkinter.Checkbutton(self.page, variable = self.additionalCheckboxes[i]["chkVar"])
                self.additionalCheckboxes[i]["chk"].grid(row = self.actualRow, column =1)
                
                self.actualRow += 1
        
    def saveFiltered(self):
        if not "filtered" in self.logData:
            return
        
        file2save = tkFileDialog.asksaveasfilename(defaultextension = ".log", filetypes = (("Log files","*.log"), ("Txt files", "*.txt") ,("CSV files","*.csv"), ("Dat files", "*.dat"), ("all files","*.*")) )
        if file2save:
            self.logData["filtered"].to_csv(file2save, sep = "\t")
    
    def applyFilter(self):
        if self.logData["logFile"] == False:
            tkMessageBox.showwarning(title="Warning", message = "Log file not selected")
            return
            
#        anythingSet = False
        actualData = self.logData["data"]
        for key in self.checkboxVars:
            if self.checkboxVars[key].get() > 0:
#                anythingSet = True
                
                if key in self.listParameters:
                    listIndex = self.listParameters[key]["listbox"].curselection()
                    if listIndex:
                        query = str(self.listParameters[key]["listbox"].get(listIndex))
                        actualData = actualData[  actualData[ self.listParameters[key]["header"] ].astype(str) == query ]

                elif key in self.numericalParameters:
                    minValue = self.numericalParameters[key]["entry_low"].get()
                    maxValue = self.numericalParameters[key]["entry_high"].get()

                    try: 
                        minValue = float(minValue)
                        actualData = actualData[  actualData[ self.numericalParameters[key]["header"]  ] > minValue ]
                    except:
                        pass
                    
                    try:
                        maxValue = float(maxValue)
                        actualData = actualData[  actualData[ self.numericalParameters[key]["header"]  ] < maxValue ]
                        
                    except:
                        pass
                    
        for adChk in self.additionalCheckboxes:
            if adChk["chkVar"].get() > 0:
                actualData = adChk["func"](actualData)
                                      
        self.printFilterResults(actualData)
            
    def printFilterResults(self, actualData ):
        recordsFound = str(len(actualData))
        
        anySort = False
        columns = []
        ascending = []
        for sortData in self.sortingMenu:
           toSort =  sortData["chk_value"].get()
           
           if toSort > 0:
               anySort = True
               itemInd = sortData["sorting_key"].get()
               sortingKey = self.sorting_keys_col1[itemInd]
               header = self.sorting_keys2header[sortingKey]
               if header in columns:
                   continue
               columns.append(header)
               
               ascendingActual = sortData["sortingTypeValue"].get()
               if ascendingActual == 0:
                   ascending.append(True)
               else:
                   ascending.append(False)
               
        
        self.ent_recordsFound.configure(state = "normal")
        self.ent_recordsFound.delete(0,"end")
        self.ent_recordsFound.insert(0, str(recordsFound))
        self.ent_recordsFound.configure(state = "readonly")
        if anySort:
            actualData = actualData.sort_values(by = columns, ascending = ascending)
        self.logData["filtered"] = actualData

        self.showRange()
    
    def showInteractions(self):
        if not "filtered" in self.logData:
            return
        currentSel = self.tree_data.focus()
        if currentSel == "" :
            return
        
        rowId = self.tree_data.item(currentSel)["values"][0] 
        data = self.logData["filtered"].iloc[[rowId]]
        pdbCode = data["PDB Code"].values[0]
        
        
        if self.currentMolecule["PdbCode"] and self.currentMolecule["PdbCode"] != pdbCode:
            cmd.delete(self.currentMolecule["PdbCode"])
            
        if self.currentMolecule["PdbCode"] != pdbCode:
            if self.logData["cifDir"] != None:
                potentialPaths = [ path.join( self.logData["cifDir"] ,  pdbCode.lower() +".cif" ), path.join( self.logData["cifDir"] ,  pdbCode.upper() +".cif" )  ]
                cifFound = False
                for filePath in potentialPaths:
                    if path.isfile(filePath):
                        cmd.load(filePath)
                        cifFound = True
                        break
                if not cifFound:
                    cmd.fetch(pdbCode)
            else:    
                cmd.fetch(pdbCode)
        
        frame = int(data["Model No"].values[0])
        if frame != 0:
            cmd.frame(frame+1)
        
        selection = self.getSelection(data)
        
        selectionName = "suprSelection"
        
        cmd.select(selectionName, selection)
        cmd.show( "sticks" , selectionName  )
        cmd.center(selectionName)
        cmd.zoom(selectionName)
        
        if self.chkvar_around.get() > 0:
            selectionAroundName = "suprAround"
            cmd.select( selectionAroundName,  "byres ( suprSelection around 5 ) " )
            cmd.show("lines", selectionAroundName)
            self.logData["displayingAround"] = True
        else:
            self.logData["displayingAround"] = False
        
        if self.logData["arrowExists"]:
            cmd.delete("Anion2Centroid")
            
        arrowBegin, arrowEnd = self.getArrow(data)
        cgo_arrow(arrowBegin, arrowEnd, 0.1, name = "Anion2Centroid")
        
        self.logData["arrowExists"] = True
        self.logData["displaying"] = True
        self.currentMolecule["PdbCode"] = pdbCode
        
        cmd.deselect()
        
    def setSelectionFunc(self, selectionFunc):
        self.getSelection = selectionFunc
        
    def setArrowFunc(self, arrowFunc):
        self.getArrow = arrowFunc
    
    def showRange(self):
        start, stop = self.privGetRange()
        if start == stop:
            return
        
        self.privShowRange(start, stop)
        
    def privGetRange(self):
        if not "filtered" in self.logData:
            return 0, 0
        
        start = self.ent_rangeStart.get()
        stop = self.ent_rangeStop.get()
        
        if not isStrInt(start) or not isStrInt(stop):
            return 0, 1000
        
        start = int(start)
        stop = int(stop)
        
        if stop < start:
            return 0, 1000
        
        if start < 0:
            return 0, 1000
        
        return start, stop
    
    def setRow2Values(self, row2Values):
        self.getValues = row2Values
    
    def privShowRange(self, start, stop):
        self.tree_data.delete(*self.tree_data.get_children())
        
        rowId = 0
        actualData = self.logData["filtered"]
        for index, row in actualData.iterrows():
            if rowId >= start and rowId < stop:
                self.tree_data.insert('', "end" , values = self.getValues(rowId, row)  )
            rowId += 1
            if rowId >= stop:
                break
    
    def showNext(self):
        start, stop = self.privGetRange()
        if start == stop:
            return
        
        diff = stop - start
        newStart = stop
        newStop = stop + diff
        
        dataLen = self.logData["filtered"].shape[0]
        if newStop > dataLen:
            newStop = dataLen
            newStart = dataLen - diff
            
        self.ent_rangeStart.delete(0, "end")
        self.ent_rangeStart.insert("end", str(newStart))
        self.ent_rangeStop.delete(0, "end")
        self.ent_rangeStop.insert("end", str(newStop))
        
        self.privShowRange(newStart, newStop)
    
    def showPrev(self):
        start, stop = self.privGetRange()
        if start == stop:
            return
        
        diff = stop - start
        newStart = start - diff
        newStop = start
        
        if newStart < 0:
            newStop = diff
            newStart = 0
            
        self.ent_rangeStart.delete(0, "end")
        self.ent_rangeStart.insert("end", str(newStart))
        self.ent_rangeStop.delete(0, "end")
        self.ent_rangeStop.insert("end", str(newStop))
        self.privShowRange(newStart, newStop)
            
    def gridSearchWidgets(self):
        self.but_apply = Tkinter.Button(self.page, width = 10, command = self.applyFilter, text = "Search")
        self.but_apply.grid(row = 22, column = 0)
        
        self.lab_data = Tkinter.Label(self.page, width = 10, text = "Records found")
        self.lab_data.grid(row = 25, column = 0)
        
        self.ent_recordsFound = Tkinter.Entry(self.page, width =20)
        self.ent_recordsFound.configure(state = "readonly")
        self.ent_recordsFound.grid(row = 25, column = 1, columnspan = 2)
        
        self.lab_range = Tkinter.Label(self.page, width = 5, text = "Range")
        self.lab_range.grid(row = 25, column = 3 )
        
        self.ent_rangeStart = Tkinter.Entry(self.page, width = 10)
        self.ent_rangeStart.grid(row = 25, column = 4, columnspan = 2)
        self.ent_rangeStart.insert("end", 0)
        
        self.ent_rangeStop = Tkinter.Entry(self.page, width = 10)
        self.ent_rangeStop.grid(row = 25, column = 6, columnspan = 2)
        self.ent_rangeStop.insert("end", 1000)
        
        self.but_showInteraction = Tkinter.Button(self.page, width = 10, command = self.showInteractions, text = "Show interact")
        self.but_showInteraction.grid(row = 50, column =0)
        
        self.but_rangeShow = Tkinter.Button(self.page, width = 6, text = "Show", command = self.showRange)
        self.but_rangeShow.grid(row = 25, column = 8, columnspan = 2)
        
        self.but_rangeNext = Tkinter.Button(self.page, width = 6, text = "Next", command = self.showNext)
        self.but_rangeNext.grid(row = 25, column = 10, columnspan = 2)
        
        self.but_rangePrev = Tkinter.Button(self.page, width = 6, text = "Prev", command = self.showPrev)
        self.but_rangePrev.grid(row = 25, column = 12, columnspan = 2)
        
    def setTreeData(self, treeData):
        self.headers = treeData
        
    def gridTree(self):
        if not self.headers:
            return
        
        self.tree_data = ttk.Treeview(self.page, columns = self.headers, show = "headings", heigh = 15 )
        for header in self.headers:
            self.tree_data.heading(header, text = header)
            self.tree_data.column(header, width = 70)
        self.tree_data.grid(row = 30, column = 0, columnspan = 40)
