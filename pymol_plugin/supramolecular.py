#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:41:18 2018

@author: michal

When I wrote this code, only God 
and I undestood what it does.

Now, only God knows.
"""

import sys
import pandas as pd
from os import path
if sys.version_info[0] < 3:
    import Tkinter
    from Tkinter import LEFT, RIGHT
    import tkMessageBox, tkFileDialog
    from pymol import cmd, plugins, CmdException, cgo
    import ttk
else:
    import tkinter as Tkinter
    from tkinter import LEFT, RIGHT
    from tkinter import filedialog as tkFileDialog
    from tkinter import messagebox as tkMessageBox
    import tkinter.ttk as ttk
    

def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
              color='blue red', name=''):
    '''
DESCRIPTION

    Create a CGO arrow between two picked atoms.

ARGUMENTS

    atom1 = string: single atom selection or list of 3 floats {default: pk1}

    atom2 = string: single atom selection or list of 3 floats {default: pk2}

    radius = float: arrow radius {default: 0.5}

    gap = float: gap between arrow tips and the two atoms {default: 0.0}

    hlength = float: length of head

    hradius = float: radius of head

    color = string: one or two color names {default: blue red}

    name = string: name of CGO object
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)
    
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
                        query = self.listParameters[key]["listbox"].get(listIndex)
                        actualData = actualData[  actualData[ self.listParameters[key]["header"] ] == query ]

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


def __init_plugin__(self=None):
    plugins.addmenuitem('Supramolecular analyser', fetchdialog)


def fetchdialog(simulation = False):
    if simulation:
        root = Tkinter.Tk()
    else:
        app = plugins.get_pmgapp()
        root = plugins.get_tk_root()
        
    appData = { "general" : {} , "anionPi" : {} , "piPi" : {}, "cationPi" : {} }
    
    self = Tkinter.Toplevel(root)
    self.title('Supramolecular analyser')
    self.minsize(1100, 700)
    self.resizable(0,0)
    
    nb = ttk.Notebook(self, height = 700, width = 1100)
    
    pageGeneral = ttk.Frame(nb)
    pageAnionPi = ttk.Frame(nb)
    pagePiPi = ttk.Frame(nb)
    pageCationPi = ttk.Frame(nb)
    pageAnionCation  = ttk.Frame(nb)
    
    nb.add(pageGeneral, text = "General")
    nb.add(pageAnionPi, text = "AnionPi")
    nb.add(pagePiPi, text = "PiPi")
    nb.add(pageCationPi, text = "CationPi")
    nb.add(pageAnionCation, text = "AnionCation")
    
    nb.grid(column = 0)
    
    guiAnionPi = SupramolecularGUI(pageAnionPi)
    guiPiPi = SupramolecularGUI(pagePiPi)
    guiCationPi = SupramolecularGUI( pageCationPi)
    guiAnionCation = SupramolecularGUI( pageAnionCation)
    
    ######################
    # GENERAL
    ######################
    
    def selectCif():
        appData["general"]["cifDir"] = tkFileDialog.askdirectory()
        ent_cifDir.configure(state = "normal")
        ent_cifDir.delete(0,"end")
        ent_cifDir.insert(0, appData["general"]["cifDir"])
        ent_cifDir.configure(state = "readonly")
    
    but_cifDir = Tkinter.Button(pageGeneral, width = 10, command = selectCif, text = "Cif dir")
    but_cifDir.grid(row = 0, column = 1)
    
    ent_cifDir = Tkinter.Entry(pageGeneral, width =45)
    ent_cifDir.configure(state = "readonly")
    ent_cifDir.grid(row = 0, column = 2, columnspan = 3)
    
    actions = [ "AND", "AND NOT" ]
    actionLabels = [ "AnionPi", "PiPi", "CationPi", "AnionCation" ]
    actionLabels2Objects = { "AnionPi" : guiAnionPi, "PiPi" : guiPiPi, "CationPi" : guiCationPi, "AnionCation" : guiAnionCation }
    actionMenu = {}
    
    column = 1
    
    lab_usePage = Tkinter.Label(pageGeneral, width = 10, text = "Use:")
    lab_usePage.grid(row = 2, column = 0)
    for label in actionLabels:
        actionMenu[label] = {}
        
        actionMenu[label]["label"] = Tkinter.Label(pageGeneral, text = label)
        actionMenu[label]["label"].grid(row = 1, column = column )
        
        actual_row = 2
        
        actionMenu[label]["checkValue"] = Tkinter.IntVar()
        actionMenu[label]["checkbox"] = Tkinter.Checkbutton(pageGeneral, variable = actionMenu[label]["checkValue"])
        actionMenu[label]["checkbox"].grid(row = actual_row, column = column)
        
        actual_row += 1
        
        actionMenu[label]["radioButton"] = {}
        actionMenu[label]["actionKey"] = Tkinter.IntVar()
        
        for value, key in enumerate(actions):
            actionMenu[label]["radioButton"][key]= Tkinter.Radiobutton(pageGeneral, text = key, value = value, variable = actionMenu[label]["actionKey"], indicatoron=0, width = 8 )
            actionMenu[label]["radioButton"][key].grid(row = actual_row, column = column)
            actual_row += 1
            
        column += 1
    
    
    def mergeResults():
        dataMerged = False
        for label in actionLabels:
            if actionMenu[label]["checkValue"].get() > 0:
                if not "filtered" in  actionLabels2Objects[label].logData:
                    continue
                elif not dataMerged:
                    dataMerged = actionLabels2Objects[label].logData
                else:
                    pass
    
    but_merge = Tkinter.Button(pageGeneral, width = 20, text = "Merge!", command = mergeResults)
    but_merge.grid(row = 10, column = 0, columnspan = 2)
    
    headers = [ "ID" , "PDB" , "Pi acid", "Pi acid id", "Anion", "Anion id", "Anion type" , "R", "alpha", "x", "h", "res", "Method" ]
    
    tree_data = ttk.Treeview(pageGeneral, columns = headers, show = "headings", heigh = 15 )
    for header in headers:
        tree_data.heading(header, text = header)
        tree_data.column(header, width = 70)
    tree_data.grid(row = 30, column = 0, columnspan = 40)
    ######################
    # HELPERS
    ######################
    
    def noAAinPiAcids(actualData):
        return actualData[(actualData["Pi acid Code"] != "TYR") & (actualData[ "Pi acid Code" ] != "PHE" ) & ( actualData["Pi acid Code"] != "HIS" ) & ( actualData["Pi acid Code"] != "TRP" ) ]
    
    def noAAinPiRes(actualData):
        return actualData[(actualData["Pi res code"] != "TYR") & (actualData[ "Pi res code" ] != "PHE" ) & ( actualData["Pi res code"] != "HIS" ) & ( actualData["Pi res code"] != "TRP" ) ]
    
    def noAAinAnions(actualData):
        return actualData[(actualData["Anion code"] != "GLU") & (actualData[ "Anion code" ] != "ASP" ) & ( actualData["Anion code"] != "TYR" ) & ( actualData["Anion code"] != "CYS") ]
   
    
    ######################
    # ANION PI
    ######################
    
    guiAnionPi.setNumericalParameters({ "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" },
                           "alpha" : { "header" : "Angle" }, "res" : { "header" : "Resolution" }  }  )
    
    guiAnionPi.setListParameters( { "Pi acid" : { "header" : "Pi acid Code" }, 
                      "Anions" : { "header" : "Anion code" } , 
                      "Groups" : { "header" : "Anion type" },
                      "Methods" : { "header" : "Method" }} )
    
    guiAnionPi.setSortingParameters({  "R" : "Distance", "Angle" : "Angle", "x" : "x" , "h" : "h",
                         "res" : "Resolution", "Pi acid" : "Pi acid Code" ,"Anion" : "Anion code" }, [  "R" , "Angle" , "x" , "h" ,
                         "res" , "Pi acid" ,"Anion" ] )
    
    guiAnionPi.setTreeData([ "ID" , "PDB" , "Pi acid", "Pi acid id", "Anion", "Anion id", "Anion type" , "R", "alpha", "x", "h", "res", "Method" ])
     
    guiAnionPi.setAdditionalCheckboxes( [ { "label" : "No AA in Pi acids", "func" : noAAinPiAcids } ,
                                         { "label" : "No AA in anions" ,  "func" : noAAinAnions }  ]  )
    
    def row2ValuesAnionPi(rowId, row):
        return ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                        row["Pi acid chain"]+str(row["Piacid id"]) , row["Anion code"], 
                                                        row["Anion chain"] + str(row["Anion id"]), row["Anion type"], 
                                                        str(row["Distance"])[:3], str(row["Angle"])[:4], str(row["x"])[:3],
                                                        str(row["h"])[:3], row["Resolution"] , row["Method"]  )
        
    guiAnionPi.setRow2Values( row2ValuesAnionPi )
    
    def getSelectionAnionPi( data):
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Anion id"].values[0]
        res2Chain = data["Anion chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    guiAnionPi.setSelectionFunc( getSelectionAnionPi )
    
    def getArrowAnionPi( data ):
        centroidCoords = [ data["Centroid x coord"].values[0] , data["Centroid y coord"].values[0] , data["Centroid z coord"].values[0] ]
        anionAtomCoords = [ data["Anion x coord"].values[0] , data["Anion y coord"].values[0] , data["Anion z coord"].values[0] ]
        
        return centroidCoords, anionAtomCoords
    
    guiAnionPi.setArrowFunc( getArrowAnionPi )
    
    ######################
    # PI PI
    ######################
    
    guiPiPi.setNumericalParameters( { "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" },
                           "alpha" : { "header" : "Angle" }, "theta" : { "header" : "theta" }  }  )
    
    guiPiPi.setListParameters( { "Pi 1" : { "header" : "Pi acid Code" }, 
                      "Pi 2" : { "header" : "Pi res code" } } )
    
    guiPiPi.setSortingParameters({  "R" : "Distance", "Angle" : "Angle", "x" : "x" , "h" : "h",
                         "theta" : "theta", "Pi 1" : "Pi acid Code" ,"Pi 2" : "Pi res code" }, [  "R" , "Angle" , "x" , "h" ,
                         "theta" , "Pi 1" ,"Pi 2" ] )
    
    guiPiPi.setTreeData([ "ID" , "PDB" , "Pi 1", "Pi 1 id", "Pi 2", "Pi 2 id" , "R", "alpha", "x", "h", "theta" ])
    
    guiPiPi.setAdditionalCheckboxes( [ { "label" : "No AA in Pi 1", "func" : noAAinPiAcids } ,
                                         { "label" : "No AA in Pi 2" ,  "func" : noAAinPiRes }  ]  )
    
    def row2ValuesPiPi(rowId, row):
        return ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                        row["Pi acid chain"]+str(row["Piacid id"]) , row["Pi res code"], 
                                                        row["Pi res chain"] + str(row["Pi res id"]), 
                                                        str(row["Distance"])[:3], str(row["Angle"])[:4], str(row["x"])[:3],
                                                        str(row["h"])[:3], row["theta"]  )
    guiPiPi.setRow2Values( row2ValuesPiPi )
    
    def getSelectionPiPi( data):
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Pi res id"].values[0]
        res2Chain = data["Pi res chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    guiPiPi.setSelectionFunc( getSelectionPiPi )
    
    def getArrowPiPi( data ):
        point1Coords = [ data["Centroid x coord"].values[0] , data["Centroid y coord"].values[0] , data["Centroid z coord"].values[0] ]
        point2Coords = [ data["Centroid 2 x coord"].values[0] , data["Centroid 2 y coord"].values[0] , data["Centroid 2 z coord"].values[0] ]
        
        return point1Coords, point2Coords
    
    guiPiPi.setArrowFunc( getArrowPiPi )
    
    ######################
    # CATION PI
    ######################
    
    guiCationPi.setNumericalParameters( { "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" },
                           "alpha" : { "header" : "Angle" } } )
    
    guiCationPi.setListParameters({ "Pi acid" : { "header" : "Pi acid Code" }, 
                      "Cation" : { "header" : "Cation code" },
                      "Element" : { "header" : "Atom symbol" }
                      })
    
    
    guiCationPi.setSortingParameters({  "R" : "Distance", "Angle" : "Angle", "x" : "x" , "h" : "h",
                          "Pi acid" : "Pi acid Code" ,"Cation" : "Cation code", "Cat. el." : "Atom symbol" }, [  "R" , "Angle" , "x" , "h" ,
                          "Pi acid" ,"Cation" , "Cat. el." ] )
    
    guiCationPi.setTreeData([ "ID" , "PDB" , "Pi acid", "Pi acid id", "Cation", "Cation id", "Cat. el." , "R", "alpha", "x", "h" ])
    
    guiCationPi.setAdditionalCheckboxes( [ { "label" : "No AA in Pi acids", "func" : noAAinPiAcids }   ]  )
    
    def getSelectionCationPi( data):
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Cation id"].values[0]
        res2Chain = data["Cation chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    guiCationPi.setSelectionFunc( getSelectionCationPi )
    
    def getArrowCationPi( data ):
        point1Coords = [ data["Centroid x coord"].values[0] , data["Centroid y coord"].values[0] , data["Centroid z coord"].values[0] ]
        point2Coords = [ data["Cation x coord"].values[0] , data["Cation y coord"].values[0] , data["Cation z coord"].values[0] ]
        
        return point1Coords, point2Coords
    
    guiCationPi.setArrowFunc( getArrowCationPi )
    
    def row2ValuesCationPi(rowId, row):
        return ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                        row["Pi acid chain"]+str(row["Piacid id"]) , row["Cation code"], 
                                                        row["Cation chain"] + str(row["Cation id"]), row["Atom symbol"], 
                                                        str(row["Distance"])[:3], str(row["Angle"])[:4], str(row["x"])[:3],
                                                        str(row["h"])[:3]  )
        
    guiCationPi.setRow2Values(row2ValuesCationPi)
    
    ######################
    # CATION ANION
    ######################
    
    guiAnionCation.setNumericalParameters({ "R" : {"header" : "Distance"} })
    
    guiAnionCation.setListParameters({ "Cation" : { "header" : "Cation Code" }, 
                      "Cat. el." : { "header" : "Cation symbol" } , 
                      "Anion" : { "header" : "Anion code" },
                      "An. at." : { "header" : "Anion symbol"}})
    
    guiAnionCation.setSortingParameters({  "R" : "Distance", "Cation" : "Cation Code", "Cat. el." : "Cation symbol" , "Anion" : "Anion code",
                         "An. el." : "Anion symbol" }, [  "R" , "Cation" , "Cat. el." , "Anion" ,
                         "An. el."  ])
    
    guiAnionCation.setTreeData([ "ID" , "PDB" , "Cation", "Cation id", "Anion", "Anion id", "Anion el.", "Cat. el." , "R"])
    
    guiAnionCation.setAdditionalCheckboxes( [ { "label" : "No AA in Pi anions", "func" : noAAinAnions }   ]  )
    
    def row2ValuesAnionCation(rowId, row):
        return ( rowId, row["PDB Code"] , row["Cation Code"], 
                                                        row["Cation chain"]+str(row["Cation id"]) , row["Anion code"], 
                                                        row["Anion chain"] + str(row["Anion id"]), row["Anion symbol"], 
                                                        row["Cation symbol"], str(row["Distance"])[:3] )
        
    guiAnionCation.setRow2Values(row2ValuesAnionCation)
    
    def getSelectionAnionCation( data):
        res1Id = data["Anion id"].values[0]
        res1Chain = data["Anion chain"].values[0]
        
        res2Id = data["Cation id"].values[0]
        res2Chain = data["Cation chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    guiAnionCation.setSelectionFunc( getSelectionAnionCation )
    
    def getArrowAnionCation( data ):
        point1Coords = [ data["Anion x coord"].values[0] , data["Anion y coord"].values[0] , data["Anion z coord"].values[0] ]
        point2Coords = [ data["Cation x coord"].values[0] , data["Cation y coord"].values[0] , data["Cation z coord"].values[0] ]
        
        return point1Coords, point2Coords
    
    guiAnionCation.setArrowFunc(getArrowAnionCation)
    
    ######################
    # ALL
    ######################
    guis = [ guiAnionPi, guiPiPi, guiCationPi, guiAnionCation]
    
    for gui in guis:
        gui.grid()
    
    if simulation:
        self.mainloop()
   
if __name__ == "__main__":
    fetchdialog(True)