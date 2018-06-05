#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 12:18:08 2018

@author: michal
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

def __init_plugin__(self=None):
    plugins.addmenuitem('Supramolecular analyser', fetchdialog)


def fetchdialog(simulation = False):
    if simulation:
        root = Tkinter.Tk()
    else:
        app = plugins.get_pmgapp()
        root = plugins.get_tk_root()
        
    
    self = Tkinter.Toplevel(root)
    self.title('Supramolecular analyser')
    self.minsize(500, 500)
    self.resizable(0,0)
    
    logData = {"logFile": False, "data" : None, "anionCodesUnique" : [], "piAcidCodesUnique" : [], 
               "cifDir" : None, "arrowExists" : False, "displaying" : False, "displayingAround" : False }
    
    checkboxVars = {}
    
    def getLogFile():
        logData["logFile"] = tkFileDialog.askopenfilename(title = "Select file", filetypes = (("Log files","*.log"), ("Txt files", "*.txt") ,("CSV files","*.csv"), ("Dat files", "*.dat"), ("all files","*.*")) )
        ent_log.configure(state = "normal")
        ent_log.delete(0,"end")
        ent_log.insert(0, logData["logFile"])
        ent_log.configure(state = "readonly")
        try:
            logData["data"] = pd.read_csv(logData["logFile"], sep = "\t")
            updateMenu()
        except:
            tkMessageBox("Pandas cannot parse this file")
            
    def updateMenu():
        for parameter in listParameters:
            parametersData = logData["data"][ listParameters[parameter]["header"]  ].unique()
            parametersData = sorted(parametersData)
            logData[parameter] = parametersData
            
            listParameters[parameter]["listbox"].delete(0, "end")
            for row in parametersData:
                listParameters[parameter]["listbox"].insert("end", row)
        
    but_log = Tkinter.Button(self, text = "Load log file", command = getLogFile, width = 10)
    but_log.grid(row = 0, column = 0)
    
    ent_log = Tkinter.Entry(self, width =20)
    ent_log.configure(state = "readonly")
    ent_log.grid(row = 0, column = 1, columnspan = 2)
    
    lab_use = Tkinter.Label(self, text = "Use filter")
    lab_use.grid(row = 0, column = 3)
    
    numericalParameters = { "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" },
                           "alpha" : { "header" : "Angle" }, "res" : { "header" : "Resolution" }  }
    
    actualRow = 1
    
    for parameter in sorted(numericalParameters.keys()):
        numericalParameters[parameter]["entry_low"] = Tkinter.Entry(self, width = 10)
        numericalParameters[parameter]["entry_low"].grid(row = actualRow, column = 0)
        
        numericalParameters[parameter]["label"] = Tkinter.Label(self, width = 10, text = "< "+parameter+" <")
        numericalParameters[parameter]["label"].grid(row=actualRow, column =1)
        
        numericalParameters[parameter]["entry_high"] = Tkinter.Entry(self, width = 10)
        numericalParameters[parameter]["entry_high"].grid(row = actualRow, column = 2)
        
        checkboxVars[parameter] = Tkinter.IntVar()
    
        numericalParameters[parameter]["checkbox"] = Tkinter.Checkbutton(self, variable = checkboxVars[parameter] )
        numericalParameters[parameter]["checkbox"].grid(row=actualRow, column =  3)
        
        actualRow+=1
        
    def showAround():
        if not logData["displaying"]:
            return
        
        if chkvar_around.get() > 0 and not logData["displayingAround"]:
            selectionAroundName = "suprAround"
            cmd.select( selectionAroundName,  "byres ( suprSelection around 5 ) " )
            cmd.show("lines", selectionAroundName)
            logData["displayingAround"] = True
        elif chkvar_around.get() == 0 and logData["displayingAround"]:
            cmd.hide("lines" , "suprAround")
            cmd.delete("suprAround")
            logData["displayingAround"] = False
            
        cmd.deselect()
        
    
    lab_around = Tkinter.Label(self, text = "Around")
    lab_around.grid(row = actualRow, column = 0)
    
    chkvar_around = Tkinter.IntVar()
    
    chk_around = Tkinter.Checkbutton(self, variable = chkvar_around, command = showAround)
    chk_around.grid(row = actualRow, column = 1)
    
    actualRow += 1
    
    lab_noAAinPiAcids = Tkinter.Label(self, text = "No AA in Pi acids")
    lab_noAAinPiAcids.grid(row = actualRow, column = 0)
    
    chkvar_noAAinPiAcids = Tkinter.IntVar()
    
    chk_noAAinPiAcids = Tkinter.Checkbutton(self, variable = chkvar_noAAinPiAcids)
    chk_noAAinPiAcids.grid(row = actualRow, column = 1)
    
    actualRow += 1
    
    lab_noAAinAnions = Tkinter.Label(self, text = "No AA in anions")
    lab_noAAinAnions.grid(row = actualRow, column = 0)
    
    chkvar_noAAinAnions = Tkinter.IntVar()
    
    chk_noAAinAnions = Tkinter.Checkbutton(self, variable = chkvar_noAAinAnions)
    chk_noAAinAnions.grid(row = actualRow, column = 1)

    def listFilter(key):
        template = listParameters[key]["entry"].get()
        template = template.upper()
        
        templateLen = len(template)
        listParameters[key]["listbox"].delete(0, "end")
        
        for row in logData[key]:
            if template == row[:templateLen]:
                listParameters[key]["listbox"].insert("end", row)
                
    def piAcidFilter():
        listFilter("Pi acid")
        
    def anionFilter():
        listFilter("Anions")
        
    def anionTypeFilter():
        listFilter("Groups")
        
    def metalFilter():
        listFilter("Metals")

    listParameters = { "Pi acid" : { "header" : "Pi acid Code", "filterFunc" : piAcidFilter }, 
                      "Anions" : { "header" : "Anion code", "filterFunc" : anionFilter } , 
                      "Groups" : { "header" : "Anion type", "filterFunc" : anionTypeFilter },
                      "Metals" : { "header" : "Metal cations", "filterFunc" : metalFilter }} 
                
    actualColumn = 4
    for parameter in sorted(listParameters.keys()):
        listParameters[parameter]["label"] = Tkinter.Label(self, text = parameter)
        listParameters[parameter]["label"].grid(row = 0, column = actualColumn)
        
        checkboxVars[parameter] = Tkinter.IntVar()
        
        listParameters[parameter]["checkbox"] = Tkinter.Checkbutton(self, variable = checkboxVars[parameter] )
        listParameters[parameter]["checkbox"].grid(row = 0, column = actualColumn+1)
        
        listParameters[parameter]["listbox"] = Tkinter.Listbox(self, width =10, height = 15, exportselection = False)
        listParameters[parameter]["listbox"].grid(row = 1, column = actualColumn, rowspan = 15, columnspan = 2)
        
        listParameters[parameter]["entry"] = Tkinter.Entry(self, width = 5)
        listParameters[parameter]["entry"].grid(row =20, column = actualColumn)
        
        listParameters[parameter]["button"] = Tkinter.Button(self, width = 1, text = "*", command = listParameters[parameter]["filterFunc"])
        listParameters[parameter]["button"].grid(row = 20, column = actualColumn + 1)
        
        actualColumn += 2
        
    #SORTING
    sorting_keys2header = {  "R" : "Distance", "Angle" : "Angle", "x" : "x" , "h" : "h",
                         "res" : "Resolution", "Pi acid" : "Pi acid Code" ,"Anion" : "Anion code" }
    sorting_keys_col1 = [  "R" , "Angle" , "x" , "h" ,
                         "res" , "Pi acid" ,"Anion" ]
    sorting_keys_col2 = [ "Ascd", "Desc" ]
    sortingMenu = []
    
    for i in range(4):
        sortingMenu.append( { }  )
        
        sortingMenu[i]["label"] = Tkinter.Label(self, text = "Sorting"+str(i))
        sortingMenu[i]["label"].grid(row = 0, column =  actualColumn)
        
        sortingMenu[i]["chk_value"] = Tkinter.IntVar()
        
        sortingMenu[i]["chk_butt"] = Tkinter.Checkbutton(self, variable = sortingMenu[i]["chk_value"])
        sortingMenu[i]["chk_butt"] .grid(row = 0, column = actualColumn+1)
        
        sortingMenu[i]["sorting_key"] = Tkinter.IntVar()
        
        actual_row = 1
        for value, key in enumerate(sorting_keys_col1):
            sortingMenu[i][key] = Tkinter.Radiobutton(self, text = key, variable = sortingMenu[i]["sorting_key"], value = value, indicatoron=0, width = 8 )
            sortingMenu[i][key].grid(row = actual_row, column = actualColumn, columnspan = 2)
            actual_row += 1
            
        sortingMenu[i]["sortingTypeValue"] = Tkinter.IntVar()
        
        for value, key in enumerate(sorting_keys_col2):
            sortingMenu[i][key]= Tkinter.Radiobutton(self, text = key, variable = sortingMenu[i]["sortingTypeValue"], value = value, indicatoron=0, width = 8 )
            sortingMenu[i][key].grid(row = actual_row, column = actualColumn, columnspan = 2)
            actual_row += 1
            
        actualColumn += 2
        
    
    def applyFilter ():
        if logData["logFile"] == False:
            tkMessageBox.showwarning(title="Warning", message = "Log file not selected")
            
#        anythingSet = False
        actualData = logData["data"]
        for key in checkboxVars:
            if checkboxVars[key].get() > 0:
#                anythingSet = True
                
                if key in listParameters:
                    listIndex = listParameters[key]["listbox"].curselection()
                    if listIndex:
                        query = listParameters[key]["listbox"].get(listIndex)
                        actualData = actualData[  actualData[ listParameters[key]["header"] ] == query ]

                elif key in numericalParameters:
                    minValue = numericalParameters[key]["entry_low"].get()
                    maxValue = numericalParameters[key]["entry_high"].get()

                    try: 
                        minValue = float(minValue)
                        actualData = actualData[  actualData[ numericalParameters[key]["header"]  ] > minValue ]
                    except:
                        pass
                    
                    try:
                        maxValue = float(maxValue)
                        actualData = actualData[  actualData[ numericalParameters[key]["header"]  ] < maxValue ]
                        
                    except:
                        pass
                
#        if not anythingSet:
#            tkMessageBox.showwarning(title="Warning", message = "Please select any filter")
                        
        if chkvar_noAAinAnions.get() > 0:
            actualData = actualData[(actualData["Anion code"] != "GLU") & (actualData[ "Anion code" ] != "ASP" ) & ( actualData["Anion code"] != "TYR" ) & ( actualData["Anion code"] != "CYS") ]
        
        if chkvar_noAAinPiAcids.get() > 0:
            actualData = actualData[(actualData["Pi acid Code"] != "TYR") & (actualData[ "Pi acid Code" ] != "PHE" ) & ( actualData["Pi acid Code"] != "HIS" ) & ( actualData["Pi acid Code"] != "TRP" ) ]

        recordsFound = str(len(actualData))
        
        anySort = False
        columns = []
        ascending = []
        for sortData in sortingMenu:
           toSort =  sortData["chk_value"].get()
           
           if toSort > 0:
               anySort = True
               itemInd = sortData["sorting_key"].get()
               sortingKey = sorting_keys_col1[itemInd]
               header = sorting_keys2header[sortingKey]
               if header in columns:
                   continue
               columns.append(header)
               
               ascendingActual = sortData["sortingTypeValue"].get()
               if ascendingActual == 0:
                   ascending.append(True)
               else:
                   ascending.append(False)
               
        
        ent_recordsFound.configure(state = "normal")
        ent_recordsFound.delete(0,"end")
        ent_recordsFound.insert(0, str(recordsFound))
        ent_recordsFound.configure(state = "readonly")
        if anySort:
            actualData = actualData.sort_values(by = columns, ascending = ascending)
        logData["filtered"] = actualData
        tree_data.delete(*tree_data.get_children())
        
        rowId = 0
        for index, row in actualData.iterrows():
            tree_data.insert('', "end" , values =  ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                    row["Pi acid chain"]+str(row["Piacid id"]) , row["Anion code"], 
                                                    row["Anion chain"] + str(row["Anion id"]), row["Anion type"], 
                                                    str(row["Distance"])[:3], str(row["Angle"])[:4], str(row["x"])[:3],
                                                    str(row["h"])[:3], row["Resolution"] , row["Metal cations"]  ) )
            rowId += 1
            if rowId >= 1000:
                break
            
            
            
    but_apply = Tkinter.Button(self, width = 10, command = applyFilter, text = "Search")
    but_apply.grid(row = 22, column = 0)
    
    ent_recordsFound = Tkinter.Entry(self, width =20)
    ent_recordsFound.configure(state = "readonly")
    ent_recordsFound.grid(row = 25, column = 1, columnspan = 2)
    
# Records list
    lab_data = Tkinter.Label(self, width = 10, text = "Records found")
    lab_data.grid(row = 25, column = 0)
    
    headers = [ "ID" , "PDB" , "Pi acid", "Pi acid id", "Anion", "Anion id", "Anion type" , "R", "alpha", "x", "h", "res", "Metal" ]
    
    tree_data = ttk.Treeview(self, columns = headers, show = "headings", heigh = 20 )
    for header in headers:
        tree_data.heading(header, text = header)
        tree_data.column(header, width = 70)
    tree_data.grid(row = 30, column = 0, columnspan = 40)
    
    currentMolecule = { "PdbCode" : None  }
    
    def showInteractions():
        if not "filtered" in logData:
            return
        currentSel = tree_data.focus()
        if currentSel == "" :
            return
        
        rowId = tree_data.item(currentSel)["values"][0] 
        data = logData["filtered"].iloc[[rowId]]
        pdbCode = data["PDB Code"].values[0]
        
        
        if currentMolecule["PdbCode"] and currentMolecule["PdbCode"] != pdbCode:
            cmd.delete(currentMolecule["PdbCode"])
            
        if currentMolecule["PdbCode"] != pdbCode:
            if logData["cifDir"] != None:
                potentialPaths = [ path.join( logData["cifDir"] ,  pdbCode.lower() +".cif" ), path.join( logData["cifDir"] ,  pdbCode.upper() +".cif" )  ]
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
        
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Anion id"].values[0]
        res2Chain = data["Anion chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        
        selectionName = "suprSelection"
        
        cmd.select(selectionName, selection)
        cmd.show( "sticks" , selectionName  )
        cmd.center(selectionName)
        cmd.zoom(selectionName)
        
        if chkvar_around.get() > 0:
            selectionAroundName = "suprAround"
            cmd.select( selectionAroundName,  "byres ( suprSelection around 5 ) " )
            cmd.show("lines", selectionAroundName)
            logData["displayingAround"] = True
        else:
            logData["displayingAround"] = False
        
        centroidCoords = [ data["Centroid x coord"].values[0] , data["Centroid y coord"].values[0] , data["Centroid z coord"].values[0] ]
        anionAtomCoords = [ data["Anion x coord"].values[0] , data["Anion y coord"].values[0] , data["Anion z coord"].values[0] ]
        
        if logData["arrowExists"]:
            cmd.delete("Anion2Centroid")
            
        cgo_arrow(anionAtomCoords, centroidCoords, 0.1, name = "Anion2Centroid")
        logData["arrowExists"] = True
        logData["displaying"] = True
        currentMolecule["PdbCode"] = pdbCode
        
        cmd.deselect()
        
    
    but_showInteraction = Tkinter.Button(self, width = 10, command = showInteractions, text = "Show interact")
    but_showInteraction.grid(row = 50, column =0)
    
    def saveFiltered():
        if not "filtered" in logData:
            return
        
        file2save = tkFileDialog.asksaveasfilename(defaultextension = ".log", filetypes = (("Log files","*.log"), ("Txt files", "*.txt") ,("CSV files","*.csv"), ("Dat files", "*.dat"), ("all files","*.*")) )
        if file2save:
            logData["filtered"].to_csv(file2save, sep = "\t")
    
    but_saveFiltered = Tkinter.Button(self, width = 10, command = saveFiltered, text = "Save filtered")
    but_saveFiltered.grid(row = 50, column = 1)
    
    def selectCif():
        logData["cifDir"] = tkFileDialog.askdirectory()
        ent_cifDir.configure(state = "normal")
        ent_cifDir.delete(0,"end")
        ent_cifDir.insert(0, logData["cifDir"])
        ent_cifDir.configure(state = "readonly")
    
    but_cifDir = Tkinter.Button(self, width = 10, command = selectCif, text = "Cif dir")
    but_cifDir.grid(row = 50, column = 2)
    
    ent_cifDir = Tkinter.Entry(self, width =17)
    ent_cifDir.configure(state = "readonly")
    ent_cifDir.grid(row = 50, column = 3, columnspan = 3)
    
#    def printInfo():
#        tkMessageBox.showinfo(title="Important", message = "Emilia Kuźniak jest najpiękniejszą kobietą na świecie co najmniej od czasów Kleopatry (nie zachowały się żadne podobizny Kleopatry pozwalają na porównanie)")
#    
#    but_info = Tkinter.Button(self, width = 10, command = printInfo , text = "Important")
#    but_info.grid(row = 50, column = 6)
    
    if simulation:
        self.mainloop()
   
if __name__ == "__main__":
    fetchdialog(True)