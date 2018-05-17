#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 12:18:08 2018

@author: michal
"""

import sys
import pandas as pd
if sys.version_info[0] < 3:
    import Tkinter
    from Tkinter import LEFT, RIGHT
    import tkMessageBox, tkFileDialog
    from pymol import cmd, plugins, CmdException
    import ttk
    from pymol import cmd
else:
    import tkinter as Tkinter
    from tkinter import LEFT, RIGHT
    from tkinter import filedialog as tkFileDialog
    from tkinter import messagebox as tkMessageBox
    import tkinter.ttk as ttk



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
    
    logData = {"logFile": False, "data" : None, "anionCodesUnique" : [], "piAcidCodesUnique" : []}
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
    
    numericalParameters = { "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" }, "alpha" : { "header" : "Angle" }, "res" : { "header" : "Resolution" }  }
    
    actualRow = 1
    
    for parameter in numericalParameters:
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

    listParameters = { "Pi acid" : { "header" : "Pi acid Code", "filterFunc" : piAcidFilter }, 
                      "Anions" : { "header" : "Anion code", "filterFunc" : anionFilter } , 
                      "Groups" : { "header" : "Anion type", "filterFunc" : anionTypeFilter }  }
                
    actualColumn = 4
    for parameter in listParameters:
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
            
        anythingSet = False
        actualData = logData["data"]
        for key in checkboxVars:
            if checkboxVars[key].get() > 0:
                anythingSet = True
                
                if key in listParameters:
                    listIndex = listParameters[key]["listbox"].curselection()
                    if listIndex:
                        query = listParameters[key]["listbox"].get(listIndex)
                        actualData = actualData[  actualData[ listParameters[key]["header"] ].str.match(query) ]

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
                
        if not anythingSet:
            tkMessageBox.showwarning(title="Warning", message = "Please select any filter")

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
                                                    str(row["h"])[:3], row["Resolution"]) )
            rowId += 1
            if rowId >= 500:
                break
            
            
            
    but_apply = Tkinter.Button(self, width = 10, command = applyFilter, text = "Search")
    but_apply.grid(row = 22, column = 0)
    
    ent_recordsFound = Tkinter.Entry(self, width =20)
    ent_recordsFound.configure(state = "readonly")
    ent_recordsFound.grid(row = 25, column = 1, columnspan = 2)
    
# Records list
    lab_data = Tkinter.Label(self, width = 10, text = "Records found")
    lab_data.grid(row = 25, column = 0)
    
    headers = [ "ID" , "PDB" , "Pi acid", "Pi acid id", "Anion", "Anion id", "Anion type" , "R", "alpha", "x", "h", "res" ]
    
    tree_data = ttk.Treeview(self, columns = headers, show = "headings" )
    for header in headers:
        tree_data.heading(header, text = header)
        tree_data.column(header, width = 70)
    tree_data.grid(row = 30, column = 0, columnspan = 40)
    
    
    def showInteractions():
        if not "filtered" in logData:
            return
        currentSel = tree_data.focus()
        rowId = tree_data.item(currentSel)["values"][0] 
        data = logData["filtered"].iloc[[rowId]]
        pdbCode = data["PDB Code"].values[0]
        cmd.fetch(pdbCode)
        
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Anion id"].values[0]
        res2Chain = data["Anion chain"].values[0]
        cmd.select("resi " + str(res1Id) + "in chain "+res1Chain +" , resi "+str(res2Id)+" in chain " + res2Chain   )
        
    
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
    
    if simulation:
        self.mainloop()
   
if __name__ == "__main__":
    fetchdialog(True)