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
        piAcidCodes = logData["data"]["Pi acid Code"].unique()
        piAcidCodes = sorted(piAcidCodes)
        logData["piAcidCodesUnique"] = piAcidCodes
        list_piAcid.delete(0, "end")
        for piAcid in piAcidCodes:
            list_piAcid.insert("end", piAcid)
            
        anionCodes = logData["data"]["Anion code"].unique()
        anionCodes = sorted(anionCodes)
        logData["anionCodesUnique"] = anionCodes
        list_anions.delete(0, "end")
        for anion in anionCodes:
            list_anions.insert("end", anion)
            
        anionTypes = logData["data"]["Anion type"].unique()
        anionTypes = sorted(anionTypes)
        logData["anionTypesUnique"] = anionCodes
        list_anionsTypes.delete(0, "end")
        for anion in anionTypes:
            list_anionsTypes.insert("end", anion)
        
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

    
# Pi acids list
    def piAcidFilter():
        piAcidTemplate = ent_piAcid.get()
        piAcidTemplate = piAcidTemplate.upper()
        
        templateLen = len(piAcidTemplate)
        list_piAcid.delete(0, "end")
#        print(templateLen)
        for piAcid in logData["piAcidCodesUnique"]:
            if piAcidTemplate == piAcid[:templateLen]:
                list_piAcid.insert("end", piAcid)
    
    lab_piAcid = Tkinter.Label(self, text = "Pi acid")
    lab_piAcid.grid(row = 0, column = 4)
    
    checkboxVars["piAcid"] = Tkinter.IntVar()
    
    check_piAcid = Tkinter.Checkbutton(self, variable = checkboxVars["piAcid"] )
    check_piAcid.grid(row = 0, column = 5)
    
    list_piAcid = Tkinter.Listbox(self, width =10, height = 10, exportselection = False)
    list_piAcid.grid( row = 1, column = 4, rowspan = 10, columnspan = 2)
    
    ent_piAcid = Tkinter.Entry(self, width = 5)
    ent_piAcid.grid(row = 15, column = 4)
    
    but_piAcid = Tkinter.Button(self, width = 1, text = "*", command = piAcidFilter)
    but_piAcid.grid(row = 15, column = 5)
    
# Anions list
    def anionFilter():
        anionTemplate = ent_anions.get()
        anionTemplate = anionTemplate.upper()
        
        list_anions.delete(0, "end")
        templateLen = len(anionTemplate)
        for anion in logData["anionCodesUnique"]:
            if anionTemplate == anion[:templateLen]:
                list_anions.insert("end", anion)
        
    
    lab_anions = Tkinter.Label(self, text = "Anions")
    lab_anions.grid(row = 0, column = 6)
    
    checkboxVars["anion"] = Tkinter.IntVar()
    
    check_anions = Tkinter.Checkbutton(self, variable = checkboxVars["anion"])
    check_anions.grid(row = 0, column = 7)
    
    list_anions = Tkinter.Listbox(self, width =10, height = 10, exportselection = False)
    list_anions.grid( row = 1, column = 6, rowspan = 10, columnspan = 2)
    
    ent_anions = Tkinter.Entry(self, width = 5)
    ent_anions.grid(row = 15, column = 6)
    
    but_anions = Tkinter.Button(self, width = 1, text = "*", command = anionFilter)
    but_anions.grid(row = 15, column = 7)
    
# Anions types list
    
    def anionTypeFilter():
        anionTemplate = ent_anionsTypes.get()
        anionTemplate = anionTemplate.upper()
        
        list_anionsTypes.delete(0, "end")
        templateLen = len(anionTemplate)
        for anion in logData["anionTypesUnique"]:
            if anionTemplate == anion[:templateLen]:
                list_anionsTypes.insert("end", anion)
        
    
    lab_anionsTypes = Tkinter.Label(self, text = "Groups")
    lab_anionsTypes.grid(row = 0, column = 8)
    
    checkboxVars["anionType"] = Tkinter.IntVar()
    
    check_anionsTypes = Tkinter.Checkbutton(self, variable = checkboxVars["anionType"])
    check_anionsTypes.grid(row = 0, column = 9)
    
    list_anionsTypes = Tkinter.Listbox(self, width =10, height = 10, exportselection = False)
    list_anionsTypes.grid( row = 1, column = 8, rowspan = 10, columnspan = 2)
    
    ent_anionsTypes = Tkinter.Entry(self, width = 5)
    ent_anionsTypes.grid(row = 15, column = 8)
    
    but_anionsTypes = Tkinter.Button(self, width = 1, text = "*", command = anionFilter)
    but_anionsTypes.grid(row = 15, column = 9)
    
# Anions types end
    
    def applyFilter ():
        if logData["logFile"] == False:
            tkMessageBox.showwarning(title="Warning", message = "Log file not selected")
            
        anythingSet = False
        actualData = logData["data"]
        for key in checkboxVars:
            if checkboxVars[key].get() > 0:
#                print(key, "jest wybrany!")
                anythingSet = True
                
                if key == "anion":
                    anionCode = list_anions.get(list_anions.curselection())
                    actualData = actualData[  actualData["Anion code"].str.match(anionCode) ]
#                    print(anionCode)
                elif key == "piAcid":
                    piAcidCode = list_piAcid.get(list_piAcid.curselection())
                    actualData = actualData[  actualData["Pi acid Code"].str.match(piAcidCode) ]
#                    print(piAcidCode)
                elif key == "anionType":
                    anionType = list_anionsTypes.get(list_anionsTypes.curselection())
                    actualData = actualData[  actualData["Anion type"].str.match(anionType) ]
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
        else:
            recordsFound = str(len(actualData))
            
            ent_recordsFound.configure(state = "normal")
            ent_recordsFound.delete(0,"end")
            ent_recordsFound.insert(0, str(recordsFound))
            ent_recordsFound.configure(state = "readonly")
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
    but_apply.grid(row = 7, column = 0)
    
    ent_recordsFound = Tkinter.Entry(self, width =20)
    ent_recordsFound.configure(state = "readonly")
    ent_recordsFound.grid(row = 15, column = 1, columnspan = 2)
    
# Records list
    lab_data = Tkinter.Label(self, width = 10, text = "Records found")
    lab_data.grid(row = 15, column = 0)
    
    headers = [ "ID" , "PDB" , "Pi acid", "Pi acid id", "Anion", "Anion id", "Anion type" , "R", "alpha", "x", "h", "res" ]
    
    tree_data = ttk.Treeview(self, columns = headers, show = "headings" )
    for header in headers:
        tree_data.heading(header, text = header)
        tree_data.column(header, width = 70)
    tree_data.grid(row = 18, column = 0, columnspan = 40)
    
    
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