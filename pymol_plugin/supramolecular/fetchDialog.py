#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 19:40:59 2018

@author: michal
"""
import sys

try:
    from pymol import plugins
except:
    pass

if sys.version_info[0] < 3:
    import Tkinter
    import ttk
    import tkFileDialog
else:
    import tkinter as Tkinter
    from tkinter import filedialog as tkFileDialog
    import tkinter.ttk as ttk

from supramolecularComposition import SupramolecularComposition

def fetchdialog(simulation = False):
    if simulation:
        root = Tkinter.Tk()
    else:
        app = plugins.get_pmgapp()
        root = plugins.get_tk_root()
        
    appData = { "cifDir" : False }
    
    self = Tkinter.Toplevel(root)
    self.title('Supramolecular analyser')
    self.minsize(1320, 800)
    self.resizable(1,1)
    
    nb = ttk.Notebook(self, height = 700, width = 1320)

    pageAnionPi = ttk.Frame(nb)
    pagePiPi = ttk.Frame(nb)
    pageCationPi = ttk.Frame(nb)
    pageAnionCation  = ttk.Frame(nb)
    pageHBonds = ttk.Frame(nb)
    pageMetalLigand = ttk.Frame(nb)
    pageLinearAnionPi = ttk.Frame(nb)
    pagePlanarAnionPi = ttk.Frame(nb)
    pageQM = ttk.Frame(nb)
    pageJobStatus = ttk.Frame(nb)
    
    nb.add(pageAnionPi, text = "AnionPi")
    nb.add(pagePiPi, text = "PiPi")
    nb.add(pageCationPi, text = "CationPi")
    nb.add(pageAnionCation, text = "AnionCation")
    nb.add(pageHBonds, text = "HBonds")
    nb.add(pageMetalLigand, text = "MetalLigand")
    nb.add(pageLinearAnionPi, text = "LinearAnionPi")
    nb.add(pagePlanarAnionPi, text = "PlanarAnionPi")
    nb.add(pageQM, text = "QM assistant")
    nb.add(pageJobStatus, text = "Job monitor")
    
    nb.grid(column = 0, row = 0, columnspan = 20)
    
    supramolecularComposition = SupramolecularComposition(pageAnionPi, pagePiPi,
                                                          pageCationPi, pageAnionCation, pageHBonds, pageMetalLigand,
                                                          pageLinearAnionPi, pagePlanarAnionPi ,pageQM, pageJobStatus)
    
    
    ######################
    # GENERAL
    ######################
    
    def readLogDir():
        appData["logDir"] = tkFileDialog.askdirectory()
        if not appData["logDir"]:
            print("nihuhu")
            return
        
        ent_logDir.configure(state = "normal")
        ent_logDir.delete(0,"end")
        ent_logDir.insert(0, appData["logDir"])
        ent_logDir.configure(state = "readonly")
        
        supramolecularComposition.readAllLogsFromDir(appData["logDir"])
        
    
    but_readLogDir = Tkinter.Button(self, width = 20, text = "Read log Dir", command = readLogDir)
    but_readLogDir.grid(row = 2, column = 0, columnspan = 2)
    
    ent_logDir = Tkinter.Entry(self, width =45)
    ent_logDir.configure(state = "readonly")
    ent_logDir.grid(row = 2, column = 2, columnspan = 4)
    
    def selectCif():
        appData["cifDir"] = tkFileDialog.askdirectory()
        if not appData["cifDir"]:
            print("nihuhu")
            return
        
        ent_cifDir.configure(state = "normal")
        ent_cifDir.delete(0,"end")
        ent_cifDir.insert(0, appData["cifDir"])
        ent_cifDir.configure(state = "readonly")
        
        supramolecularComposition.selectCifDir(appData["cifDir"])
        
        
    
    but_cifDir = Tkinter.Button(self, width = 10, command = selectCif, text = "Cif dir")
    but_cifDir.grid(row = 2, column = 6)
    
    ent_cifDir = Tkinter.Entry(self, width =45)
    ent_cifDir.configure(state = "readonly")
    ent_cifDir.grid(row = 2, column = 7, columnspan = 3)
    
    def setParellelSelection():
        state = var_parallelSelection.get()
        supramolecularComposition.setParallelSelection(state)
    
    lab_parallelSelection = Tkinter.Label(self, text = "Parallel selection")
    lab_parallelSelection.grid(row = 2, column = 10)
    
    var_parallelSelection = Tkinter.IntVar()
    chk_parallelSelection = Tkinter.Checkbutton(self, variable = var_parallelSelection ,command = setParellelSelection)
    chk_parallelSelection.grid(row =2, column = 11)
    
    actionMenu = {}
    
    column = 1
    
    lab_usePage = Tkinter.Label(self, width = 10, text = "Use:")
    lab_usePage.grid(row = 4, column = 0)
    
    lab_doNotUsePage = Tkinter.Label(self, width = 10, text = "Exclude:" )
    lab_doNotUsePage.grid(row = 5, column = 0 )
    
    for label in supramolecularComposition.actionLabels:
        actionMenu[label] = {}
        
        actionMenu[label]["label"] = Tkinter.Label(self, text = label)
        actionMenu[label]["label"].grid(row = 3, column = column )
        
        actual_row = 4
        
        actionMenu[label]["checkValue"] = Tkinter.IntVar()
        actionMenu[label]["checkbox"] = Tkinter.Checkbutton(self, variable = actionMenu[label]["checkValue"])
        actionMenu[label]["checkbox"].grid(row = actual_row, column = column)
        
        actual_row += 1
        
        actionMenu[label]["checkValueExclude"] = Tkinter.IntVar()
        actionMenu[label]["checkboxExclude"] = Tkinter.Checkbutton(self, variable = actionMenu[label]["checkValueExclude"])
        actionMenu[label]["checkboxExclude"].grid(row = actual_row, column = column)
        
        
            
        column += 1
    
    lab_showInt = Tkinter.Label(self, width = 10, text = "Show:" )
    lab_showInt.grid(row = 6, column = 0)
    
    showMenu = {}
    column = 1
    
    for label in supramolecularComposition.actionLabels:
        showMenu[label] = {}
        
        showMenu[label]["checkValue"] = Tkinter.IntVar()
        showMenu[label]["checkbox"] = Tkinter.Checkbutton(self, variable = showMenu[label]["checkValue"])
        showMenu[label]["checkbox"].grid(row = 6, column = column)
            
        column += 1
    
    
    def mergeResults():
        supramolecularComposition.merge(actionMenu)
            
    
    but_merge = Tkinter.Button(self, width = 20, text = "Merge!", command = mergeResults)
    but_merge.grid(row = 5, column = 9, columnspan = 2)
    
    def showAllInteractions():
        supramolecularComposition.showAll(showMenu)
    
    but_showMany = Tkinter.Button(self, width = 20, text = "Show", command = showAllInteractions)
    but_showMany.grid(row = 6, column = 9, columnspan = 2)
    
    but_saveState = Tkinter.Button(self, width = 20, text = "Save GUI state", command = supramolecularComposition.saveState)
    but_saveState.grid(row = 5, column = 11, columnspan = 2)
    
    but_loadState = Tkinter.Button(self, width = 20, text = "Load GUI state", command = supramolecularComposition.loadState)
    but_loadState.grid(row = 6, column = 11, columnspan = 2)
    
    ######################
    # ALL
    ######################
    
    supramolecularComposition.grid()
    
    if simulation:
        self.mainloop()