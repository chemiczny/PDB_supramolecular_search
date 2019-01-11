#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 19:40:59 2018

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
    self.minsize(1500, 800)
    self.resizable(0,0)
    
    nb = ttk.Notebook(self, height = 700, width = 1400)

    pageAnionPi = ttk.Frame(nb)
    pagePiPi = ttk.Frame(nb)
    pageCationPi = ttk.Frame(nb)
    pageAnionCation  = ttk.Frame(nb)
    pageHBonds = ttk.Frame(nb)
    pageMetalLigand = ttk.Frame(nb)
    
    nb.add(pageAnionPi, text = "AnionPi")
    nb.add(pagePiPi, text = "PiPi")
    nb.add(pageCationPi, text = "CationPi")
    nb.add(pageAnionCation, text = "AnionCation")
    nb.add(pageHBonds, text = "HBonds")
    nb.add(pageMetalLigand, text = "MetalLigand")
    
    nb.grid(column = 0, row = 0, columnspan = 20)
    
    supramolecularComposition = SupramolecularComposition(pageAnionPi, pagePiPi,
                                                          pageCationPi, pageAnionCation, pageHBonds, pageMetalLigand)
    
    
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
    but_readLogDir.grid(row = 2, column = 1, columnspan = 1)
    
    ent_logDir = Tkinter.Entry(self, width =45)
    ent_logDir.configure(state = "readonly")
    ent_logDir.grid(row = 2, column = 2, columnspan = 2)
    
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
    but_cifDir.grid(row = 2, column = 4)
    
    ent_cifDir = Tkinter.Entry(self, width =45)
    ent_cifDir.configure(state = "readonly")
    ent_cifDir.grid(row = 2, column = 5, columnspan = 3)
    
    actionMenu = {}
    
    column = 2
    
    lab_usePage = Tkinter.Label(self, width = 10, text = "Use:")
    lab_usePage.grid(row = 4, column = 1)
    
    lab_doNotUsePage = Tkinter.Label(self, width = 10, text = "Exclude:" )
    lab_doNotUsePage.grid(row = 5, column = 1 )
    
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
    lab_showInt.grid(row = 6, column = 1)
    
    showMenu = {}
    column = 2
    
    for label in supramolecularComposition.actionLabels:
        showMenu[label] = {}
        
        showMenu[label]["checkValue"] = Tkinter.IntVar()
        showMenu[label]["checkbox"] = Tkinter.Checkbutton(self, variable = showMenu[label]["checkValue"])
        showMenu[label]["checkbox"].grid(row = 6, column = column)
            
        column += 1
    
    
    def mergeResults():
        supramolecularComposition.merge(actionMenu)
            
        
#    def countSubstringInList( list2check, substring ):
#        found = 0
#        
#        for el in list2check:
#            if substring in list2check:
#                found +=1
#                
#        return found
    
    but_merge = Tkinter.Button(self, width = 20, text = "Merge!", command = mergeResults)
    but_merge.grid(row = 5, column = 8, columnspan = 2)
    
    def showAllInteractions():
        supramolecularComposition.showAll(showMenu)
    
    but_showMany = Tkinter.Button(self, width = 20, text = "Show", command = showAllInteractions)
    but_showMany.grid(row = 6, column = 8, columnspan = 2)
    
    ######################
    # ALL
    ######################
    
    supramolecularComposition.grid()
    
    if simulation:
        self.mainloop()