#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:21:16 2019

@author: michal
"""
import sys

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
    
def getAllSelectionNames():
    return cmd.get_names("selections", 0)
    
class QMGUI:
    def __init__(self, page):
        self.page = page
        
    def gridBasicLabels(self):
        chargeLabel = Tkinter.Label(self.page, text = "charge")
        chargeLabel.grid(row=1, column=0)
        
        spinLabel = Tkinter.Label(self.page, text = "spin")
        spinLabel.grid(row=2, column=0)
        
        seleLabel = Tkinter.Label(self.page, text = "sele")
        seleLabel.grid(row=3, column=0)
        
        seleAroundLabel = Tkinter.Label(self.page, text = "around sele")
        seleAroundLabel.grid(row=4, column=0)
        
        self.seleRadius = Tkinter.Entry(self.page, width=5)
        self.seleRadius.grid(row=5, column=0)
        self.seleRadius.insert('end', "5")
        
        
    def gridHost(self):
        hostLabel = Tkinter.Label(self.page, text = "Host")
        hostLabel.grid(row=0, column = 1)
        
        self.hostCharge = Tkinter.Entry(self.page, width = 5)
        self.hostCharge.grid(row=1, column=1)
        
        self.hostSpin = Tkinter.Entry(self.page, width = 5)
        self.hostSpin.grid(row=2, column =1)
        self.hostSpin.insert( 'end', "1")
        
        self.hostSele = Tkinter.Button(self.page, width = 5, text = "get", command = self.readHost)
        self.hostSele.grid(row=3, column =1)
        
        self.hostAroundSele = Tkinter.Button(self.page, width = 5, text = "get", command = self.readHostFromSeleAround)
        self.hostAroundSele.grid(row=4, column = 1)
        
    def readHost(self):
        try:
            stateNo = cmd.get_state()
            cmd.create("host", "sele", stateNo)
        except:
            print("lo kurla")
            
    def readHostFromSeleAround(self):
        try:
            stateNo = cmd.get_state()
            radius = self.seleRadius.get()
            cmd.create("host", "byres ( sele around "+radius+")", stateNo)
        except:
            print("lo kurla")
        
    def gridGuest(self):
        guestLabel = Tkinter.Label(self.page, text = "Guest")
        guestLabel.grid(row=0, column = 2)
        
        self.guestCharge = Tkinter.Entry(self.page, width = 5)
        self.guestCharge.grid(row=1, column=2)
        
        self.guestSpin = Tkinter.Entry(self.page, width = 5)
        self.guestSpin.grid(row=2, column =2)
        self.guestSpin.insert('end', '1')
        
        self.guestSele = Tkinter.Button(self.page, text = "get" ,  width = 5, command = self.readGuest)
        self.guestSele.grid(row=3, column =2)
        
        self.guestAround = Tkinter.Button(self.page, text = "get", width = 5, command = self.readGuestFromSeleAround)
        self.guestAround.grid(row=4, column=2)
        
    def readGuest(self):
        try:
            stateNo = cmd.get_state()
            cmd.create("guest", "sele", stateNo)
        except:
            print("lo kurla")
            
    def readGuestFromSeleAround(self):
        try:
            stateNo = cmd.get_state()
            radius = self.seleRadius.get()
#            print(radius)
            cmd.create("guest", "byres ( sele around "+radius+")", stateNo)
        except:
            print("lo kurla")
    
    def gridComplex(self):
        complexLabel = Tkinter.Label(self.page, text = "Complex")
        complexLabel.grid(row=0, column = 3)
        
        self.complexCharge = Tkinter.Entry(self.page, width = 5)
        self.complexCharge.grid(row=1, column=3)
        
        self.complexSpin = Tkinter.Entry(self.page, width = 5)
        self.complexSpin.grid(row=2, column =3)
        self.complexSpin.insert('end', '1')
        
    
    def gridWrite(self):
        writeButton = Tkinter.Button(self.page, text = "write")
        
    
    def grid(self):
        self.gridBasicLabels()
        self.gridHost()
        self.gridGuest()
        self.gridComplex()