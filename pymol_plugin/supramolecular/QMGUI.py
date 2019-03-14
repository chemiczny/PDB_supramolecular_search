#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:21:16 2019

@author: michal
"""
import sys
import os

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
        
#        seleAroundLabel = Tkinter.Label(self.page, text = "around sele")
#        seleAroundLabel.grid(row=4, column=0)
#        
#        self.seleRadius = Tkinter.Entry(self.page, width=5)
#        self.seleRadius.grid(row=5, column=0)
#        self.seleRadius.insert('end', "5")
        
        
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
        
#        self.hostAroundSele = Tkinter.Button(self.page, width = 5, text = "get", command = self.readHostFromSeleAround)
#        self.hostAroundSele.grid(row=4, column = 1)
        
    def readHost(self):
        try:
            stateNo = cmd.get_state()
            cmd.create("host", "sele", stateNo)
            cmd.select( "hostFrozen", "none")
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
        
#        self.guestAround = Tkinter.Button(self.page, text = "get", width = 5, command = self.readGuestFromSeleAround)
#        self.guestAround.grid(row=4, column=2)
        
    def readGuest(self):
        try:
            stateNo = cmd.get_state()
            cmd.create("guest", "sele", stateNo)
            cmd.select("guestFrozen", "none")
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
        writeButton = Tkinter.Button(self.page, text = "write", command = self.write)
        writeButton.grid(row=10, column = 0)
        
    def write(self):
        options = { "mustexist" : False , "title" : "Job directory selection"}
        directory = tkFileDialog.askdirectory(**options)
        if not directory:
            return
        
        if not os.path.isdir(directory):
            os.makedirs(directory)
            
        basename = os.path.basename(directory)
        
        slurmFile = os.path.join( directory, basename + ".slurm" )
        inpFile = os.path.join(directory, basename + ".inp")
        inpFileBasename = os.path.basename(inpFile)
        
        slurmF = open(slurmFile, 'w')
        
        slurmHead = self.slurmTextG16.get("1.0", "end")
        
        slurmF.write(slurmHead)
        slurmF.write("\n")
        
        slurmF.write("\nmodule add plgrid/apps/gaussian/g16.A.03\n\n")
        slurmF.write("g16 "+inpFileBasename+ "\n\n")
        
        slurmF.close()
        
        inpF = open(inpFile, 'w')
        
        inpF.write("%Chk="+inpFileBasename.replace(".inp",".chk")+"\n")
        routeSection = self.routeSection.get("1.0", "end")
        inpF.write(routeSection)
        
        inpF.write("\nEmilka jest najpiekniejsza!\n\n")
        
        chargeTotal = self.complexCharge.get()
        spinTotal = self.complexSpin.get()
        chargeGuest = self.guestCharge.get()
        spinGuest = self.guestSpin.get()
        chargeHost = self.hostCharge.get()
        spinHost = self.hostSpin.get()
        
        inpF.write(chargeTotal+","+spinTotal+" "+chargeGuest+","+spinGuest+" "+chargeHost+","+spinHost+"\n")
        model = cmd.get_model("guest and guestFrozen")
        
        for atom in model.atom:
            resnum = atom.resi
            resname = atom.resn
            pdbname = atom.name
            element = atom.symbol
            
            x = str(atom.coord[0])
            y = str(atom.coord[1])
            z = str(atom.coord[2])
            
            inpF.write(element+"(Fragment=1, PDBName="+pdbname+", ResName="+resname+", ResNum="+str(resnum)+") -1 "+x+" "+y+" "+z+"\n" )
            
            
        model = cmd.get_model("guest and not guestFrozen")
        
        for atom in model.atom:
            resnum = atom.resi
            resname = atom.resn
            pdbname = atom.name
            element = atom.symbol
            
            x = str(atom.coord[0])
            y = str(atom.coord[1])
            z = str(atom.coord[2])
            
            inpF.write(element+"(Fragment=1, PDBName="+pdbname+", ResName="+resname+", ResNum="+str(resnum)+") "+x+" "+y+" "+z+"\n" )
            
            
        model = cmd.get_model("host and hostFrozen")
        
        for atom in model.atom:
            resnum = atom.resi
            resname = atom.resn
            pdbname = atom.name
            element = atom.symbol
            
            x = str(atom.coord[0])
            y = str(atom.coord[1])
            z = str(atom.coord[2])
            
            inpF.write(element+"(Fragment=2, PDBName="+pdbname+", ResName="+resname+", ResNum="+str(resnum)+") -1 "+x+" "+y+" "+z+"\n" )
            
        model = cmd.get_model("host and not hostFrozen")
        
        for atom in model.atom:
            resnum = atom.resi
            resname = atom.resn
            pdbname = atom.name
            element = atom.symbol
            
            x = str(atom.coord[0])
            y = str(atom.coord[1])
            z = str(atom.coord[2])
            
            inpF.write(element+"(Fragment=2, PDBName="+pdbname+", ResName="+resname+", ResNum="+str(resnum)+") "+x+" "+y+" "+z+"\n" )
            
        inpF.write("\nstoichiometry=H2O1\n")
        inpF.write("solventname=Water2\n")
        inpF.write("eps=4\n")
        inpF.write("epsinf=1.77556\n")
        inpF.write("\n\n")
        
        inpF.close()
        
    def gridGaussianRouteSection(self):
        routeSectionLabel = Tkinter.Label(self.page, text = "Gaussian route section")
        routeSectionLabel.grid(row = 1, column=5, columnspan = 5)
        
        self.routeSection =Tkinter.Text(self.page, width =50, height = 10 )
        self.routeSection.grid(row = 2, column = 5, columnspan = 5, rowspan = 5)
        self.routeSection.insert("end", "%Mem=100GB\n#P B3LYP/6-31G(d,p)\n# Opt Counterpoise=2\n# SCRF(Solvent=Generic, Read)\n# Gfinput IOP(6/7=3)  Pop=full  Density  Test \n# Units(Ang,Deg)")
        
    def gridSlurmSection(self):
        slurmSection = Tkinter.Label(self.page, text = "Slurm input")
        slurmSection.grid(row = 1, column = 10, columnspan = 10)
        
        self.slurmTextG16 = Tkinter.Text(self.page, width = 50, height = 10)
        self.slurmTextG16.grid(row = 2, column = 10 , columnspan = 5, rowspan = 5)
        self.slurmTextG16.insert("end" , "#!/bin/env bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=24\n#SBATCH --time=70:00:00\n##### Nazwa kolejki\n#SBATCH -p plgrid\n" )
                  
    
    def grid(self):
        self.gridBasicLabels()
        self.gridHost()
        self.gridGuest()
        self.gridComplex()
        self.gridWrite()
        self.gridSlurmSection()
        
        self.gridGaussianRouteSection()