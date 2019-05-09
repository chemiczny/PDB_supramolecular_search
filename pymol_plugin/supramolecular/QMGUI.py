#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:21:16 2019

@author: michal
"""
import sys
import os
from copy import deepcopy

try:
    from pymol import cmd
except:
    pass

if sys.version_info[0] < 3:
    import Tkinter
    import tkMessageBox, tkFileDialog
else:
    import tkinter as Tkinter
    from tkinter import filedialog as tkFileDialog
    from tkinter import messagebox as tkMessageBox
    
def getAllSelectionNames():
    return cmd.get_names("selections", 0)
    
class QMGUI:
    def __init__(self, page):
        self.page = page
        self.keywordSet = {}
        
        self.chargeTotal = None
        self.spinTotal = None
        self.chargeGuest = None
        self.spinGuest = None
        self.chargeHost = None
        self.spinHost = None
        
    def gridBasicLabels(self):
        actualRow = 1
        chargeLabel = Tkinter.Label(self.page, text = "charge")
        chargeLabel.grid(row=actualRow, column=0)
        
        actualRow+=1
        
        spinLabel = Tkinter.Label(self.page, text = "spin")
        spinLabel.grid(row=actualRow, column=0)
        
        actualRow +=1
        
        seleLabel = Tkinter.Label(self.page, text = "sele")
        seleLabel.grid(row=actualRow, column=0, rowspan=2)
        
        actualRow +=2
        
        frozenLabel = Tkinter.Label(self.page, text = "frozen")
        frozenLabel.grid(row=actualRow, column=0, rowspan=2)
        
        actualRow+=2
        
        addAtomsLabel = Tkinter.Label(self.page, text = "Add atoms")
        addAtomsLabel.grid(row=actualRow, column=0, rowspan=2)
        
        actualRow+=2
        
        self.countAtoms = Tkinter.Button(self.page, width=5, text = "count at", command = self.countAtoms)
        self.countAtoms.grid(row=actualRow, column=0)
        
#        seleAroundLabel = Tkinter.Label(self.page, text = "around sele")
#        seleAroundLabel.grid(row=4, column=0)
#        
#        self.seleRadius = Tkinter.Entry(self.page, width=5)
#        self.seleRadius.grid(row=5, column=0)
#        self.seleRadius.insert('end', "5")
    def countAtoms(self):
        try:
            modelHost = cmd.get_model("host")
            hostAtomNo =  len(modelHost.atom)
        except:
            hostAtomNo = 0
        
        try:
            modelGuest = cmd.get_model("guest")
            guestAtomNo = len(modelGuest.atom)
        except:
            guestAtomNo = 0
        
        complexAtomNo = hostAtomNo + guestAtomNo
        
        self.guestAtomsNo.delete(0, "end")
        self.guestAtomsNo.insert("end", str(guestAtomNo))
        
        self.hostAtomsNo.delete(0, "end")
        self.hostAtomsNo.insert("end", str(hostAtomNo))
        
        self.complexAtomsNo.delete(0, "end")
        self.complexAtomsNo.insert("end", str(complexAtomNo))


        
    def gridHost(self):
        actualRow = 0
        
        hostLabel = Tkinter.Label(self.page, text = "Host")
        hostLabel.grid(row=actualRow, column = 1)
        
        actualRow +=1
        
        self.hostCharge = Tkinter.Entry(self.page, width = 5)
        self.hostCharge.grid(row=actualRow, column=1)
        
        actualRow +=1
        
        self.hostSpin = Tkinter.Entry(self.page, width = 5)
        self.hostSpin.grid(row=actualRow, column =1)
        self.hostSpin.insert( 'end', "1")
        
        actualRow +=1
        
        self.hostSele = Tkinter.Button(self.page, width = 5, text = "get", command = self.readHost)
        self.hostSele.grid(row=actualRow, column =1)
        
        actualRow+=1
        
        self.hostAddSele = Tkinter.Button(self.page, width = 5, text = "add", command = self.addSele2Host)
        self.hostAddSele.grid(row=actualRow, column =1)
        
        actualRow+=1
        
        self.hostFrozenClear = Tkinter.Button(self.page, width =5, text = "clear", command = self.clearHostFrozen)
        self.hostFrozenClear.grid(row=actualRow, column=1)
        
        actualRow+=1
        
        self.hostFrozenDefault = Tkinter.Button(self.page, width = 5, text = "default", command = self.defaultHostFrozen)
        self.hostFrozenDefault.grid(row=actualRow, column=1)
        
        actualRow+=1
        
        self.hostAddH = Tkinter.Button(self.page, width = 5, text = "addH", command = self.addHHost)
        self.hostAddH.grid(row=actualRow, column=1)
        
        actualRow+=1
        
        self.hostAddNCA = Tkinter.Button(self.page, width = 5, text = "chop", command = self.hostChop)
        self.hostAddNCA.grid(row=actualRow, column=1)
        
        actualRow+=1
        
        self.hostAtomsNo = Tkinter.Entry(self.page, width = 5)
        self.hostAtomsNo.grid(row=actualRow, column=1)
        
        actualRow+=1
        
#        self.hostAroundSele = Tkinter.Button(self.page, width = 5, text = "get", command = self.readHostFromSeleAround)
#        self.hostAroundSele.grid(row=4, column = 1)
        
    def addSele2Host(self):
        try:
            stateNo = cmd.get_state()
            cmd.create("host", "%host or sele", stateNo)
        except:
            stateNo = cmd.get_state()
            cmd.create("host", "sele", stateNo)
        
    def clearHostFrozen(self):
        cmd.select( "hostFrozen", "none")
        
    def defaultHostFrozen(self):
        cmd.select("hostFrozen", " %hostFrozen or ( host and ( name CA ) ) ")
        
    def addHHost(self):
        cmd.h_add("host")
        
    def hostChop(self):
        self.chop("host")
        
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
        actualRow = 0
        guestLabel = Tkinter.Label(self.page, text = "Guest")
        guestLabel.grid(row=actualRow, column = 2)
        
        actualRow+=1
        
        self.guestCharge = Tkinter.Entry(self.page, width = 5)
        self.guestCharge.grid(row=actualRow, column=2)
        
        actualRow+=1
        
        self.guestSpin = Tkinter.Entry(self.page, width = 5)
        self.guestSpin.grid(row=actualRow, column =2)
        self.guestSpin.insert('end', '1')
        
        actualRow+=1
        
        self.guestSele = Tkinter.Button(self.page, text = "get" ,  width = 5, command = self.readGuest)
        self.guestSele.grid(row=actualRow, column =2)
        
        actualRow+=1
        
        self.guestSeleAdd = Tkinter.Button(self.page, text = "add", width = 5, command = self.addSele2Guest)
        self.guestSeleAdd.grid(row=actualRow, column=2)
        
        actualRow+=1
        
        self.guestFrozenClear = Tkinter.Button(self.page, width = 5, text = "clear", command = self.clearGuestFrozen)
        self.guestFrozenClear.grid(row=actualRow, column=2)
        
        actualRow+=1
        
        self.guestFrozenDefault = Tkinter.Button(self.page, width = 5, text = "default", command = self.defaultGuestFrozen)
        self.guestFrozenDefault.grid(row=actualRow, column=2)
        
        actualRow+=1
        
        self.guestAddH = Tkinter.Button(self.page, width = 5, text = "addH", command = self.addHGuest)
        self.guestAddH.grid(row=actualRow, column=2)
        
        actualRow+=1
        
        self.guestAddNCA = Tkinter.Button(self.page, width = 5, text = "chop", command = self.guestChop)
        self.guestAddNCA.grid(row=actualRow, column=2)
        
        actualRow+=1
        
        self.guestAtomsNo = Tkinter.Entry(self.page, width = 5)
        self.guestAtomsNo.grid(row=actualRow, column=2)
        
        actualRow+=1
        
#        self.guestAround = Tkinter.Button(self.page, text = "get", width = 5, command = self.readGuestFromSeleAround)
#        self.guestAround.grid(row=4, column=2)
        
    def addSele2Guest(self):
        try:
            stateNo = cmd.get_state()
            cmd.create("guest", "%guest or sele", stateNo)
        except:
            stateNo = cmd.get_state()
            cmd.create("guest", "sele", stateNo)
            
    def chop(self, sele):
        try:
            model = cmd.get_model( sele )
        except:
            return
        
        chain2resId2atomNames = {}
        
        for atom in model.atom:
            resnum = int(atom.resi)
            pdbname = atom.name
            chain = atom.chain
            
            if not chain in chain2resId2atomNames:
                chain2resId2atomNames[chain] = { resnum : set([ pdbname ]) }
            else:
                if resnum in chain2resId2atomNames[chain]:
                    chain2resId2atomNames[chain][resnum].add(pdbname)
                else:
                    chain2resId2atomNames[chain][resnum] = set([ pdbname ])
                
        for chain in chain2resId2atomNames:
            resId2atomNames = chain2resId2atomNames[chain]
            
            resIds2addN = set([])
            resIds2addC = set([])
            
            for resnum in resId2atomNames:
                if resId2atomNames[resnum] != set( [ "CA", "C", "O" ]) and resId2atomNames[resnum] != set( [ "CA", "N" ]):
                    resIds2addC.add(resnum+1)
                    resIds2addN.add(resnum-1)
                    
            resIds2addN-= set(resId2atomNames.keys())
            resIds2addC-= set(resId2atomNames.keys())
            
            stateNo = cmd.get_state()
            for resnum in resIds2addC:
                cmd.create( sele, " %"+sele+" or ( ( resi "+str(resnum)+ " and chain "+chain+ " ) and ( name CA or name N) ) ", stateNo)
                cmd.bond( "%"+sele+" and name N and resi "+str(resnum), "%guest and name C and resi "+str(resnum-1), 1 )
                
            for resnum in resIds2addN:
                cmd.create( sele, " %"+sele+" or ( ( resi "+str(resnum)+ " and chain "+chain+ " ) and ( name CA or name C or name O) ) ", stateNo)
                cmd.bond( "%"+sele+" and name C and resi "+str(resnum), "%guest and name N and resi "+str(resnum+1), 1 )
            
        cmd.show("sticks", sele )
        
    def guestChop(self):
        self.chop("guest")
        
    def clearGuestFrozen(self):
        cmd.select( "guestFrozen", "none")
        
    def defaultGuestFrozen(self):
        cmd.select("guestFrozen", " %guestFrozen or ( guest and ( name CA  ) ) ")
        
    def addHGuest(self):
        cmd.h_add("guest")
        
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
        
        self.complexAtomsNo = Tkinter.Entry(self.page, width = 5)
        self.complexAtomsNo.grid(row=9, column=3)
        
    
    def gridWrite(self):
        writeButton = Tkinter.Button(self.page, text = "write", command = self.write)
        writeButton.grid(row=10, column = 0)
        
    def write(self):
        self.chargeTotal = self.complexCharge.get()
        self.spinTotal = self.complexSpin.get()
        self.chargeGuest = self.guestCharge.get()
        self.spinGuest = self.guestSpin.get()
        self.chargeHost = self.hostCharge.get()
        self.spinHost = self.hostSpin.get()
        
        if not self.chargeTotal or not self.chargeGuest or not self.chargeHost:
            tkMessageBox.showwarning(title = "Error!", message = "Please fill the charge data")
            return
        
        if not self.spinTotal or not self.spinGuest or not self.spinHost:
            tkMessageBox.showwarning(title = "Error!", message = "Please fill the spin data")
            return
        
        options = { "mustexist" : False , "title" : "Job directory selection"}
        directory = tkFileDialog.askdirectory(**options)
        if not directory:
            return
        
        if not os.path.isdir(directory):
            os.makedirs(directory)
            
        basename = os.path.basename(directory)
        
        cmd.save( os.path.join(directory, basename+".pdb" ), "host or guest"  )
        
        additional = open( os.path.join(directory, basename+".dat"), 'w' )
        objects = cmd.get_object_list()
        for obj in objects:
            if obj != "host" and obj != "guest":
                additional.write(obj+"\n")
        
        
        additional.close()
        
        if self.keywordSet:
            for setName in self.keywordSet:
                self.combineSet(directory, setName)
                
        else:
            self.writeInputJob(directory, basename, {})
            
    def combineSet(self, directory, setName):
        queue = [{}]
        
        for key in self.keywordSet[setName]:
            newQueue = []
            
            for value in self.keywordSet[setName][key]:
            
                for element in queue:
                    newElement = deepcopy(element)
                    newElement[key] = value
                    newQueue.append(newElement)
                    
            queue = newQueue
            
        for jobDict in queue:
            basename = setName
            for key in jobDict:
                newPart = jobDict[key]
                newPart = newPart.replace(" ","")
                newPart = newPart.replace("=","")
                newPart = newPart.replace("(","")
                newPart = newPart.replace(")", "")
                
                basename += "_" + newPart
            
            newDirectory = os.path.join(directory, basename)
            os.makedirs(newDirectory)
            
            self.writeInputJob( newDirectory, basename, jobDict)
        
    def writeInputJob(self, directory, basename, keywordDict):
        
        slurmFile = os.path.join( directory, basename + ".slurm" )
        inpFile = os.path.join(directory, basename + ".inp")
        inpFileBasename = os.path.basename(inpFile)
        
        slurmF = open(slurmFile, 'w')
        
        slurmHead = self.slurmTextG16.get("1.0", "end")
        slurmHead = self.updateText(slurmHead, keywordDict)
        
        slurmF.write(slurmHead)
        slurmF.write("\n")
        
        slurmF.write("\nmodule add plgrid/apps/gaussian/g16.A.03\n\n")
        slurmF.write("g16 "+inpFileBasename+ "\n\n")
        
        slurmF.close()
        
        inpF = open(inpFile, 'w')
        
        inpF.write("%Chk="+inpFileBasename.replace(".inp",".chk")+"\n")
        routeSection = self.routeSection.get("1.0", "end")
        routeSection = self.updateText(routeSection, keywordDict)
        inpF.write(routeSection)
        
        inpF.write("\nEmilka jest najpiekniejsza!\n\n")
        
        inpF.write(self.chargeTotal+","+self.spinTotal+" "+self.chargeGuest+","+self.spinGuest+" "+self.chargeHost+","+self.spinHost+"\n")
        
        model = cmd.get_model("guest and guestFrozen")
        self.writeModelToFile(model, inpF, 2, True)
            
        model = cmd.get_model("guest and not guestFrozen")
        self.writeModelToFile(model,inpF, 1)
            
        model = cmd.get_model("host and hostFrozen")
        self.writeModelToFile(model, inpF, 2, True)
            
        model = cmd.get_model("host and not hostFrozen")
        self.writeModelToFile(model, inpF, 2)
            
        inpF.write("\n")
        
        additionaInput = self.additionalSection.get("1.0", "end")
        additionaInput = self.updateText(additionaInput, keywordDict)
        inpF.write(additionaInput)
        inpF.write("\n\n")
        
        inpF.close()
        
    def writeModelToFile(self, model, file2append, fragmentNo, frozen = False):
        for atom in model.atom:
            resnum = atom.resi
            resname = atom.resn
            pdbname = atom.name
            element = atom.symbol
            
            x = str(atom.coord[0])
            y = str(atom.coord[1])
            z = str(atom.coord[2])
            
            if not frozen:
                file2append.write(element+"(Fragment="+str(fragmentNo)+", PDBName="+pdbname+", ResName="+resname+", ResNum="+str(resnum)+") "+x+" "+y+" "+z+"\n" )
            else:
                file2append.write(element+"(Fragment="+str(fragmentNo)+", PDBName="+pdbname+", ResName="+resname+", ResNum="+str(resnum)+") -1"+x+" "+y+" "+z+"\n" )
                
    def getKeywordsFromText(self, text):
        keywords = []
        
        insideKeyword = False
        
        for letter in text:
            if letter == "{":
                insideKeyword = True
                keywords.append("")
                continue
            elif letter == "}":
                insideKeyword = False
                
            if insideKeyword:
                keywords[-1] += letter
                
        return keywords
    
    def updateText(self, text, keyDict):
        keywords = self.getKeywordsFromText(text)
        
        for key in keywords:
            if key in keyDict:
                text = text.replace("{"+key+"}", keyDict[key])
            else:
                text = text.replace("{"+key+"}", " ")
                
        return text
        
    def gridGaussianRouteSection(self):
        routeSectionLabel = Tkinter.Label(self.page, text = "Gaussian route section")
        routeSectionLabel.grid(row = 1, column=5, columnspan = 5)
        
        self.routeSection =Tkinter.Text(self.page, width =50, height = 10 )
        self.routeSection.grid(row = 2, column = 5, columnspan = 5, rowspan = 5)
        self.routeSection.insert("end", "%Mem=100GB\n#P B3LYP/6-31G(d,p)\n# Opt Counterpoise=2\n# SCRF(Solvent=Water, Read)\n# Gfinput IOP(6/7=3)  Pop=full  Density  Test \n# Units(Ang,Deg)")
        
    def gridGaussianAdditionalInput(self):
        additionalInputLabel = Tkinter.Label(self.page, text="Additional input")
        additionalInputLabel.grid(row=7, column = 5, columnspan = 5)
        
        self.additionalSection = Tkinter.Text(self.page, width = 50, height = 10)
        self.additionalSection.grid(row=8, column=5, columnspan = 5, rowspan = 5)
        self.additionalSection.insert("end", "eps=4\n")
                                 
    def gridSlurmSection(self):
        slurmSection = Tkinter.Label(self.page, text = "Slurm input")
        slurmSection.grid(row = 1, column = 10, columnspan = 10)
        
        self.slurmTextG16 = Tkinter.Text(self.page, width = 50, height = 10)
        self.slurmTextG16.grid(row = 2, column = 10 , columnspan = 5, rowspan = 5)
        self.slurmTextG16.insert("end" , "#!/bin/env bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=24\n#SBATCH --time=70:00:00\n##### Nazwa kolejki\n#SBATCH -p plgrid\n" )
                  
    def gridSetKeywordsSection(self):
        setLabel = Tkinter.Label(self.page, text = "Set")
        setLabel.grid(row = 7, column = 10, columnspan = 2)
        
        self.setListbox = Tkinter.Listbox(self.page, width =12, height = 9, exportselection = False)
        self.setListbox.grid(row= 8, column = 10, columnspan = 2, rowspan = 5)
        self.setListbox.bind("<<ListboxSelect>>", self.selectSet )
        
        setNewButton = Tkinter.Button(self.page, width = 4, text = "New:", command = self.addSet)
        setNewButton.grid(row = 20, column = 10, columnspan = 1)
        
        setCopyButton = Tkinter.Button(self.page, width = 4, text = "Copy:", command = self.copySet)
        setCopyButton.grid(row = 20, column = 11, columnspan = 1)
        
        self.setNewEntry = Tkinter.Entry(self.page, width = 8)
        self.setNewEntry.grid(row = 21, column = 10, columnspan = 2)
        
        setDeleteButton = Tkinter.Button(self.page, width = 8, text = "Delete:", command = self.deleteSet)
        setDeleteButton.grid(row = 22, column = 10, columnspan = 2)
        
        #########################################################
        
        keysLabel = Tkinter.Label(self.page, text = "Keys")
        keysLabel.grid(row = 7, column = 12, columnspan = 2)
        
        self.keysListbox = Tkinter.Listbox(self.page, width =12, height = 9, exportselection = False)
        self.keysListbox.grid(row= 8, column = 12, columnspan = 2, rowspan = 5)
        self.keysListbox.bind("<<ListboxSelect>>", self.selectKey)
        
        keysNewButton = Tkinter.Button(self.page, width = 4, text = "New:", command = self.addKey)
        keysNewButton.grid(row = 20, column = 12, columnspan = 1)
        
        keysCopyButton = Tkinter.Button(self.page, width = 4, text = "Copy:", command = self.copyKey)
        keysCopyButton.grid(row = 20, column = 13, columnspan = 1)
        
        self.keysNewEntry = Tkinter.Entry(self.page, width = 8)
        self.keysNewEntry.grid(row = 21, column = 12, columnspan = 2)
        
        keysDeleteButton = Tkinter.Button(self.page, width = 8, text = "Delete:", command = self.deleteKey)
        keysDeleteButton.grid(row = 22, column = 12, columnspan = 2)
        
        #########################################################
        
        valuesLabel = Tkinter.Label(self.page, text = "Values")
        valuesLabel.grid(row = 7, column = 14, columnspan = 2)
        
        self.valuesListbox = Tkinter.Listbox(self.page, width =12, height = 9, exportselection = False)
        self.valuesListbox.grid(row= 8, column = 14, columnspan = 2, rowspan = 5)
        
        valuesNewButton = Tkinter.Button(self.page, width = 8, text = "New:", command = self.addValues)
        valuesNewButton.grid(row = 20, column = 14, columnspan = 2)
        
        self.valuesNewEntry = Tkinter.Entry(self.page, width = 8)
        self.valuesNewEntry.grid(row = 21, column = 14, columnspan = 2)
        
        valuesDeleteButton = Tkinter.Button(self.page, width = 8, text = "Delete:", command = self.deleteValue)
        valuesDeleteButton.grid(row = 22, column = 14, columnspan = 2)
        
    def addSet(self):
        newRecord = self.setNewEntry.get()
        
        if not newRecord:
            return
        
        self.keywordSet[newRecord] = {}
        self.setListbox.insert(0, newRecord)
        self.setListbox.selection_clear(0, "end")
        self.setListbox.selection_set(0)
        self.selectSet(0)
        
        self.setNewEntry.delete(0, "end")
        
    def copySet(self):
        newRecord = self.setNewEntry.get()
        
        if not newRecord:
            return
        
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
            tkMessageBox.showwarning(title = "Error!", message = "Please select the set")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        self.keywordSet[newRecord] = deepcopy( self.keywordSet[selectedSet] )
        self.setListbox.insert(0, newRecord)
        self.setListbox.selection_clear(0, "end")
        self.setListbox.selection_set(0)
        self.selectSet(0)
        
        self.setNewEntry.delete(0, "end")
    
    def addKey(self):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
            tkMessageBox.showwarning(title = "Error!", message = "Please select the set")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        newRecord = self.keysNewEntry.get()
        
        if not newRecord:
            return
        
        self.keywordSet[selectedSet][newRecord] = []
        self.keysListbox.insert(0, newRecord)
        self.keysListbox.selection_clear(0, "end")
        self.keysListbox.selection_set(0)
        self.selectKey(0)
        
        self.keysNewEntry.delete(0, "end")
        
    def copyKey(self):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
            tkMessageBox.showwarning(title = "Error!", message = "Please select the set")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        selectedKey = self.keysListbox.curselection()
        if not selectedKey:
            tkMessageBox.showwarning(title = "Error!", message = "Please select the key")
            return
        
        selectedKey = self.keysListbox.get(selectedKey)
        
        newRecord = self.keysNewEntry.get()
        
        if not newRecord:
            return
        
        self.keywordSet[selectedSet][newRecord] = deepcopy( self.keywordSet[selectedSet][selectedKey] )
        self.keysListbox.insert(0, newRecord)
        self.keysListbox.selection_clear(0, "end")
        self.keysListbox.selection_set(0)
        self.selectKey(0)
        
        self.keysNewEntry.delete(0, "end")
    
    def addValues(self):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
            tkMessageBox.showwarning(title = "Error!", message = "Please select the set")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        selectedKey = self.keysListbox.curselection()
        if not selectedKey:
            tkMessageBox.showwarning(title = "Error!", message = "Please select the key")
            return
        
        selectedKey = self.keysListbox.get(selectedKey)
        
        newRecord = self.valuesNewEntry.get()
        
        if not newRecord:
            return
        
        if newRecord in self.keywordSet[selectedSet][selectedKey]:
            return
        
        self.keywordSet[selectedSet][selectedKey] = [ newRecord ] + self.keywordSet[selectedSet][selectedKey]
        self.valuesListbox.insert(0, newRecord)
        self.valuesListbox.selection_clear(0, "end")
        self.valuesListbox.selection_set(0)
        
        self.valuesNewEntry.delete(0, "end")
        
    def selectSet(self, arg):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
            print("nihuhu")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        self.keysListbox.delete(0,"end")
        self.valuesListbox.delete(0, "end")
        
        for key in self.keywordSet[selectedSet]:
            self.keysListbox.insert("end", key)
    
    def selectKey(self, arg):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
#            tkMessageBox.showwarning(title = "Error!", message = "Please select the set")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        selectedKey = self.keysListbox.curselection()
        if not selectedKey:
            print("nihuhu")
            return
        
        selectedKey = self.keysListbox.get(selectedKey)
        self.valuesListbox.delete(0, "end")
        
        for value in self.keywordSet[selectedSet][selectedKey]:
            self.valuesListbox.insert("end", value)
        
    def deleteSet(self):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
            print("nihuhu")
            return
        
        selectedSetName = self.setListbox.get(selectedSet)
        self.setListbox.delete(selectedSet)
        
        self.keysListbox.delete(0,"end")
        self.valuesListbox.delete(0, "end")
        
        del self.keywordSet[selectedSetName]
    
    def deleteKey(self):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
#            tkMessageBox.showwarning(title = "Error!", message = "Please select the set")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        selectedKey = self.keysListbox.curselection()
        if not selectedKey:
            print("nihuhu")
            return
        
        selectedKeyName = self.keysListbox.get(selectedKey)
        self.keysListbox.delete(selectedKey)
        
        self.valuesListbox.delete(0, "end")
        del self.keywordSet[selectedSet][selectedKeyName]
        
    
    def deleteValue(self):
        selectedSet = self.setListbox.curselection()
        if not selectedSet:
#            tkMessageBox.showwarning(title = "Error!", message = "Please select the set")
            return
        
        selectedSet = self.setListbox.get(selectedSet)
        
        selectedKey = self.keysListbox.curselection()
        if not selectedKey:
            print("nihuhu")
            return
        
        selectedKey = self.keysListbox.get(selectedKey)
        
        selectedValue = self.valuesListbox.curselection()
        if not selectedValue:
            return
        
        selectedValueName = self.valuesListbox.get(selectedValue)
        self.valuesListbox.delete(selectedValue)
        
        self.keywordSet[selectedSet][selectedKey].remove(selectedValueName)
        
        
    def insertKeywordSet(self):
        self.setListbox.delete(0, "end")
        self.keysListbox.delete(0,"end")
        self.valuesListbox.delete(0, "end")
        
        for setName in self.keywordSet:
            self.setListbox.insert("end", setName)
    
    def grid(self):
        self.gridBasicLabels()
        self.gridHost()
        self.gridGuest()
        self.gridComplex()
        self.gridWrite()
        self.gridSlurmSection()
        self.gridSetKeywordsSection()
        
        self.gridGaussianRouteSection()
        self.gridGaussianAdditionalInput()
        
    def getState(self):
        state = {  "inputBegin" : "", "additionalInput" : "", "slurmConfig" : "" , "keywordSet" : {} }
        
        state["inputBegin"] = self.routeSection.get("1.0", "end")
        state["additionalInput"] = self.additionalSection.get("1.0", "end")
        state["slurmConfig"] = self.slurmTextG16.get("1.0", "end")
        state["keywordSet"] = self.keywordSet
        
        return state
    
    def loadState(self, state):
        if "inputBegin" in state:
            self.routeSection.delete("1.0", "end")
            self.routeSection.insert("end", state["inputBegin"] )
        
        if "additionalInput" in state:
            self.additionalSection.delete("1.0", "end")
            self.additionalSection.insert("end", state["additionalInput"])
        
        if "slurmConfig" in state:
            self.slurmTextG16.delete("1.0", "end")
            self.slurmTextG16.insert("end", state["slurmConfig"])
        
        if "keywordSet" in state:
            self.keywordSet = state["keywordSet"]
            self.insertKeywordSet()
        