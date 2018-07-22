#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 19:40:59 2018

@author: michal
"""
import pandas as pd
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

from supramolecularGUI import SupramolecularGUI

def fetchdialog(simulation = False):
    if simulation:
        root = Tkinter.Tk()
    else:
        app = plugins.get_pmgapp()
        root = plugins.get_tk_root()
        
    appData = { "cifDir" : False }
    
    self = Tkinter.Toplevel(root)
    self.title('Supramolecular analyser')
    self.minsize(1100, 700)
    self.resizable(0,0)
    
    nb = ttk.Notebook(self, height = 700, width = 1100)

    pageAnionPi = ttk.Frame(nb)
    pagePiPi = ttk.Frame(nb)
    pageCationPi = ttk.Frame(nb)
    pageAnionCation  = ttk.Frame(nb)
    
    nb.add(pageAnionPi, text = "AnionPi")
    nb.add(pagePiPi, text = "PiPi")
    nb.add(pageCationPi, text = "CationPi")
    nb.add(pageAnionCation, text = "AnionCation")
    
    nb.grid(column = 0, row = 0, columnspan = 20)
    
    guiAnionPi = SupramolecularGUI(pageAnionPi)
    guiPiPi = SupramolecularGUI(pagePiPi)
    guiCationPi = SupramolecularGUI( pageCationPi)
    guiAnionCation = SupramolecularGUI( pageAnionCation)
    
    guis = [ guiAnionPi, guiPiPi, guiCationPi, guiAnionCation]
    
    ######################
    # GENERAL
    ######################
    
    def selectCif():
        appData["cifDir"] = tkFileDialog.askdirectory()
        if not appData["cifDir"]:
            print("nihuhu")
            return
        
        ent_cifDir.configure(state = "normal")
        ent_cifDir.delete(0,"end")
        ent_cifDir.insert(0, appData["cifDir"])
        ent_cifDir.configure(state = "readonly")
        
        for gui in guis:
            gui.logData["cifDir"] = appData["cifDir"]
        
        
    
    but_cifDir = Tkinter.Button(self, width = 10, command = selectCif, text = "Cif dir")
    but_cifDir.grid(row = 2, column = 1)
    
    ent_cifDir = Tkinter.Entry(self, width =45)
    ent_cifDir.configure(state = "readonly")
    ent_cifDir.grid(row = 2, column = 2, columnspan = 3)
    
    actionLabels = [ "AnionPi", "PiPi", "CationPi", "AnionCation" ]
    actionLabels2Objects = { "AnionPi" : guiAnionPi, "PiPi" : guiPiPi, "CationPi" : guiCationPi, "AnionCation" : guiAnionCation }
    actionMenu = {}
    
    column = 1
    
    lab_usePage = Tkinter.Label(self, width = 10, text = "Use:")
    lab_usePage.grid(row = 4, column = 0)
    for label in actionLabels:
        actionMenu[label] = {}
        
        actionMenu[label]["label"] = Tkinter.Label(self, text = label)
        actionMenu[label]["label"].grid(row = 3, column = column )
        
        actual_row = 4
        
        actionMenu[label]["checkValue"] = Tkinter.IntVar()
        actionMenu[label]["checkbox"] = Tkinter.Checkbutton(self, variable = actionMenu[label]["checkValue"])
        actionMenu[label]["checkbox"].grid(row = actual_row, column = column)
            
        column += 1
    
    
    def mergeResults():
        headersId ={ 
                    "Pi" : [ "Pi acid Code" , "Pi acid chain" , "Piacid id"] , 
                    "Anion" : ["Anion code", "Anion chain" , "Anion id"] , 
                    "Cation" : ["Cation code", "Cation chain", "Cation id"]  }
        
        selectedData = []
        
        
        for label in actionLabels:
            if actionMenu[label]["checkValue"].get() > 0:
                if not "filtered" in  actionLabels2Objects[label].logData:
                    continue
                else:
                    selectedData.append(label)
                    
        if len(selectedData) < 2:
            return
        
        uniqueData = []
        actualKeys = []
        for key in selectedData:
            headers = [ "PDB Code", "Model No" ] 
            
            for headerKey in headersId:
                if headerKey in key:
                    headers += headersId[headerKey]
                    
            if len(uniqueData) == 0:
                uniqueData = actionLabels2Objects[ key ].logData["filtered"][ headers ].drop_duplicates()
            else:
                uniqueData = pd.merge( uniqueData,  actionLabels2Objects[ key ].logData["filtered"][ headers ], on = list( set(actualKeys) & set(headers) ))
                uniqueData = uniqueData.drop_duplicates()
                
            actualKeys = list(set( actualKeys + headers ))
            
            
        for key in selectedData:
            headers = [ "PDB Code", "Model No" ] 
            
            for headerKey in headersId:
                if headerKey in key:
                    headers += headersId[headerKey]
                    
            mergingKeys = list(set(actualKeys) & set(headers) )
            tempDataFrame = uniqueData[ mergingKeys   ].drop_duplicates()
            mergedData = pd.merge( actionLabels2Objects[ key ].logData["filtered"], tempDataFrame, on = mergingKeys )
            actionLabels2Objects[ key ].printFilterResults( mergedData )
            
            
        
    def countSubstringInList( list2check, substring ):
        found = 0
        
        for el in list2check:
            if substring in list2check:
                found +=1
                
        return found
    
    but_merge = Tkinter.Button(self, width = 20, text = "Merge!", command = mergeResults)
    but_merge.grid(row = 4, column = 6, columnspan = 2)
    
    def showAllInteractions():
        pass
    
    but_showAll = Tkinter.Button(self, width = 20, text = "All inter.", command = showAllInteractions)
    but_showAll.grid(row = 4, column = 8, columnspan = 2)

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
                      "Element" : { "header" : "Atom symbol" },
                      "Chain" : {"header" : "RingChain" }
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
    
    guiAnionCation.setListParameters({ "Cation" : { "header" : "Cation code" }, 
                      "Cat. el." : { "header" : "Cation symbol" } , 
                      "Anion" : { "header" : "Anion code" },
                      "An. at." : { "header" : "Anion symbol"}})
    
    guiAnionCation.setSortingParameters({  "R" : "Distance", "Cation" : "Cation code", "Cat. el." : "Cation symbol" , "Anion" : "Anion code",
                         "An. el." : "Anion symbol" }, [  "R" , "Cation" , "Cat. el." , "Anion" ,
                         "An. el."  ])
    
    guiAnionCation.setTreeData([ "ID" , "PDB" , "Cation", "Cation id", "Anion", "Anion id", "Anion el.", "Cat. el." , "R"])
    
    guiAnionCation.setAdditionalCheckboxes( [ { "label" : "No AA in Pi anions", "func" : noAAinAnions }   ]  )
    
    def row2ValuesAnionCation(rowId, row):
        return ( rowId, row["PDB Code"] , row["Cation code"], 
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
    
    for gui in guis:
        gui.grid()
    
    if simulation:
        self.mainloop()