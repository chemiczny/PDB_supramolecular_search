#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:33:21 2018

@author: michal
"""
try:
    from pymol import cmd
except:
    pass
from supramolecularGUI import SupramolecularGUI
#from simpleFilters import noAAinAnions, noNUinAnions

def onlyAnions(actualData):
    return actualData[ actualData["isAnion"] == True ]

def onlyComplexes(actualData):
    return actualData[ actualData["Complex"] == True ]

class MetalLigandGUI(SupramolecularGUI):
    def __init__(self, page):
        SupramolecularGUI.__init__(self, page)
        self.setNumericalParameters({ "R" : {"header" : "Distance"}  })

        self.setListParameters({ "Cation" : { "header" : "Cation code" }, 
                  "Cat. el." : { "header" : "Cation element" } , 
                  "Ligand" : { "header" : "Anion code" },
                  "Lig. el.." : { "header" : "Ligand element"},
                  "Anion type" : { "header" : "anionType"} ,
                  "coordNo" : {"header" : "CoordNo"},
                  "summary" : { "header" : "Summary" },
                  "PDB" : { "header" : "PDB Code" } })

        self.setSortingParameters({  "R" : "Distance", "Cation" : "Cation code", "Cat. el." : "Cation element" , "Ligand" : "Anion code",
                     "Lig. el." : "Ligand element", "Anion type" : "anionType" }, [  "R" , "Cation" , "Cat. el." , "Ligand" ,
                     "Lig. el." , "Anion type"  ])

        self.setTreeData([ "ID" , "PDB" , "Cation", "Cation id", "Cat. el." , "Ligand", "Ligand id", "Ligand el.", "Is Anion", "Anion type"  , "Anion gr. id",
                          "R", "Complex", "Summary", "Coord No" ])

        self.setAdditionalCheckboxes( [ { "label" : "Only anions as ligands", "func" : onlyAnions } ,
                                          { "label" : "Only complexes", "func" : onlyComplexes } ]  )
    
        self.setUniqueParameters( { "PDB" : ["PDB Code"], "Metal" : [ "Cation Code" , "Cation chain", "Cation id"  ],
                                   "Ligand" : [ "Anion code" , "Anion chain", "Anion id" ] , "Model" :  [ "Model No" ] } , [ "PDB" , "Metal" , "Ligand", "Model" ] )
    
        self.interactionButtons = [ { "name" : "All metal int." , "headers" : ["PDB Code", "Cation code", "Cation chain", "Cation id", "Model No"] }]
 
    
        self.arrowName = "MetalLigandArrow"
        self.arrowColor = "red orange"
        self.numOfSortingMenu = 1
        
    def getValues(self, rowId, row):
        return ( rowId, row["PDB Code"] , row["Cation code"], 
                                                    row["Cation chain"]+str(row["Cation id"]) ,row["Cation element"], 
                                                    row["Anion code"], 
                                                    row["Anion chain"] + str(row["Anion id"]), row["Ligand element"], 
                                                    str(row["isAnion"]), 
                                                    row["anionType"] , row["Anion group id"],
                                                    str(row["Distance"])[:3] , str(row["Complex"]) , str(row["Summary"]) , str(row["CoordNo"])) 
        
    def getSelection(self,  data):
        res1Id = data["Anion id"].values[0]
        res1Chain = data["Anion chain"].values[0]
        
        res2Id = data["Cation id"].values[0]
        res2Chain = data["Cation chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getSelectionFromRow(self,  data):
        res1Id = data["Anion id"]
        res1Chain = data["Anion chain"]
        
        res2Id = data["Cation id"]
        res2Chain = data["Cation chain"]
#        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getArrow(self, data ):
        point1Coords = [ data["Anion x coord"].values[0] , data["Anion y coord"].values[0] , data["Anion z coord"].values[0] ]
        point2Coords = [ data["Cation x coord"].values[0] , data["Cation y coord"].values[0] , data["Cation z coord"].values[0] ]
        
        return point1Coords, point2Coords
    
    def getArrowFromRow(self, data ):
        point1Coords = [ data["Anion x coord"] , data["Anion y coord"] , data["Anion z coord"] ]
        point2Coords = [ data["Cation x coord"] , data["Cation y coord"] , data["Cation z coord"] ]
        
        return point1Coords, point2Coords