#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:32:43 2018

@author: michal
"""
try:
    from pymol import cmd
except:
    pass
from supramolecularGUI import SupramolecularGUI
from simpleFilters import noAAinPiAcids, noAAinAnions, noNUinPiAcids, noNUinAnions

class AnionPiGUI(SupramolecularGUI):
    def __init__(self, page, parallelSelectionFunction, name):
        SupramolecularGUI.__init__(self, page, parallelSelectionFunction, name)
        self.setNumericalParameters({ "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" },
                           "alpha" : { "header" : "Angle" }, "res" : { "header" : "Resolution" }  }  )
    
        self.setListParameters( { "Pi acid" : { "header" : "Pi acid Code" }, 
                      "Anions" : { "header" : "Anion code" } , 
                      "Groups" : { "header" : "Anion type" },
                      "Methods" : { "header" : "Method" } ,
                      "PDB" : { "header" : "PDB Code" },
                      "Ring" : { "header" : "Ring elements" },
                      "Type" : { "header" : "Structure type" } }  )
    
        self.setSortingParameters({  "R" : "Distance", "Angle" : "Angle", "x" : "x" , "h" : "h",
                         "res" : "Resolution", "Pi acid" : "Pi acid Code" ,"Anion" : "Anion code" }, [  "R" , "Angle" , "x" , "h" ,
                         "res" , "Pi acid" ,"Anion" ] )
    
        self.setTreeData([ "ID" , "PDB" , "Pi acid", "Pi acid id", "Anion", "Anion id", "Anion type" ,"Anion gr. id", "R", "alpha", "x", "h", "res", "Method", "Type" ])
     
        self.setAdditionalCheckboxes( [ { "label" : "No AA in Pi acids", "func" : noAAinPiAcids } ,
                                         { "label" : "No AA in anions" ,  "func" : noAAinAnions } ,
                                         { "label" : "No NU in Pi acids", "func" : noNUinPiAcids } ,
                                         { "label" : "No NU in anions", "func" : noNUinAnions } ]  )
    
        self.setUniqueParameters( { "PDB" : ["PDB Code"], "Pi acid" : [ "Pi acid Code" , "Pi acid chain", "Piacid id"  ],
                                   "Anion" : [ "Anion code" , "Anion chain", "Anion id" ], "Anion id" : ["Anion group id"] ,
                                   "Ring id" : ["CentroidId"] ,  "Model" : [ "Model No" ] } , [ "PDB" , "Pi acid" , "Ring id", "Anion" , "Anion id" , "Model" ]  )
    
        self.interactionButtons = [ { "name" : "All anion int." , "headers" : ["PDB Code", "Anion code", "Anion chain", "Anion id", "Model No"] },
                                     { "name" : "All pi acid int." , "headers" : [ "PDB Code", "Pi acid Code" , "Pi acid chain", "Piacid id", "Model No" ]  }]
    
        self.arrowName = "anionPiArrow"
        self.arrowColor = "blue red"
    
    def getValues(self, rowId, row):
        return ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                        row["Pi acid chain"]+str(row["Piacid id"]) , row["Anion code"], 
                                                        row["Anion chain"] + str(row["Anion id"]), row["Anion type"], row["Anion group id"],
                                                        str(row["Distance"])[:3], str(row["Angle"])[:4], str(row["x"])[:3],
                                                        str(row["h"])[:3], row["Resolution"] , row["Method"], row["Structure type"]  )
        
    def getSelection(self,  data):
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Anion id"].values[0]
        res2Chain = data["Anion chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getSelectionFromRow(self,  data):
        res1Id = data["Piacid id"]
        res1Chain = data["Pi acid chain"]
        
        res2Id = data["Anion id"]
        res2Chain = data["Anion chain"]
#        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getArrow(self, data ):
        centroidCoords = [ data["Centroid x coord"].values[0] , data["Centroid y coord"].values[0] , data["Centroid z coord"].values[0] ]
        anionAtomCoords = [ data["Anion x coord"].values[0] , data["Anion y coord"].values[0] , data["Anion z coord"].values[0] ]
        
        return centroidCoords, anionAtomCoords
    
    def getArrowFromRow(self, data ):
        centroidCoords = [ data["Centroid x coord"] , data["Centroid y coord"] , data["Centroid z coord"] ]
        anionAtomCoords = [ data["Anion x coord"] , data["Anion y coord"] , data["Anion z coord"] ]
        
        return centroidCoords, anionAtomCoords