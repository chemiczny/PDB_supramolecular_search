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

class PlanarAnionPiGUI(SupramolecularGUI):
    def __init__(self, page):
        SupramolecularGUI.__init__(self, page)
        self.setNumericalParameters({ "alpha" : { "header" : "Angle" } }  )
    
        self.setListParameters( { "Pi acid" : { "header" : "Pi acid Code" }, 
                      "Anions" : { "header" : "Anion code" } , 
                      "PDB" : { "header" : "PDB Code" } }  )
    
        self.setSortingParameters({  "Angle" : "Angle", "Pi acid" : "Pi acid Code" ,"Anion" : "Anion code" }, [   "Angle" ,  "Pi acid" ,"Anion" ] )
    
        self.setTreeData([ "ID" , "PDB" , "Pi acid", "Pi acid id", "Anion", "Anion id" ,"Anion gr. id",  "alpha" ])
     
        self.setAdditionalCheckboxes( [ { "label" : "No AA in Pi acids", "func" : noAAinPiAcids } ,
                                         { "label" : "No AA in anions" ,  "func" : noAAinAnions } ,
                                         { "label" : "No NU in Pi acids", "func" : noNUinPiAcids } ,
                                         { "label" : "No NU in anions", "func" : noNUinAnions } ]  )
    
        self.arrowName = "planarAnionPiArrow"
        self.arrowColor = "blue red"
    
    def getValues(self, rowId, row):
        return ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                        row["Pi acid chain"]+str(row["Piacid id"]) , row["Anion code"], 
                                                        row["Anion chain"] + str(row["Anion id"]), row["Anion group id"],
                                                        str(row["Angle"])[:4]  )
        
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
        anionAtomCoords = [ data["Anion group x coord"].values[0] , data["Anion group y coord"].values[0] , data["Anion group z coord"].values[0] ]
        
        return centroidCoords, anionAtomCoords
    
    def getArrowFromRow(self, data ):
        centroidCoords = [ data["Centroid x coord"] , data["Centroid y coord"] , data["Centroid z coord"] ]
        anionAtomCoords = [ data["Anion group x coord"] , data["Anion group y coord"] , data["Anion group z coord"] ]
        
        return centroidCoords, anionAtomCoords