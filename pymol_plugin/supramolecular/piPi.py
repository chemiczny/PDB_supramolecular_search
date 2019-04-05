#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:33:21 2018

@author: michal
"""
import sys
if sys.version_info[0] < 3:
    from pymol import cmd
from supramolecularGUI import SupramolecularGUI
from simpleFilters import noAAinPiAcids, noAAinPiRes, noNUinPiAcids, noNUinPiRes

class PiPiGUI(SupramolecularGUI):
    def __init__(self, page):
        SupramolecularGUI.__init__(self, page)
        self.setNumericalParameters( { "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" },
                           "alpha" : { "header" : "Angle" }, "theta" : { "header" : "theta" } , "omega" : { "header" : "omega" }  }  )
    
        self.setListParameters( { "Pi 1" : { "header" : "Pi acid Code" }, 
                      "Pi 2" : { "header" : "Pi res code" }  ,
                      "PDB" : { "header" : "PDB Code" }} )
    
        self.setSortingParameters({  "R" : "Distance", "Angle" : "Angle", "x" : "x" , "h" : "h",
                         "theta" : "theta", "Pi 1" : "Pi acid Code" ,"Pi 2" : "Pi res code" , "omega" : "omega" }, [  "R" , "Angle" , "x" , "h" ,
                         "theta" , "Pi 1" ,"Pi 2" , "omega" ] )
    
        self.setTreeData([ "ID" , "PDB" , "Pi 1", "Pi 1 id", "Pi 2", "Pi 2 id" , "R", "alpha", "x", "h", "theta", "omega" ])
    
        self.setAdditionalCheckboxes( [ { "label" : "No AA in Pi 1", "func" : noAAinPiAcids } ,
                                         { "label" : "No AA in Pi 2" ,  "func" : noAAinPiRes },
                                         { "label" : "No NU in Pi 1", "func" : noNUinPiAcids },
                                         { "label" : "No NU in Pi 2", "func" : noNUinPiRes } ]  )
    
        self.arrowName = "piPiArrow"
        self.arrowColor = "blue green"
    
    def getValues(self, rowId, row):
        return ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                        row["Pi acid chain"]+str(row["Piacid id"]) , row["Pi res code"], 
                                                        row["Pi res chain"] + str(row["Pi res id"]), 
                                                        str(row["Distance"])[:3], str(row["Angle"])[:4], str(row["x"])[:3],
                                                        str(row["h"])[:3], str(row["theta"])[:4] , str(row["omega"])[:4]  )
        
    def getSelection(self,  data):
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Pi res id"].values[0]
        res2Chain = data["Pi res chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getSelectionFromRow(self,  data):
        res1Id = data["Piacid id"]
        res1Chain = data["Pi acid chain"]
        
        res2Id = data["Pi res id"]
        res2Chain = data["Pi res chain"]
#        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getArrow(self, data ):
        point1Coords = [ data["Centroid x coord"].values[0] , data["Centroid y coord"].values[0] , data["Centroid z coord"].values[0] ]
        point2Coords = [ data["Centroid 2 x coord"].values[0] , data["Centroid 2 y coord"].values[0] , data["Centroid 2 z coord"].values[0] ]
        
        return point1Coords, point2Coords
    
    def getArrowFromRow(self, data ):
        point1Coords = [ data["Centroid x coord"] , data["Centroid y coord"] , data["Centroid z coord"] ]
        point2Coords = [ data["Centroid 2 x coord"] , data["Centroid 2 y coord"] , data["Centroid 2 z coord"] ]
        
        return point1Coords, point2Coords