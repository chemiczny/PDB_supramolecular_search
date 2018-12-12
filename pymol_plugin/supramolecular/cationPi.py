#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:32:57 2018

@author: michal
"""
import sys
if sys.version_info[0] < 3:
    from pymol import cmd
from supramolecularGUI import SupramolecularGUI
from simpleFilters import noAAinPiAcids, noNUinPiAcids

class CationPiGUI(SupramolecularGUI):
    def __init__(self, page):
        SupramolecularGUI.__init__(self, page)
        self.setNumericalParameters( { "R" : {"header" : "Distance"}, "h" : {"header" : "h"}, "x" : { "header" : "x" },
                           "alpha" : { "header" : "Angle" } , "Chain size" : {"header" : "RingChain" } } )
    
        self.setListParameters({ "Pi acid" : { "header" : "Pi acid Code" }, 
                      "Cation" : { "header" : "Cation code" },
                      "Element" : { "header" : "Atom symbol" },
                      "Chain" : {"header" : "RingChain" },
                      "Complex" : {"header" : "Complex" },
                      "CoordNo" : {"header" : "CoordNo"} ,
                      "PDB" : { "header" : "PDB Code" }
                      })
    
    
        self.setSortingParameters({  "R" : "Distance", "Angle" : "Angle", "x" : "x" , "h" : "h",
                          "Pi acid" : "Pi acid Code" ,"Cation" : "Cation code", "Cat. el." : "Atom symbol" }, [  "R" , "Angle" , "x" , "h" ,
                          "Pi acid" ,"Cation" , "Cat. el." ] )
    
        self.setTreeData([ "ID" , "PDB" , "Pi acid", "Pi acid id", "Cation", "Cation id", "Cat. el." , "R", "alpha", "x", "h", "chain", "chain Flat",
                          "complex", "coord No" ])
    
        self.setAdditionalCheckboxes( [ { "label" : "No AA in Pi acids", "func" : noAAinPiAcids } ,
                                           { "label" : "No NU in Pi acids", "func" : noNUinPiAcids } ]  )
    
        self.arrowName = "cationPiArrow"
        self.arrowColor = "blue orange"
    
    def getValues(self, rowId, row):
        return ( rowId, row["PDB Code"] , row["Pi acid Code"], 
                                                        row["Pi acid chain"]+str(row["Piacid id"]) , row["Cation code"], 
                                                        row["Cation chain"] + str(row["Cation id"]), row["Atom symbol"], 
                                                        str(row["Distance"])[:3], str(row["Angle"])[:4], str(row["x"])[:3],
                                                        str(row["h"])[:3] , str(row["RingChain"]), str(row["ChainFlat"]) , str(row["Complex"]), str(row["CoordNo"]) )
        
    def getSelection(self,  data):
        res1Id = data["Piacid id"].values[0]
        res1Chain = data["Pi acid chain"].values[0]
        
        res2Id = data["Cation id"].values[0]
        res2Chain = data["Cation chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getSelectionFromRow(self,  data):
        res1Id = data["Piacid id"]
        res1Chain = data["Pi acid chain"]
        
        res2Id = data["Cation id"]
        res2Chain = data["Cation chain"]
#        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getArrow(self, data ):
        point1Coords = [ data["Centroid x coord"].values[0] , data["Centroid y coord"].values[0] , data["Centroid z coord"].values[0] ]
        point2Coords = [ data["Cation x coord"].values[0] , data["Cation y coord"].values[0] , data["Cation z coord"].values[0] ]
        
        return point1Coords, point2Coords
    
    def getArrowFromRow(self, data ):
        point1Coords = [ data["Centroid x coord"] , data["Centroid y coord"] , data["Centroid z coord"] ]
        point2Coords = [ data["Cation x coord"] , data["Cation y coord"] , data["Cation z coord"] ]
        
        return point1Coords, point2Coords