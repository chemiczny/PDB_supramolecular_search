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
from simpleFilters import noAAinAnions, noNUinAnions, noNUinHDonors

class HBondsGUI(SupramolecularGUI):
    def __init__(self, page):
        SupramolecularGUI.__init__(self, page)
        self.setNumericalParameters({ "R" : {"header" : "Distance"}  })

        self.setListParameters({ "Acceptor" : { "header" : "Anion code" }, 
                  "Acceptor el." : { "header" : "Acceptor atom" } , 
                  "Donor" : { "header" : "Donor code" },
                  "Donor el." : { "header" : "Donor atom"} })

        self.setSortingParameters({  "R" : "Distance", "Donor" : "Donor code", "Don. el." : "Donor atom" , "Acceptor" : "Anion code",
                     "Acc. el." : "Acceptor atom" }, [  "R" , "Donor" , "Don. el." , "Acceptor" ,
                     "Acc. el."  ])

        self.setTreeData([ "ID" , "PDB" , "Acceptor", "Anion id", "Acc. el." , "Donor", "Donor id", "Donor el.", 
                          "R" ])

        self.setAdditionalCheckboxes( [ { "label" : "No AA in acceptors", "func" : noAAinAnions } ,
                                      # { "label" : "No AA in donors", "func" : noAAinHDonors } ,
                                       { "label" : "No NU in acceptors", "func" : noNUinAnions } ,
                                          { "label" : "No NU in donors", "func" : noNUinHDonors } ]  )
    
        self.arrowName = "HBondArrow"
        self.arrowColor = "red violet"
        
    def getValues(self, rowId, row):
        return ( rowId, row["PDB Code"] , row["Anion code"], 
                                                    row["Anion chain"]+str(row["Anion id"]) ,row["Acceptor atom"], 
                                                    row["Donor code"], 
                                                    row["Donor chain"] + str(row["Donor id"]), row["Donor atom"], 
                                                    str(row["Distance"])[:3] ) 
        
    def getSelection(self,  data):
        res1Id = data["Anion id"].values[0]
        res1Chain = data["Anion chain"].values[0]
        
        res2Id = data["Donor id"].values[0]
        res2Chain = data["Donor chain"].values[0]
        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getSelectionFromRow(self,  data):
        res1Id = data["Anion id"]
        res1Chain = data["Anion chain"]
        
        res2Id = data["Donor id"]
        res2Chain = data["Donor chain"]
#        cmd.hide("everything")
        selection =  "( " + "chain "+res1Chain +" and resi "+ str(res1Id) +" ) or ( " +" chain " + res2Chain + " and resi "+str(res2Id)+")"
        return selection
    
    def getArrow(self, data ):
        point1Coords = [ data["Acceptor x coord"].values[0] , data["Acceptor y coord"].values[0] , data["Acceptor z coord"].values[0] ]
        point2Coords = [ data["Donor x coord"].values[0] , data["Donor y coord"].values[0] , data["Donor z coord"].values[0] ]
        
        return point1Coords, point2Coords
    
    def getArrowFromRow(self, data ):
        point1Coords = [ data["Acceptor x coord"] , data["Acceptor y coord"] , data["Acceptor z coord"] ]
        point2Coords = [ data["Donor x coord"] , data["Donor y coord"] , data["Donor z coord"] ]
        
        return point1Coords, point2Coords