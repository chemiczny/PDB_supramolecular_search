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
from simpleFilters import noAAinAnions, noNUinAnions, noAAinCations

class AnionCationGUI(SupramolecularGUI):
    def __init__(self, page):
        SupramolecularGUI.__init__(self, page)
        self.setNumericalParameters({ "R" : {"header" : "Distance"} ,
                                     "Lat dif" : { "header" : "Latitude diff" } })

        self.setListParameters({ "Cation" : { "header" : "Cation code" }, 
                  "Cat. el." : { "header" : "Cation symbol" } , 
                  "Anion" : { "header" : "Anion code" },
                  "An. at." : { "header" : "Anion symbol"},
                  "Pi acid" : { "header" : "Pi acid Code" },
                  "Same hem." : { "header" : "Same semisphere"} ,
                  "PDB" : { "header" : "PDB Code" } })

        self.setSortingParameters({  "R" : "Distance", "Cation" : "Cation code", "Cat. el." : "Cation symbol" , "Anion" : "Anion code",
                     "An. el." : "Anion symbol", "Lat dif" : "Latitude diff" }, [  "R" , "Cation" , "Cat. el." , "Anion" ,
                     "An. el." , "Lat dif"  ])

        self.setTreeData([ "ID" , "PDB" , "Cation", "Cation id", "Cat. el." , "Anion", "Anion id", "Anion el.", "Anion gr. id", "Pi acid", "Pi acid id"  , 
                          "R", "Same semisphere", "Lat dif" ])

        self.setAdditionalCheckboxes( [ { "label" : "No AA in anions", "func" : noAAinAnions } ,
                                          { "label" : "No NU in anions", "func" : noNUinAnions },
                                          { "label" : "No AA in cations", "func" : noAAinCations }]  )
    
        self.setUniqueParameters( { "PDB" : ["PDB Code"], "Cation" : [ "Cation code" , "Cation chain", "Cation id"  ],
                                   "Anion" : [ "Anion code" , "Anion chain", "Anion id" ] , "Model" : [ "Model No" ] } ,
                                [ "PDB" , "Cation" , "Anion", "Model" ] )
   
        self.interactionButtons = [ { "name" : "All cat. int." , "headers" : ["PDB Code", "Cation code", "Cation chain", "Cation id", "Model No"] },
                                     { "name" : "All anion int." , "headers" : [ "PDB Code", "Anion code" , "Anion chain", "Anion id", "Model No" ]  }]
    
    
        self.arrowName = "AnionCationArrow"
        self.arrowColor = "red orange"
        
    def getValues(self, rowId, row):
        return ( rowId, row["PDB Code"] , row["Cation code"], 
                                                    row["Cation chain"]+str(row["Cation id"]) ,row["Cation symbol"], 
                                                    row["Anion code"], 
                                                    row["Anion chain"] + str(row["Anion id"]), row["Anion symbol"], row["Anion group id"],
                                                    row["Pi acid Code"], 
                                                    row["Pi acid chain"]+str(row["Piacid id"]) ,
                                                    str(row["Distance"])[:3] , str(row["Same semisphere"]) , str(row["Latitude diff"])) 
        
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