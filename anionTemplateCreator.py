#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 12:34:08 2018

@author: michal
"""

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.readwrite.json_graph import node_link_data
from os.path import isdir, join, isfile
from os import mkdir
import json
from glob import glob
import shutil

class anionMatcher(GraphMatcher):
    def semantic_feasibility(self, G1_node, G2_node):
        if "charged" in self.G1.node[G1_node]:
            if self.G1.node[G1_node]["charged"] != self.G2.node[G2_node]["charged"]:
                return False
        elif self.G2.node[G2_node]["charged"]:
            return False
        
        if self.G2.node[G2_node]["terminating"]:
            if len(list(self.G2.neighbors(G2_node))) != len(list(self.G1.neighbors(G1_node))):
                return False
            
        if ( self.G2.node[G2_node]["element"] == "X" or "X" in self.G2.node[G2_node]["aliases"] ) and  not self.G1.node[G1_node]["element"] in self.G2.node[G2_node]["notAliases"] :
            return True
        
        return self.G1.node[G1_node]["element"] == self.G2.node[G2_node]["element"] or self.G1.node[G1_node]["element"] in self.G2.node[G2_node]["aliases"]

def addAtribute( graph, nodes, key ):
    if isinstance(nodes, list):
        for nodeId in nodes:
            graph.node[nodeId][key] = True
            
    else:
        graph.node[nodes][key] = True 
        

def saveAnion( atoms, bonds, charged, name, priority , terminating = [], aliases = {},
              notAliases = {}, geometry = {}, fullIsomorphism = False, nameMapping = {}  ):
    graph = nx.Graph()
    
    for i, el in enumerate(atoms):
        graph.add_node(i, element = el, terminating = False, bonded = False, aliases = [], charged = False )
        
    graph.add_edges_from(bonds)
    
    addAtribute( graph, terminating, "terminating")
        
    for nodeId in aliases:
        graph.node[nodeId]["aliases"] = aliases[nodeId]
    
    for nodeId in notAliases:
        graph.node[nodeId]["notAliases"] = notAliases[nodeId]
            
    if not geometry:
        graph.graph["geometry"]= "no restrictions"
    else:
        graph.graph["geometry"]= geometry
        
    graph.graph["fullIsomorphism"] = fullIsomorphism
    graph.graph["name"] = name
    graph.graph["nameMapping"] = nameMapping
    graph.graph["priority"] = priority
        
    fileName = str(priority)+"_"+name
    
    if isinstance( charged , list ) :
        for nodeId in charged:
            saveAnionJson(graph, fileName, nodeId)
    else:
       saveAnionJson(graph, fileName, charged)
        
def saveAnionJson( graph, fileName, charged):
    mainElement = graph.node[charged]["element"]
    elements = [ mainElement ] 
    
    if "aliases" in graph.node[charged]:
        elements += graph.node[charged]["aliases"]
        graph.node[charged]["aliases"] = []
    
    graph.node[charged]["charged"] = True
    graph.graph["charged"] = charged 
    
    oldName = ""
    nameMapping = False
    if "X" in graph.graph["name"] and charged in graph.graph["nameMapping"]:
        oldName = graph.graph["name"]
        nameMapping = graph.graph["nameMapping"][charged]
        graph.graph["nameMapping"].pop(charged)
        
    
    for element in elements:
        graph.node[charged]["element"] = element
        
        if nameMapping:
            graph.graph["name"] = oldName.replace( nameMapping , element)
        
        dir_path = join("anion_templates", element)
        if not isdir( dir_path ):
            mkdir( dir_path )
            
        path2save = getUniquePath( dir_path , fileName)
        output = open(path2save, 'w')
        
        json.dump(node_link_data(graph), output )
        output.close()
            
    graph.node[charged]["charged"] = False
    
def getUniquePath(dirPath, fileName):
    path2save = join( dirPath , fileName+".json")
    if not isfile(path2save):
        return path2save
    
    similarFiles = glob( join(dirPath, fileName)+"_*.json" )
    if not similarFiles:
        return join( dirPath , fileName+"_0.json")
    
    maxNumber = -1
    for s in similarFiles:
        newNumber = int( s[:-5].split("_")[-1] )
        maxNumber = max(maxNumber, newNumber)
        
    return join( dirPath , fileName+"_"+str(maxNumber+1)+".json")
    
def clearAnionTemplates():
    if isdir("anion_templates"):
        shutil.rmtree("anion_templates")
    mkdir("anion_templates")

if __name__ == "__main__":
    clearAnionTemplates()
        
# atoms, bonds, charged, name, priority, terminating = [], aliases = {}, notAliases = {}, geometry = {}, fullIsomorphism = False
    #OXYGEN
    
#    #RCOOH
    saveAnion( [ "C" , "C", "O", "O" ], [ (0,1), (1,2), (1,3) ], 
              2, "RCOO", 0, terminating = [2, 3], 
              geometry = "planar" )
    
    #ClO, BrO, IO, 
    saveAnion([ "CL",  "O" ], [(0, 1)], 
              1, "XO", 5, fullIsomorphism = True,
              aliases = { 0 : [ "BR", "I" ] }, nameMapping = { 0 : "X"} )
    
    #NO2, ClO2, BRO2, 
    saveAnion([ "N",  "O" , "O" ], [(0, 1), (0,2)], 
              1, "XO2", 10, fullIsomorphism = True, aliases = { 0 : ["CL", "BR"]}, nameMapping = { 0 : "X" })
    
    #NO3, CO3, PO3, SO3, AsO3, BO3, ClO3, BRO3
    saveAnion( ["N", "O", "O", "O"], [(0,1), (0,2), (0,3)],
              1, "XO3", 15, fullIsomorphism = True, 
              aliases = { 0 : [ "C", "P", "B", "S", "AS", "CL", "BR" ] }, nameMapping = { 0 : "X" } )
    
    #PO4, SO4, AsO4, ClO4, BRO4
    saveAnion( ["P", "O", "O", "O", "O"], [(0,1), (0,2), (0,3), (0, 4)],
              1, "XO4", 20, fullIsomorphism = True, 
              aliases = { 0 : [ "S", "AS", "CL", "BR" ] }, nameMapping = { 0 : "X" } )
    
#   Ph-OH
    saveAnion( [ "C" , "C" , "C" , "C" , "C", "C" , "O" ], [(0,1),(1,2), (2,3), (3,4),( 4, 5), (5, 0), (5,6)],
              6, "PhOH", 25, terminating = [6], geometry = "planarWithSubstituents")
    
    #    #RBOOH
    saveAnion( [ "X" , "B", "O", "O" ], [ (0,1), (1,2), (1,3) ], 
              2, "RBOO", 30, terminating = [2, 3], 
              notAliases = {0 : [ "O" ] } )
        
    #COO
    saveAnion( [  "C", "O", "O" ], [ (0,1), (0,2) ], 
              1, "COO", 35, terminating = [1, 2],  )
    
    #R-PO4, R-SO4, R-AsO4
    saveAnion( ["P", "O", "O", "O", "O"], [(0,1), (0,2), (0,3), (0, 4)],
              1, "R-XO4", 45, terminating = [ 1, 2, 3 ] ,
              aliases = { 0 : [ "S", "AS" ] }, nameMapping = { 0 : "X" } )
    
    #R2-PO4, R2-SO4, R2-AsO4
    saveAnion( ["P", "O", "O", "O", "O"], [(0,1), (0,2), (0,3), (0, 4)],
              1, "R2-XO4", 47, terminating = [ 1, 2 ] ,
              aliases = { 0 : [ "S", "AS" ] }, nameMapping = { 0 : "X" } )
    
    #R3-PO4, R3-SO4, R3-AsO4
    saveAnion( ["P", "O", "O", "O", "O"], [(0,1), (0,2), (0,3), (0, 4)],
              1, "R2-XO4", 48, terminating = [ 1 ] ,
              aliases = { 0 : [ "S", "AS" ] }, nameMapping = { 0 : "X" } )
    
    #RAsO3, RPO3, RSO3
    saveAnion( ["P", "O", "O", "O", "C"], [(0,1), (0,2), (0,3), (0, 4)],
              1, "RXO3", 50,  terminating = [1, 2, 3] , 
              aliases = { 0 : [ "S", "AS" ] }, nameMapping = { 0 : "X" } )
    
    #R2AsO2, R2PO2, RRSO2
    saveAnion( ["P", "O", "O", "C", "C"], [(0,1), (0,2), (0,3), (0, 4)],
              1, "R2XO2", 55, terminating = [1, 2],
              aliases = { 0 : [ "S", "AS" ] }, nameMapping = { 0 : "X" } )
    
    #F, CL, BR, I, S
    saveAnion( [ "F" ], [], 0, "X", 55, aliases = { 0 : [ "CL", "BR", "I", "S"] },
              fullIsomorphism = True, nameMapping = { 0 : "X"})
    
    #SCN
    saveAnion([ "S",  "C" , "N" ], [(0, 1), (0,2)], 
              [0,1,2], "SCN", 62, fullIsomorphism = True)
    
#    #RSH
    saveAnion( [ "X" , "S" ], [ (0,1)], 
              1, "RSH", 60, terminating = [1], 
              notAliases = {0 : [ "O" ] } )
    
    
    #N3
    saveAnion([ "N",  "N" , "N" ], [(0, 1), (0,2)], 
              [0,1], "N3", 70, fullIsomorphism = True)
    
    #CN
    saveAnion([  "C" , "N" ], [(0, 1)], 
              [0,1], "CN", 75, fullIsomorphism = True)

    
   #    #RSSR
    saveAnion( [ "X" , "S", "S" ], [ (0,1), (1,2)], 
              1, "RSS", 80 , 
              notAliases = {0 : [ "O" ] } )
