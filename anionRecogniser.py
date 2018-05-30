#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 14:04:23 2018

@author: michal
"""
from ringDetection import molecule2graph, getSubstituents, isFlat, isFlatPrimitive
from anionTemplateCreator import anionMatcher
import json
from os.path import join
from glob import glob
from networkx.readwrite.json_graph import node_link_graph
from copy import copy

def extractNeighbours( atomList, ligandCode, ns ):
    """
    Wydziel atomy, ktore moga byc anionami w sasiedztwie liganda
    
    Wejscie:
    atomList - lista obiektow Atom (Biopython)
    ligandCode - kod liganda
    
    Wysjcie:
    extractedAtoms - lista slownikow z informacjami o potencjalnym anionie,
        klucze: Atom - atom, AnionType - rodzaj anionu
    """
    extractedAtoms = []
    
    for atom in atomList:
        if not isItWorthAnalyzing(atom, ligandCode):
            continue
        
        potentiallySupramolecular = False
        anionType = ""
        
        potentiallySupramolecular, anionType = searchInAnionTemplates(atom, ns)            
            
        if potentiallySupramolecular:
#            print("cos watergo uwagi :D", atom.get_fullname())
            extractedAtoms.append({ "Atom" : atom, "AnionType" : anionType})
            
    return extractedAtoms

def isItWorthAnalyzing(atom, ligandCode):
    """
    Szybka, wstepna weyfikacja atomu jako potencjalnego anionu
    Odrzucane sa te atomy, ktore naleza do niekwasnych aminokwasow,
    czasteczki wody oraz ligandu
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    ligandCode - kod liganda
    
    Wyjscie:
    True/False - czy atom moze byc potencjalnym anionem?
    """
    parent_name = atom.get_parent().get_resname()
    
#    bad_parents = [ "HOH", "ALA", "ARG", "ASN", "CYS", "GLN", "GLY", "HIS",
#        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
        
    bad_parents = [ "HOH", "DOD", "OXY"]
        
    if parent_name == ligandCode:
        return False
    
    if parent_name in bad_parents:
        return False
        
    return True

def graph2Composition( graph ):
    composition = {}
    for node in list(graph):
        element = graph.node[node]["element"]
        if not element in composition:
            composition[element] = 1
        else:
            composition[element] += 1
            
    return composition

def searchInAnionTemplates( atom, ns):
    if not hasattr(searchInAnionTemplates, "templates"):
        searchInAnionTemplates.templates =  getAllTemplates()
    element = atom.element.upper()
    
    if not element in searchInAnionTemplates.templates:
        return False, element
    
    atomParent = atom.get_parent()
    
    aminoacids = [ "ALA", "ARG", "ASN", "CYS", "GLN", "GLY", "HIS", "GLU", "ASP",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
    
    if atomParent.get_resname().upper() in aminoacids and not element in [ "O", "S" ]:
        return False, element
    
    atoms = list(atomParent.get_atoms())
    
    neighbors = ns.search( atom.get_coord(), 2.2 , 'A')
    anionId = atomParent.get_id()[1]
    
    for neighbor in neighbors:
        if neighbor.get_parent().get_id()[1] != anionId:
            atoms.append(neighbor)
            
    graph, atomInd = molecule2graph( atoms, atom )
    composition = graph2Composition(graph)
    
    priority2template = searchInAnionTemplates.templates[element]

    for priority in sorted(priority2template.keys()):
        for guess in priority2template[priority]:
            if not dummyCompare(copy(composition), guess):
                 continue
            matchResult, anionGroup = try2matchTemplate(graph, atomInd, guess, atoms)
            if matchResult:
                return matchResult, anionGroup
            
    return False, element

def dummyCompare(composition, template):
    for node in list(template):
        if not "aliases" in template.node[node]:
            element = template.node[node]["element"]
            if element == "X":
                continue
            
            if not element in composition:
                return False
            else:
                composition[element] -= 1
                if composition[element] < 0:
                    return False
        elif not template.node[node]["aliases"]:
            element = template.node[node]["element"]
            if element == "X":
                continue
            if not element in composition:
                return False
            else:
                composition[element] -= 1
                if composition[element] < 0:
                    return False
    return True

def getAllTemplates():
    templates = join("anion_templates", "*", "*.json")
    anionsTemplates = glob(  templates )
    
    graphTemplates = {}
    
    for template in anionsTemplates:
        jsonF = open(template)
        graphTemplate = node_link_graph(json.load(jsonF))
        jsonF.close()
        
        priority = graphTemplate.graph["priority"]
        element = graphTemplate.node[ graphTemplate.graph["charged"] ]["element"]
        
        if not element in graphTemplates:
            graphTemplates[element] = { priority : [ graphTemplate ] }
        elif not priority in graphTemplates[element]:
            graphTemplates[element][priority] = [ graphTemplate ]
        else:
            graphTemplates[element][priority].append(graphTemplate)
            
    return graphTemplates
    
def try2matchTemplate(moleculeGraph, atomId, graphTemplate, atoms):
    
    moleculeGraph.node[atomId]["charged"] = True
    anMatcher = anionMatcher(moleculeGraph, graphTemplate)
    
    if graphTemplate.graph["fullIsomorphism"]:
        result = anMatcher.is_isomorphic()
    else:
        result = anMatcher.subgraph_is_isomorphic()
        
    if not result:
        return result, moleculeGraph.node[atomId]["element"]
    
    matching = anMatcher.mapping
    
    if graphTemplate.graph["geometry"] == "planarWithSubstituents":
        flatAnalysis = isFlat(atoms, list(matching.keys()), getSubstituents( moleculeGraph, list(matching.keys()) ) )
        if not flatAnalysis["isFlat"]:
            return False, moleculeGraph.node[atomId]["element"]
    elif graphTemplate.graph["geometry"] == "planar":
#        print("Sprawdzam płaskosc")
        flatAnalysis = isFlatPrimitive(atoms, list(matching.keys() ))
        if not flatAnalysis["isFlat"]:
            print("nie je plaski", graphTemplate.graph["name"])
            return False, moleculeGraph.node[atomId]["element"]
#        else:
#            print("jest plaski", graphTemplate.graph["name"])
    
    anionGroup = graphTemplate.graph["name"]
    
    if "X" in graphTemplate.graph["name"] and  graphTemplate.graph["nameMapping"]:
        matching = anMatcher.mapping
        for node in graphTemplate.graph["nameMapping"]:
            element = getElementFromMatch( matching, int(node), moleculeGraph)
            anionGroup = anionGroup.replace( "X", element )
            break
    
    return True, anionGroup

def getElementFromMatch( matching, node, graph):
    for match in matching:
        if matching[match] == node:
            return graph.node[match]["element"]

if __name__ == "__main__":
    pass