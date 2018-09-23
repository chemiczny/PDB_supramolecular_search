#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 14:04:23 2018

@author: michal
"""
#fast ang ugly code especially for prometheus
import sys
from os.path import isdir
if isdir("/net/people/plgglanow/pythonPackages") and not "/net/people/plgglanow/pythonPackages" in sys.path :
    sys.path.insert(0, "/net/people/plgglanow/pythonPackages" )
    
from Bio.PDB import Selection
from ringDetection import getSubstituents, isFlat, isFlatPrimitive, findInGraph, moleculeFragment2graph
from anionTemplateCreator import anionMatcher
import json
from os.path import join
from glob import glob
from networkx.readwrite.json_graph import node_link_graph
from copy import copy
from biopythonUtilities import createResId, createResIdFromAtom

def extractAnionAtoms( atomList, ligand, ns ):
    """
    Wydziel atomy, ktore moga byc anionami w sasiedztwie liganda
    
    Wejscie:
    atomList - lista obiektow Atom (Biopython)
    ligandCode - kod liganda
    
    Wysjcie:
    extractedAtoms - lista slownikow z informacjami o potencjalnym anionie,
        klucze: Atom - atom, AnionType - rodzaj anionu
    """
#    print("wyszukiwanie anionow start ", ligand.get_resname(), ligand.get_id())
    extractedAtoms = []
#    print("wyciagam somsiasow")
    residuesNeighbor = Selection.unfold_entities(atomList, 'R')  
#    print("wyciagnente")
    resId2atoms, resId2res = createResDicts(atomList, residuesNeighbor, ligand)
#    print("atomy przetransformowane")
    aminoacids = [ "ALA", "ARG", "ASN", "CYS", "GLN", "GLY", "HIS", "GLU", "ASP",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
    
    for resId in resId2atoms:
        resAtoms = resId2atoms[resId]
        elements = atomsList2ElementsList( resAtoms )
        resName = resId2res[ resId ].get_resname()
        
        if resName.upper() in aminoacids and not ( "O" in elements or "S" in elements ):
            continue
#        print("Wyciagam polaczenia")
        atoms = getResidueWithConnections( resAtoms, ns )
#        initGraph = molecule2graph(atoms)
#        print("Polaczenia wyciagniete")
        for atom in resAtoms:
            potentiallySupramolecular = False
            anionType = ""
            
            potentiallySupramolecular, anionType = searchInAnionTemplates(atom, atoms)            
                
            if potentiallySupramolecular:
    #            print("cos watergo uwagi :D", atom.get_fullname())
                extractedAtoms.append({ "Atom" : atom, "AnionType" : anionType})
            
#    print("Wyszukiwanie anionow stop")
    return extractedAtoms

def getResidueWithConnections( atoms, ns ):
    atomParent = atoms[0].get_parent()
    
    parentAtoms = list(atomParent.get_atoms())
    neighbors = []
    for atom in atoms:
        neighbors += ns.search( atom.get_coord(), 2 , 'A')
#    neighbors = set(neighbors)
    neighbors = list(set(neighbors))
    newNeighbors = []
    for atom in neighbors:
#        if atom.element != "C" :
        newNeighbors += ns.search( atom.get_coord(), 2 , 'A')
        
    neighbors = list(set(newNeighbors))
    
    newNeighbors = []
    for atom in neighbors:
        if atom.element != "C":
            newNeighbors += ns.search( atom.get_coord(), 2 , 'A')
        
    neighbors = list(set(newNeighbors))
    parent = atomParent.get_parent()
    for neighbor in neighbors:
        if parent != neighbor.get_parent():
            parentAtoms.append(neighbor)
    
    return list(set(parentAtoms))

def atomsList2ElementsList( atomList ):
    elements = set()
    for a in atomList:
        elements.add(a.element)
        
    return elements

def createResDicts(atoms, residues, ligand):
    ligandId = createResId(ligand)
    resId2res = {}
    for res in residues:
        resId2res[ createResId(res) ] = res
        
    bad_parents = [ "HOH", "DOD", "OXY"]
    resId2atoms = {}
    
    for a in atoms:
        resId = createResIdFromAtom(a)
        if resId == ligandId:
            continue
        resname = a.get_parent().get_resname()
        
        if resname in bad_parents:
            continue
        if not resId in resId2atoms:
            resId2atoms[resId] = [a]
        else:
            resId2atoms[resId].append(a)
            
    return resId2atoms, resId2res
            

def graph2Composition( graph ):
    composition = {}
    for node in list(graph):
        element = graph.node[node]["element"]
        if not element in composition:
            composition[element] = 1
        else:
            composition[element] += 1
            
    return composition

def searchInAnionTemplates( atom, atoms ):
    if not hasattr(searchInAnionTemplates, "templates"):
        searchInAnionTemplates.templates =  getAllTemplates()
    element = atom.element.upper()
    
    if not element in searchInAnionTemplates.templates:
        return False, element
            
#    print("transformuje na graf")
    graph, atomInd = moleculeFragment2graph( atoms, atom, 5.0 )
#    print("rozmiar grafu: ", len(graph.nodes()))
#    print("Przetransformowalem")
    composition = graph2Composition(graph)
#    print("Kompozycja zrobiena")
    
    priority2template = searchInAnionTemplates.templates[element]
#    print("templaty dla atomu znaleziene")
    for priority in sorted(priority2template.keys()):
        for guess in priority2template[priority]:
            if not dummyCompare(copy(composition), guess):
                 continue
#            print("maczuje")
            matchResult, anionGroup = try2matchTemplate(graph, atomInd, guess, atoms)
#            print("pomaczowane")
            if matchResult:
                return matchResult, anionGroup
            
    return False, element

def searchInAnionTemplatesBis( atom, atoms, initGraph ):
    if not hasattr(searchInAnionTemplates, "templates"):
        searchInAnionTemplates.templates =  getAllTemplates()
    element = atom.element.upper()
    
    if not element in searchInAnionTemplates.templates:
        return False, element
            
    graph, atomInd = findInGraph( initGraph, atom, atoms )
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
#    print("Tworze obiekta")
    anMatcher = anionMatcher(moleculeGraph, graphTemplate)
#    print("Potworzylem")
    
    if graphTemplate.graph["fullIsomorphism"]:
#        print("pelny izo")
        result = anMatcher.is_isomorphic()
    else:
#        print("czastkowy izo")
        result = anMatcher.subgraph_is_isomorphic()
#    print("mom rezultat")
    if not result:
        moleculeGraph.node[atomId]["charged"] = False
        return result, moleculeGraph.node[atomId]["element"]
    
    matching = anMatcher.mapping
    
    if graphTemplate.graph["geometry"] == "planarWithSubstituents":
        flatAnalysis = isFlat(atoms, list(matching.keys()), getSubstituents( moleculeGraph, list(matching.keys()) ) )
        if not flatAnalysis["isFlat"]:
            moleculeGraph.node[atomId]["charged"] = False
            return False, moleculeGraph.node[atomId]["element"]
    elif graphTemplate.graph["geometry"] == "planar":
#        print("Sprawdzam płaskosc")
        flatAnalysis = isFlatPrimitive(atoms, list(matching.keys() ))
        if not flatAnalysis["isFlat"]:
            print("nie je plaski", graphTemplate.graph["name"])
            moleculeGraph.node[atomId]["charged"] = False
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
    
    moleculeGraph.node[atomId]["charged"] = False
    return True, anionGroup

def getElementFromMatch( matching, node, graph):
    for match in matching:
        if matching[match] == node:
            return graph.node[match]["element"]

if __name__ == "__main__":
    pass