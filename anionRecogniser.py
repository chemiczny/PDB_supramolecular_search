#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 14:04:23 2018

@author: michal
"""
#fast ang ugly code especially for prometheus
import sys
from os.path import isdir

if isdir("/net/archive/groups/plggsuprm/pythonPackages") and not "/net/archive/groups/plggsuprm/pythonPackages" in sys.path :
    sys.path.insert(0, "/net/archive/groups/plggsuprm/pythonPackages" )
    
from Bio.PDB import Selection, NeighborSearch
from ringDetection import getSubstituents, isFlat, isFlatPrimitive, findInGraph, molecule2graph
from anionTemplateCreator import anionMatcher
import json
from os.path import join
from glob import glob
from networkx.readwrite.json_graph import node_link_graph
from copy import copy
from biopythonUtilities import createResId, createResIdFromAtom
#from supramolecularLogging import writeAdditionalInfo
#from time import time

class AnionData:
    def __init__(self, anionType, charged, anionId):
        self.anionType = anionType
        self.charged = charged
        self.anionId = anionId
        self.hBondsAnalyzed = False
        
        
class Property:
    def __init__(self, kind, atomsInvolved):
        self.kind = kind
        self.atomsInvolved = atomsInvolved


class AnionRecogniser:
    def __init__(self):
        self.freeAnionId = 0
        self.templates = getAllTemplates()
        
        self.properties2calculatePack = {}
        
    def cleanPropertiesPack(self):
        self.properties2calculatePack = {}

    def extractAnionAtoms(self, atomListOrig, ligand, ns ):
        """
        Wydziel atomy, ktore moga byc anionami w sasiedztwie liganda
        
        Wejscie:
        atomList - lista obiektow Atom (Biopython)
        ligandCode - kod liganda
        
        Wysjcie:
        extractedAtoms - lista slownikow z informacjami o potencjalnym anionie,
            klucze: Atom - atom, AnionType - rodzaj anionu
        """
    #    print("ekstrackja start")
    #    print("wyszukiwanie anionow start ", ligand.get_resname(), ligand.get_id())
        extractedAtoms = []
        atomList = []
        
        for atom in atomListOrig:
            if atom.get_parent() == ligand:
                continue
            if hasattr(atom, "anionData"):
                if atom.anionData.charged:
                    extractedAtoms.append({ "Atom" : atom, "AnionType" : atom.anionData.anionType, "AnionId" : atom.anionData.anionId})
            else:
                atomList.append(atom)
                
    #    print("wyciagam somsiasow")
        residuesNeighbor = Selection.unfold_entities(atomList, 'R')  
    #    print("wyciagnente")
    #    dictStart = time()
        resId2atoms, resId2res = createResDicts(atomList, residuesNeighbor, ligand)
    #    dictTime = time() - dictStart
    #    print("atomy przetransformowane")
        aminoacids = [ "ALA", "ARG", "ASN", "CYS", "GLN", "GLY", "HIS", "GLU", "ASP",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
        
    #    findConnectionsTime = 0
    #    graphSearch = 0
    #    loopStart = time()
        for resId in resId2atoms:
            resAtoms = resId2atoms[resId]
            elements = atomsList2ElementsList( resAtoms )
            resName = resId2res[ resId ].get_resname()
            
            if resName.upper() in aminoacids and not ( "O" in elements or "S" in elements ):
                continue
    #        print("Wyciagam polaczenia")
    #        getResStart = time()
            atoms = getResidueWithConnections( resAtoms, ns )                                        
            
            nsRes = NeighborSearch(atoms)
    #        findConnectionsTime += time() - getResStart
    #        initGraph = molecule2graph(atoms)
    #        print("Polaczenia wyciagniete")
            for atom in resAtoms:
                potentiallySupramolecular = False
                anionType = ""
    #            startSearch = time()
                potentiallySupramolecular, anionType = self.searchInAnionTemplates(atom, atoms, nsRes)            
    #            graphSearch += time() - startSearch
                    
                if potentiallySupramolecular:
        #            print("cos watergo uwagi :D", atom.get_fullname())
                    extractedAtoms.append({ "Atom" : atom, "AnionType" : atom.anionData.anionType, "AnionId" : atom.anionData.anionId})
                
    #    print("Wyszukiwanie anionow stop")
    #    loopTime = time() - loopStart
    #    print("preparacja :", dictTime)
    #    print("petelka: ", loopTime)
    #    
    #    print("wycinanie reziduum ", findConnectionsTime)
    #    print("przeszukiwanie grafu", graphSearch)
        return extractedAtoms
    
    def searchInAnionTemplates(self,  atom, atoms, ns ):
    #    fullTime = 0
    #    graphConversion = 0
    #    matchingTemplate = 0
    #    
    #    start = time()
        if hasattr(atom, "anionData"):
            if atom.anionData.charged:
                return True, atom.anionData.anionType
            else:
                return False, atom.anionData.anionType
                
        element = atom.element.upper()
        
        
        if not element in self.templates:
            return False, element
                
    #    print("transformuje na graf")
    #    graphStart = time()
        atoms5 = ns.search(atom.get_coord(), 5.0, 'A')
        graph, atomInd = molecule2graph( atoms5, atom )
    #    graphConversion += time() - graphStart
    #    print("rozmiar grafu: ", len(graph.nodes()))
    #    print("Przetransformowalem")
        composition = graph2Composition(graph)
    #    print("Kompozycja zrobiena")
        
        priority2template = self.templates[element]
    #    print("templaty dla atomu znaleziene")
        for priority in sorted(priority2template.keys()):
            for guess in priority2template[priority]:
                if not dummyCompare(copy(composition), guess):
                     continue
                 
                matchResult, anionGroup = self.try2matchTemplate(graph, atomInd, guess, atoms5)
                if matchResult:
                    return matchResult, anionGroup
                
        return False, element
    
    def try2matchTemplate(self, moleculeGraph, atomId, graphTemplate, atoms):
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
        reverseMapping = {}
        for key in matching:
            reverseMapping[ matching[key] ] = key
        
        if graphTemplate.graph["geometry"] == "planarWithSubstituents":
            flatAnalysis = isFlat(atoms, list(matching.keys()), getSubstituents( moleculeGraph, list(matching.keys()) ) )
            if not flatAnalysis["isFlat"]:
                moleculeGraph.node[atomId]["charged"] = False
                return False, moleculeGraph.node[atomId]["element"]
        elif graphTemplate.graph["geometry"] == "planar":
    #        print("Sprawdzam płaskosc")
            flatAnalysis = isFlatPrimitive(atoms, list(matching.keys() ), 0.5)
            if not flatAnalysis["isFlat"]:
    #            print("nie je plaski", graphTemplate.graph["name"])
                moleculeGraph.node[atomId]["charged"] = False
                return False, moleculeGraph.node[atomId]["element"]
    #        else:
    #            print("jest plaski", graphTemplate.graph["name"])
        
        anionGroup = graphTemplate.graph["name"]
        
        if "X" in graphTemplate.graph["name"] and  graphTemplate.graph["nameMapping"]:
            matching = anMatcher.mapping
            for node in graphTemplate.graph["nameMapping"]:
#                element =  moleculeGraph.node[reverseMapping[node]]["element"]
                element = getElementFromMatch( matching, int(node), moleculeGraph)
                anionGroup = anionGroup.replace( "X", element )
                break
        
        moleculeGraph.node[atomId]["charged"] = False
        
        for aId in matching:
            templateAtomId = matching[aId]
            
            if templateAtomId in graphTemplate.graph["otherCharges"]:
                atoms[aId].anionData = AnionData(anionGroup, True, self.freeAnionId)
            else:
                atoms[aId].anionData = AnionData(anionGroup, False, self.freeAnionId)
                
        atoms[atomId].anionData = AnionData(anionGroup, True, self.freeAnionId)
                        
        if graphTemplate.graph["properties2measure"] and not self.freeAnionId in self.properties2calculatePack:
            self.properties2calculatePack[ self.freeAnionId ] = []
            
            for property2calc in graphTemplate.graph["properties2measure"]:
                kind = list(property2calc.keys())[0]
                
                atomsInvolvedIndexes = property2calc[kind]
                atomsInvolved = []
                
                for aInd in atomsInvolvedIndexes:
                    atomsInvolved.append( atoms[reverseMapping[aInd]] )
                    
                self.properties2calculatePack[ self.freeAnionId ].append( Property(kind, atomsInvolved) )
        
        self.freeAnionId += 1
        return True, anionGroup

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


#
#def searchInAnionTemplatesBis( atom, atoms, initGraph ):
#    if not hasattr(searchInAnionTemplates, "templates"):
#        searchInAnionTemplates.templates =  getAllTemplates()
#    element = atom.element.upper()
#    
#    if not element in searchInAnionTemplates.templates:
#        return False, element
#            
#    graph, atomInd = findInGraph( initGraph, atom, atoms )
#    composition = graph2Composition(graph)
#    
#    priority2template = searchInAnionTemplates.templates[element]
#
#    for priority in sorted(priority2template.keys()):
#        for guess in priority2template[priority]:
#            if not dummyCompare(copy(composition), guess):
#                 continue
#            matchResult, anionGroup = try2matchTemplate(graph, atomInd, guess, atoms)
#            if matchResult:
#                return matchResult, anionGroup
#            
#    return False, element


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

def getElementFromMatch( matching, node, graph):
    for match in matching:
        if matching[match] == node:
            return graph.node[match]["element"]
        

if __name__ == "__main__":
    pass