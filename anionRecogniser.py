#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 14:04:23 2018

@author: michal
"""
from ringDetection import molecule2graph

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
        element_symbol = atom.element
        
        if not isItWorthAnalyzing(atom, ligandCode):
            continue
        
        potentiallySupramolecular = False
        anionType = ""
        
        if element_symbol == "O":
            potentiallySupramolecular, anionType = handleOxygen(atom)
        elif element_symbol in [ "F", "CL", "BR", "I" ]:
            potentiallySupramolecular, anionType = handleHalogens(atom)
        elif element_symbol in [ "S", "SE" ]:
            potentiallySupramolecular, anionType = handleChalcogens(atom, ns)
        elif element_symbol in [ "N", "P" ]:
            potentiallySupramolecular, anionType = handlePnictogens(atom)
        else:
            potentiallySupramolecular, anionType = handleOthers(atom)            
            
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

def handleOxygen( atom ):
    """
    Zweryfikuj atom tlenu jako potencjalny anion.
    W tej chwili jako aniony identyfikujemy wszystkie atomy tlenu, ktore:
    -nie naleza do czasteczki wody
    -jest sa zwiazane z atomem C, to atom C jest zwiazany z co najmniej dwoma
        atomami O (wykluczamy alkohole, etery, aldehydy, ketony...)
        
    TODO: a estry? a C-koniec bialka?
    
    Wejscie:
    atom - obiekt Atom (Biopython), tlen
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    graph, oxygenInd = molecule2graph( atoms, atom )
    
    oxygen_neighbors = []
    oxygen_neighbors = list(graph.neighbors(oxygenInd))
    
    if len(oxygen_neighbors) == 0:
        return True, "Complex-O?"
    elif len(oxygen_neighbors) > 1:
#        print("Prawdopodobnie ester lub eter")
        return False, "O"
        
    oxygen_neighbor_index = oxygen_neighbors[0]
    oxygen_neighbor_symbol = atoms[ oxygen_neighbor_index ].element
    
    center_neighbors = list(graph.neighbors( oxygen_neighbor_index ))
    
    oxygens_found = 0
    for unknown_atom in center_neighbors:
        if atoms[unknown_atom].element == "O":
            oxygens_found+=1
            
    if oxygen_neighbor_symbol == "C" and oxygens_found < 2:
        return False, oxygen_neighbor_symbol+"O"+str(oxygens_found)
    elif oxygen_neighbor_symbol == "C" and oxygens_found == 2:
        acid = oxygenInCarboxylicGroup( oxygenInd, graph )
        return acid, oxygen_neighbor_symbol+"O"+str(oxygens_found)
    elif oxygen_neighbor_symbol == "N" and oxygens_found == 2:
        isNitro = oxygenInNitroGroup( oxygenInd, graph )
        return not isNitro, oxygen_neighbor_symbol+"O"+str(oxygens_found)
    
    return True, oxygen_neighbor_symbol+"O"+str(oxygens_found)
    
def oxygenInCarboxylicGroup( atomNode , graph ):
    neighbors = list(graph.neighbors(atomNode))

    if len(neighbors) != 1:
        return False
        
    if graph.node[ neighbors[0] ]["element"] != "C" :
        return False
        
    neighborsC = list(graph.neighbors( neighbors[0] ))
    
    if len(neighborsC) < 2:
        return False
        
    for potentiallyOxygen in neighborsC:
        if graph.node[potentiallyOxygen]["element"] == "O":
            oxygenNeigh = list(graph.neighbors( potentiallyOxygen ))
            if len( oxygenNeigh ) != 1:
                return False
                
    return True

def oxygenInNitroGroup( atomNode, graph ):
    neighbors = list(graph.neighbors(atomNode))
    
    if len(neighbors) != 1:
        return False
        
    if graph.node[ neighbors[0] ]["element"] != "N" :
        return False
        
    neighborsN = list(graph.neighbors( neighbors[0] ))
    
    if len(neighborsN) < 3:
        return False
        
    elements = []
    for n in neighborsN:
        elements.append(graph.node[n]["element"])
                
    if elements.count("O") == 2 and elements.count("C") == 1:
        return True 
        
    return False

def handleHalogens( atom):
    """
    Zweryfikuj atom fluorowca jako potencjalny anion.
    W chwili obecnej przyjmujemy, ze atom fluorowca jest anionem
    wtedy gdy jest nie zwiazany z zadna wieksza czasteczka (czyli wykrywamy
    aniony proste: F-, CL-, BR-, I-)
    
    Wejscie:
    atom - obiekt Atom (Biopython), fluorowiec
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) == 1:
        return True, atom.element
        
    graph, halogenInd = molecule2graph( atoms, atom )
    halogen_neighbors = list(graph.neighbors(halogenInd))
    
    if len(halogen_neighbors) == 0:
        return True, "Complex-"+atom.element+"?"
        
    return False, atom.element
    

def handleChalcogens(atom, ns):
    """
    Zweryfikuj pozostale tlenowce (oprocz tlenu) jako potencjalne aniony.
    Podobnie jak w przypadku fluorowcow. W ten sposob ogarniamy tylko S2-
    
    Wejscie:
    atom - obiekt Atom (Biopython), tlenowiec
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) == 1:
        return True, atom.element
        
    graph, chalcogenInd = molecule2graph( atoms, atom )
    chalcogen_neighbors = list(graph.neighbors(chalcogenInd))
    
    if len(chalcogen_neighbors) == 0:
        return True, "Complex-"+atom.element+"?"
    elif len(chalcogen_neighbors) == 1:
        #tiol
        atom_neighbor = atoms[chalcogen_neighbors[0]]
        if atom.element == "S" and atom_neighbor.element == "C":
            neighbors = ns.search( atom.get_coord(), 2.15, 'A')
            neighborsElements = [ neighbor.element.upper() for neighbor in neighbors ]
            carbons = neighborsElements.count("C")
            sulfurs = neighborsElements.count("S")
            if len(neighbors) == 2 and carbons == 1 and sulfurs == 1 :
                return True, "RSH"
            elif len(neighbors) == 3 and carbons == 1 and sulfurs == 2 :
                return True, "S~S"
                
        elif atom.element == "S" and atom_neighbor.element == "S":  
            #mostek disulfidowy
            return True, "S~S?"
        
            
    elif len(chalcogen_neighbors) == 2:
        #mostek disulfidowy
        neighbors_element = []
        neighbors_element.append( atoms[chalcogen_neighbors[0]].element )
        neighbors_element.append( atoms[chalcogen_neighbors[1]].element )
        
        if atom.element == "S" and "S" in neighbors_element and "C" in neighbors_element:
            return True, "RS~SR"
        elif atom.element == "S" and "C" in neighbors_element and "N" in neighbors_element:
            return True, "SCN"
            
    return False, atom.element

def handlePnictogens(atom):
    """
    Zweryfikuj atom azotowca jako potencjalny anion.
    W chwili obecnej rozwazamy tylko jony CN- oraz N3-, NCS-
    
    Wejscie:
    atom - obiekt Atom (Biopython), azotowiec
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) == 2:
    
        is_carbon = False
        
        for unknown in atoms:
            if unknown.element == "C":
                is_carbon = True
                
        return is_carbon, "CN"
        
    elif len(atoms) ==3:
        onlyNitrogen = True
        sulphfurExists = False
        carbonExists = False
        
        for unknown in atoms:
            if unknown.element != "N":
                onlyNitrogen = False
            elif unknown.element == "S":
                sulphfurExists = True
            elif unknown.element == "C":
                carbonExists = True
                
        if onlyNitrogen:
            return True, "N3"
        elif carbonExists and sulphfurExists:
            return True, "NCS"
    
    return False, atom.element

def handleOthers(atom):
    """
    Zweryfikuj niezakwalifikowany atom jako potencjalny anion.
    Obecnie nie wiemy jakie kryteria tu wrzucic wiec ta procedura zawsze zwraca
    False
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    return False, atom.element