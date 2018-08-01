#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 13:58:24 2018

@author: michal
"""
import numpy as np
from numpy_utilities import normalize
import networkx as nx
from supramolecularLogging import writeAdditionalInfo
from multiprocessing import current_process

def getAverageCoords ( allAtomsList, atomsIndList ):
    """
    Oblicz srednie wspolrzedne x, y, z wybranych atomow
    
    Wejscie:
    allAtomsList - lista obiektow Atom (cala czasteczka)
    atomsIndList - lista indeksow atomow, ktorych usrednione polozenie
                    ma byc obliczone
                    
    Wyjscie:
    averageCoords - 3-elementowa lista zawierajaca usrednione wspolrzedne 
                    x, y, z wybranych atomow
    """
    averageCoords = [ 0. , 0. , 0. ]
    
    for atomInd in atomsIndList:
        new_coords = allAtomsList[atomInd].get_coord()
        
        for coord in range(3):
            averageCoords[coord] += new_coords[coord]
            
    atomsNo = float(len(atomsIndList))
    
    for coord in range(3):
        averageCoords[coord] /= atomsNo
        
    return averageCoords
    
def isFlat(allAtomsList, atomsIndList, substituents):
    """
    Sprawdz czy wybrane atomy leza w jednej plaszczyznie.
    Procedura: na podstawie polozen trzech pierwszych atomow wyznacza sie 
    wektor normalnych do wyznaczonej przez nich plaszczyzny. Nastepnie 
    sprawdzane jest czy kolejne wiazania tworza wektory prostopadle
    do wektora normalnego. Dopuszczalne jest odchylenie 5 stopni.
    
    TODO: Nie lepiej byloby obliczyc tensor momentu bezwladnosci, 
    zdiagonalizowac go i rozstrzygnac na podstawie jego wartosci wlasnych?
    
    Wejscie:
    allAtomsList - lista obiektow Atom (cala czasteczka)
    atomsIndList - lista indeksow atomow, ktore maja byc zweryfikowane
                    pod wzgledem lezenia w jednej plaszczyznie
                    
    Wyjscie:
    verdict - slownik, posiada klucze: isFlat (zmienna logiczna, True jesli
                struktura jest plaska), normVec (3-elementowa lista float,
                wspolrzedne wektora normalnego plaszczyzny, jesli struktura nie
                jest plaska wszystkie jego wspolrzedne sa rowne 0)
    """
        
    verdict = { 'isFlat' : False, 'normVec' : [ 0, 0, 0] }
    
    if len(atomsIndList) < 3:
        print("Znaleziono 2-wu elementowy cykl!!!")
        return verdict
    
    norm_vec = getNormVec( allAtomsList, atomsIndList)
    
    if len(atomsIndList) <= 3:
        verdict['normVec'] = norm_vec
        return verdict    
    
    centroid = getAverageCoords(allAtomsList, atomsIndList) 
    for i in range(1, len(atomsIndList)  ):
        atomInd = atomsIndList[i]
        D = np.array(allAtomsList[ atomInd ].get_coord() )
        new_vec = normalize( centroid - D )
        
        if abs( np.inner( new_vec, norm_vec ) ) > 0.09:
            return verdict
            
    for substituentKey in substituents:
        D = np.array( allAtomsList[ substituentKey ].get_coord() )
        E = np.array( allAtomsList[ substituents[substituentKey ] ].get_coord() )
        new_vec = normalize( E - D )
        
        if abs( np.inner( new_vec, norm_vec ) ) > 0.20:
            return verdict
            
    verdict['isFlat'] = True
    verdict['normVec'] = norm_vec
    verdict['coords'] = centroid
    return verdict
    
def isFlatPrimitive(allAtomsList, atomsIndList, maxDist = 0.15):
    """
    Sprawdz czy wybrane atomy leza w jednej plaszczyznie.
    Procedura: na podstawie polozen trzech pierwszych atomow wyznacza sie 
    wektor normalnych do wyznaczonej przez nich plaszczyzny. Nastepnie 
    sprawdzane jest czy kolejne wiazania tworza wektory prostopadle
    do wektora normalnego. Dopuszczalne jest odchylenie 5 stopni.
    
    TODO: Nie lepiej byloby obliczyc tensor momentu bezwladnosci, 
    zdiagonalizowac go i rozstrzygnac na podstawie jego wartosci wlasnych?
    
    Wejscie:
    allAtomsList - lista obiektow Atom (cala czasteczka)
    atomsIndList - lista indeksow atomow, ktore maja byc zweryfikowane
                    pod wzgledem lezenia w jednej plaszczyznie
                    
    Wyjscie:
    verdict - slownik, posiada klucze: isFlat (zmienna logiczna, True jesli
                struktura jest plaska), normVec (3-elementowa lista float,
                wspolrzedne wektora normalnego plaszczyzny, jesli struktura nie
                jest plaska wszystkie jego wspolrzedne sa rowne 0)
    """
        
    verdict = { 'isFlat' : False, 'normVec' : [ 0, 0, 0] }
    
    if len(atomsIndList) < 3:
        print("Znaleziono 2-wu elementowy cykl!!!")
        return verdict
    
    norm_vec = getNormVec( allAtomsList, atomsIndList)
    
    if len(atomsIndList) <= 3:
        verdict['normVec'] = norm_vec
        return verdict    
    
    centroid = getAverageCoords(allAtomsList, atomsIndList) 
    D = -np.inner(centroid, norm_vec)
    
    for atomInd in atomsIndList:
        atomCoord = np.array(allAtomsList[ atomInd ].get_coord() )
        atomDist = abs( np.inner(norm_vec, atomCoord) + D  )
        
        if atomDist > maxDist:
            return verdict
            
    verdict['isFlat'] = True
    verdict['normVec'] = norm_vec
    verdict['coords'] = centroid
    return verdict

def getNormVec( allAtomsList, atomsIndList ):
    norm_vec = np.array([0. , 0. , 0.])
    expanded_list = atomsIndList + atomsIndList[:2]
    for i in range( len(expanded_list) - 2 ):
        A = np.array( allAtomsList[ expanded_list[i] ].get_coord() )
        B = np.array( allAtomsList[ expanded_list[i+1] ].get_coord() )
        C = np.array( allAtomsList[ expanded_list[i+2] ].get_coord() )
        
        vec1 = A-B
        vec2 = B-C
        
        norm_vec += normalize(np.cross(vec1, vec2))
        
    return normalize(norm_vec)

def getRingsCentroids( molecule, returnGraph = False ):
    """
    Znajdz pierscienie w czasteczce i wyznacz wspolrzedne ich srodkow jesli
    sa one aromatyczne.
    
    Procedura: 
    - utworz graf na podstawie danych o atomach, wszystkie atomy lezace blizej 
        niz 1.8 A sa traktowane jako wierzcholki grafu polaczone krawedzia
    - znajdz fundamentalne cykle w grafie (Stosowany algorytm:
    Paton, K. An algorithm for finding a fundamental set of cycles of a graph. 
    Comm. ACM 12, 9 (Sept 1969), 514-518.)
    - odrzuc cykle, ktore skladaja sie z wiecej niz 6 wierzcholkow
    - sprawdz czy znalezione pierscienie sa plaskie
    
    Wejscie:
    -molecule - obiekt Residue (Biopython)
    
    Wyjscie:
    -centroids - lista slownikow z danymi o znalezionych (potencjalnie) 
        aromatycznych pierscienieniach. Slownik zawiera klucze:
        coords (wspolrzedne srodka pierscienia), 
        normVec (wektor normalnych plaszczyzny pierscienia)
    """
    atoms = list(molecule.get_atoms())
    G = molecule2graph(atoms)
                
    cycles = list(nx.cycle_basis(G))
#    print("Znalazlem cykli: ", len(cycles))
    centroids = []
    
    for cycleId, cycle in enumerate(cycles):
        if not onlyLighAtomsInCycle(cycle, atoms):
            continue
        
        if len(cycle) > 6:
            continue
        
        substituents = getSubstituents( G, cycle )        
        flatAnalyse = isFlat(atoms, cycle, substituents)
        if not flatAnalyse['isFlat']:
#            print("cykl nie jest plaski")
            continue
#        else:
#            print("cykl jest plaski")
        
        centroids.append({ "coords" : flatAnalyse["coords"], "normVec" : flatAnalyse["normVec"], "ringSize" : len(cycle) ,
                          "cycleAtoms" : cycle, "cycleId" : cycleId })
        
    if returnGraph:
        return centroids, G
    else:
        return centroids

def onlyLighAtomsInCycle( cycle, atoms ):
    for atomInd in cycle:
        if not atoms[atomInd].element in [ "C", "O", "N", "S" ]:
            return False
    return True
    
def getSubstituents( graphMolecule, cycle ):
    substituents = {}
    
    for atom in cycle:
        candidates = graphMolecule.neighbors( atom )
        for candidate in candidates:
            if not candidate in cycle:
                substituents[atom] = candidate
                
    return substituents

def createAtomId(atom):
    atomCoords = atom.get_coord()
    
    return str(atomCoords[0])+str(atomCoords[1])+str(atomCoords[2])

def molecule2graph( atoms, atom = None ):
    """
    Konwersja czasteczki na graf (networkx)
    
    Wejscie:
    atom - obiekt Atom (Biopython), ktorego polozenie w grafie jest istotne
    atoms - lista wszystkich atomow
    
    Wyjscie:
    G, atomInd - graf (networkx), indeks wejsciowego atomu (wierzcholek w grafie)
    """
#        print("szukam ", atomName, parentId, atom.get_parent().get_resname())
        
    thresholds = { "C" : 1.8, "O" : 1.8, "N" : 1.8, "S" : 2.2,
                  "F" : 1.6, "CL" : 2.0, "BR" : 2.1, "I" : 2.2 }
    
    G = nx.Graph()
    atoms_found = []
    for atom1Ind, atom1 in enumerate(atoms):
        threshold1 = 2.2

        if atom1.element == "H":
            continue
        
        G.add_node(atom1Ind, element = atom1.element)
        
        if atom1.element in thresholds.keys():
            threshold1 = thresholds[atom1.element]
        
        if atom != None:
            if atom == atom1:
                atoms_found.append(atom1Ind)
        
        for atom2Ind, atom2 in enumerate(atoms[atom1Ind+1:], atom1Ind+1):
            threshold2 = 2.2
          
            if atom2.element == "H":
                continue
            
            if atom2.element in thresholds.keys():
                threshold2 = thresholds[atom2.element]
            
            distance = atom1 - atom2
            
            threshold = max( threshold1, threshold2 )
            if distance < threshold :
                G.add_edge(atom1Ind, atom2Ind)
                
    if atom != None:
        if len(atoms_found) != 1 :
            fileId = current_process()
            writeAdditionalInfo( "Many atoms in graph found!!!! "+str(atoms_found)  , fileId)
            print("WTF!? ", atoms_found)
            
        if not atoms_found[0] in G.nodes():
            print("Nie znalazlem "+ atom.element+" w grafie!")
            G.add_node( atoms_found[0] )
            
        nodes2stay = [ ]
        for node in G.nodes():
            if nx.has_path(G, node, atoms_found[0]):
                nodes2stay.append(node)
            
#        if atoms_found:
        return G.subgraph(nodes2stay), atoms_found[0]
#        else:
#            return G.subgraph(nodes2stay), None
#        return G, atoms_found[0]
        
    return G

def moleculeFragment2graph( atoms, atom , maxDist ):
    """
    Konwersja czasteczki na graf (networkx)
    
    Wejscie:
    atom - obiekt Atom (Biopython), ktorego polozenie w grafie jest istotne
    atoms - lista wszystkich atomow
    
    Wyjscie:
    G, atomInd - graf (networkx), indeks wejsciowego atomu (wierzcholek w grafie)
    """        
    thresholds = { "C" : 1.8, "O" : 1.8, "N" : 1.8, "S" : 2.2,
                  "F" : 1.6, "CL" : 2.0, "BR" : 2.1, "I" : 2.2 }
    
    G = nx.Graph()
    atoms_found = []
    for atom1Ind, atom1 in enumerate(atoms):
        threshold1 = 2.2

        if atom1.element == "H":
            continue
        
        if atom1 - atom > maxDist:
            continue
        
        G.add_node(atom1Ind, element = atom1.element)
        
        if atom1.element in thresholds.keys():
            threshold1 = thresholds[atom1.element]
        
        if atom != None:
            if atom == atom1:
                atoms_found.append(atom1Ind)
        
        for atom2Ind, atom2 in enumerate(atoms[atom1Ind+1:], atom1Ind+1):
            threshold2 = 2.2
          
            if atom2.element == "H":
                continue
            
            if atom2 - atom > maxDist:
                continue
            
            if atom2.element in thresholds.keys():
                threshold2 = thresholds[atom2.element]
            
            distance = atom1 - atom2
            
            threshold = max( threshold1, threshold2 )
            if distance < threshold :
                G.add_edge(atom1Ind, atom2Ind)
                
    if atom != None:
        if len(atoms_found) != 1 :
            fileId = current_process()
            writeAdditionalInfo( "Many atoms in graph found!!!! "+str(atoms_found)  , fileId)
            print("WTF!? ", atoms_found)
            
        if not atoms_found[0] in G.nodes():
            print("Nie znalazlem "+ atom.element+" w grafie!")
            G.add_node( atoms_found[0] )
            
        nodes2stay = [ ]
        for node in G.nodes():
            if nx.has_path(G, node, atoms_found[0]):
                nodes2stay.append(node)
            
        return G.subgraph(nodes2stay), atoms_found[0]
#        return G, atoms_found[0]
        
    return G

def findInGraph( G, atom, atomList ):
    nodes2stay = [ ]
    atoms_found = []

    for node in G.nodes:
        if  atomList[node]  == atom:
            atoms_found.append(node)
            
    if len(atoms_found) != 1 :
        print("WTF!? ", atoms_found)
        
    if atoms_found:
        for node in G.nodes():
            
            if nx.has_path(G, node, atoms_found[0]):
                nodes2stay.append(node)
        
        return G.subgraph(nodes2stay), atoms_found[0]
    else:
        return G.subgraph(nodes2stay), None