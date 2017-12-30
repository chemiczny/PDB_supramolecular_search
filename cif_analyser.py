# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 13:56:41 2017

@author: michal
"""

from Bio.PDB import MMCIFParser, NeighborSearch, Selection
import networkx as nx
import numpy as np

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def getAverageCoords ( allAtomsList, atomsIndList ):
    averageCoords = [ 0. , 0. , 0. ]
    
    for atomInd in atomsIndList:
        new_coords = allAtomsList[atomInd].get_coord()
        
        for coord in range(3):
            averageCoords[coord] += new_coords[coord]
            
    atomsNo = float(len(atomsIndList))
    
    for coord in range(3):
        averageCoords[coord] /= atomsNo
        
    return averageCoords
    
def isFlat(allAtomsList, atomsIndList):
    if len(atomsIndList) <= 3:
        print("LOL strasznie krotki pierscien")
        return True
        
    A = np.array( allAtomsList[ atomsIndList[0] ].get_coord() )
    B = np.array( allAtomsList[ atomsIndList[1] ].get_coord() )
    C = np.array( allAtomsList[ atomsIndList[2] ].get_coord() )
    
    vec1 = A-B
    vec2 = B-C
    
    np.cross(vec1, vec2)
    norm_vec = normalize(np.cross(vec1, vec2))
    
    lastAtom = C    
    for i in range( 3, len(atomsIndList)  ):
        atomInd = atomsIndList[i]
        D = np.array( allAtomsList[ atomInd ].get_coord() )
        new_vec = normalize( lastAtom - D )
        
        if abs( np.inner( new_vec, norm_vec ) ) > 0.09:
            return False
            
    
    return True

def getRingsCentroids( molecule ):
    atoms = list(molecule.get_atoms())
    G = nx.Graph()
    
    for atom1Ind in range(len(atoms)):
        atom1 = atoms[atom1Ind]
        if atom1.get_name() == "H":
            continue
        
        for atom2Ind in range(atom1Ind+1, len(atoms)):
            atom2 = atoms[atom2Ind]            
            if atom2.get_name() == "H":
                continue
            
            distance = atom1 - atom2
            
            if distance < 1.8 :
                G.add_edge(atom1Ind, atom2Ind)
                
    cycles = list(nx.cycle_basis(G))
    centroids = []
    for cycle in cycles:
        if len(cycle) > 6:
            continue
        
        if not isFlat(atoms, cycle):
            continue
        
        centroids.append(getAverageCoords( atoms, cycle))
        
    return centroids
        

def findSupramolecularAnionPiLigand( ligandCode, cifFile ):
    parser = MMCIFParser()
    structure = parser.get_structure('temp', cifFile)
    atoms = Selection.unfold_entities(structure, 'A')    
    
    ligands = []
    for residue in structure.get_residues():
        if ligandCode == residue.get_resname():
            print("Znalazlem!", ligandCode)
            ligands.append(residue)
    
    for ligand in ligands:
        centroids = getRingsCentroids( ligand )
        #print(centroids)
        
        ns = NeighborSearch(atoms)
        for centroid in centroids:
            neighbors = ns.search(np.array(centroid), 5, 'A')
            extractNeighbours( neighbors, ligandCode )
            
def extractNeighbours( atomList, ligandCode ):
    
    for atom in atomList:
        element_symbol = atom.element
        
        if not isItWorthAnalyzing(atom, ligandCode):
            continue
        
        potentiallySupramolecular = False
        
        if element_symbol == "O":
            potentiallySupramolecular = handleOxygen(atom)
        elif element_symbol in [ "F", "CL", "BR", "I" ]:
            potentiallySupramolecular = handleHalogens(atom)
        elif element_symbol in [ "S", "Se" ]:
            potentiallySupramolecular = handleChalcogens(atom)
        elif element_symbol in [ "N", "P" ]:
            potentiallySupramolecular = handlePnictogens(atom)
        else:
            potentiallySupramolecular = handleOthers(atom)            
            
        if potentiallySupramolecular:
            print("Znaleziono cos do obliczen")
            
            print(atom.get_parent().get_resname(), atom.get_fullname(), atom.element)

def isItWorthAnalyzing(atom, ligandCode):
    parent_name = atom.get_parent().get_resname()
    
    bad_parents = [ "HOH", "ALA", "ARG", "ASN", "CYS", "GLN", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
        
    if parent_name == ligandCode:
        return False
    
    if parent_name in bad_parents:
        return False
        
    return True

def handleOxygen( atom ):
    atoms = list(atom.get_parent().get_atoms())
    graph, oxygenInd = molecule2graph(atom, atoms )
    
    oxygen_neighbors = []
    oxygen_neighbors = graph.neighbors(oxygenInd)
    
    if len(oxygen_neighbors) > 1:
        print("Prawdopodobnie ester lub eter")
        return
        
    oxygen_neighbor_index = oxygen_neighbors[0]
    oxygen_neighbor_symbol = atoms[ oxygen_neighbor_index ].element
    
    print("Moj somsiad: ", oxygen_neighbor_symbol)
    
    center_neighbors = graph.neighbors( oxygen_neighbor_index )
    
    oxygens_found = 0
    for unknown_atom in center_neighbors:
        if atoms[unknown_atom].element == "O":
            oxygens_found+=1
            
    if oxygen_neighbor_symbol == "C" and oxygens_found < 2:
        return False
    
    return True
    
    

def handleHalogens( atom):
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) > 1:
        return False
    
    return True

def handleChalcogens(atom):
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) > 1:
        return False
    
    return True

def handlePnictogens(atom):
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) != 2:
        return False
    
    is_carbon = False
    
    for unknown in atoms:
        if unknown.element == "C":
            is_carbon = True
            
    return is_carbon

def handleOthers(atom):
    return False

def molecule2graph( atom, atoms ):
    atomName = atom.get_fullname()    
    
    G = nx.Graph()
    atoms_found = []
    for atom1Ind in range(len(atoms)):
        atom1 = atoms[atom1Ind]
        if atom1.element == "H":
            continue
        
        if atom1.get_fullname() == atomName:
            atoms_found.append(atom1Ind)
        
        for atom2Ind in range(atom1Ind+1, len(atoms)):
            atom2 = atoms[atom2Ind]            
            if atom2.element == "H":
                continue
            
            distance = atom1 - atom2
            
            if distance < 1.8 :
                G.add_edge(atom1Ind, atom2Ind)
                
    if len(atoms_found) != 1:
        print("WTF!? ", atoms_found)
        
    return G, atoms_found[0]

if __name__ == "__main__":
    findSupramolecularAnionPiLigand( "LUM", "cif/1he5.cif" )
    findSupramolecularAnionPiLigand( "NAP", "cif/3bcj.cif" )
    findSupramolecularAnionPiLigand( "NAP", "cif/4lbs.cif" )
#findSupramolecular( "LUM", 666 )