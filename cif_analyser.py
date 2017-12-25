# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 13:56:41 2017

@author: michal
"""

from Bio.PDB import *
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
        

def findSupramolecular( ligandCode, cifFile ):
    parser = MMCIFParser()
    structure = parser.get_structure('temp', cifFile)
    
    ligands = []
    for residue in structure.get_residues():
        if ligandCode == residue.get_resname():
            print("Znalazlem!")
            ligands.append(residue)
    
    for ligand in ligands:
        centroids = getRingsCentroids( ligand )
        print(centroids)
    

if __name__ == "__main__":
    findSupramolecular( "LUM", "1he5.cif" )
#findSupramolecular( "LUM", 666 )