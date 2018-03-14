# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 20:19:28 2018

@author: michal
"""

import sys
babel_path = "/home/michal/Tools//openbabel-install/lib/python3.5/site-packages"
if not babel_path in sys.path:
    sys.path.append(babel_path)
import pybel
import networkx as nx
import numpy as np
import math

def sdf_get_molecule( sdf_file, molecule_code ):
    sdf = pybel.readfile("sdf", sdf_file )
    
    molecules = []
    for molecule in sdf:
        if molecule_code == molecule.data["field_0"]:
            molecules.append(molecule)
            
    return molecules
    

def pybel_molecule2graph(molecule):  
    atoms = molecule.atoms     
    thresholds = { "C" : 1.8, "O" : 1.8, "N" : 1.8, "S" : 2.2,
                  "F" : 1.6, "CL" : 2.0, "BR" : 2.1, "I" : 2.2 }
                  
    elements = [ None,
        "H", "He",
        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]
    
    G = nx.Graph()

    for atom1Ind, atom1 in enumerate(atoms):
        threshold1 = 2.2
        
        element1 = elements[ atom1.atomicnum ]
        if element1 == "H":
            continue
        
        if element1 in thresholds.keys():
            threshold1 = thresholds[element1]
        
        
        for atom2Ind, atom2 in enumerate(atoms[atom1Ind+1:], atom1Ind+1):
            threshold2 = 2.2      
            element2 = elements[ atom2.atomicnum ]
            if element2 == "H":
                continue
            
            if element2 in thresholds.keys():
                threshold2 = thresholds[element2]
            
            distance = 0
            dist_vec = np.array(atom1.coords) - np.array(atom2.coords)
            for element in dist_vec:
                distance += element*element
                
            dist = math.sqrt(distance)
            
            threshold = max( threshold1, threshold2 )
            if dist < threshold :
                print(dist)
                G.add_edge(atom1Ind, atom2Ind)
                G.node[atom1Ind]["element"] = element1
                G.node[atom2Ind]["element"] = element2
        
    return G

if __name__ == "__main__":
    ligprepFile = "sdf/ligprep_2-out_cutted.sdf"
    
    molecules = sdf_get_molecule(ligprepFile, "PYR")
    print("znalazlem: ", len(molecules))
    G = pybel_molecule2graph(molecules[0])
    for node in G.nodes():
        print(G.node[node]["element"])
