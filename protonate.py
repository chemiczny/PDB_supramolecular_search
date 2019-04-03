#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 15:16:17 2019

@author: michal
"""
from ringDetection import molecule2graph
import networkx as nx
import math
import numpy as np
from numpy_utilities import normalize, rotateVector


class HydrogenAtom(object): 
    def __init__(self, coord): 
        self.coord = coord 
        self.element = "H"
    
    def get_coord(self): 
        return self.coord 
    
    def __sub__(self, atom):
        r2 = 0
        for x1, x2 in zip(self.coord, atom.coord):
            r2 += (x1-x2)*(x1-x2)
            
        return math.sqrt(r2)

class Protonate:
    """ Protonates atoms using VSEPR theory """

    def __init__(self, verbose=False):
        self.verbose=verbose

        self.valence_electrons = {
                                  'N': 5,
                                  'O': 6,
                                }

        self.standard_charges= {'ARG-NH1':1.0,
                                'ASP-OD2':-1.0,
                                'GLU-OE2':-1.0,
                                'HIS-ND1':1.0,
                                'LYS-NZ':1.0,
                                'ARG-NE' :1,
                                'N+':1.0,
                                'C-':-1.0}


        self.sybyl_charges = {'N.pl3':+1,
                              'N.3':+1,
                              'N.4':+1,
                              'N.ar':+1,
                              'O.co2-':-1}


        self.bond_lengths = {'C':1.09,
                             'N':1.01,
                             'O':0.96,
                             'F':0.92,
                             'Cl':1.27,
                             'Br':1.41,
                             'I':1.61,
                             'S':1.35}
        
        self.number_of_pi_electrons_in_bonds_in_backbone = { 'O':1}

        self.number_of_pi_electrons_in_conjugate_bonds_in_backbone = {'N':1}

        self.number_of_pi_electrons_in_bonds_in_sidechains = {'ARG-NH1':1,
                                                              'ASN-OD1':1,
                                                              'ASP-OD1':1,
                                                              'GLU-OE1':1,
                                                              'GLN-OE1':1,
                                                              'HIS-ND1':1}
        
        self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains = {
                                                                'ARG-NH2':1,
                                                                'ASN-ND2':1,
                                                                'GLN-NE2':1,
                                                                'HIS-NE2':1,
                                                                'TRP-NE1':1,
                                                                'GLU-OE2':1,
                                                                'ASP-OD2':1,
                                                                'ARG-NE' :1}


        self.number_of_pi_electrons_in_bonds_ligands = {'N.pl3':0,
                                                        'O.2':1,
                                                        'O.co2':1,
                                                        'N.ar':1,
                                                        'N.1':2}

        self.number_of_pi_electrons_in_conjugate_bonds_in_ligands = {'N.am':1,'N.pl3':1}

        self.protonation_methods = {4:self.tetrahedral,
                                    3:self.trigonal}

        self.moleculeGraph = {}
        
        self.anionId = -1
        self.atomList = []
        self.hydrogenAtomsList = []
        self.connected2H = []
        
        return




    def protonate(self, atomList, anionAtom):
        self.atomList = atomList
        self.moleculeGraph, self.anionId = molecule2graph(atomList, anionAtom, False)
        self.moleculeGraph = self.moleculeGraph.copy()
        
        # protonate all atoms
        for atomId in self.moleculeGraph.nodes():
            atom = self.atomList[atomId]
            element = atom.element
            
            if element in ["O", "N" ] and atom != anionAtom:
#                print("protonuje: ", atom.get_parent().get_resname(), atom.get_name())
                self.protonate_atom(atomId)
                
        hInd = len(self.atomList)
        for connInd in self.connected2H:
            self.moleculeGraph.add_edge(hInd, connInd)
            hInd +=1
        self.atomList += self.hydrogenAtomsList

        return


    def set_charge(self, atomId):
        atom = self.atomList[atomId]
        # atom is a protein atom

        key = '%3s-%s'%(atom.get_parent().get_resname(), atom.get_name())
        self.moleculeGraph.nodes[atomId]["key"] = key
        if key in list(self.standard_charges.keys()):
            self.moleculeGraph.nodes[atomId]["charge"] = self.standard_charges[key]
        else:
            self.moleculeGraph.nodes[atomId]["charge"] = 0

#        print("ustalono ladunek: ", self.moleculeGraph.nodes[atomId]["charge"])
        return 

    def protonate_atom(self, atom):

        self.set_charge(atom)
        self.set_number_of_pi_electrons(atom)
        self.set_number_of_protons_to_add(atom)
        self.set_steric_number_and_lone_pairs(atom)
        self.add_protons(atom)
        return


    def set_number_of_pi_electrons(self, atom):
        atomKey = self.moleculeGraph.nodes[atom]["key"]
        
        atomObj = self.atomList[atom]
        aminoacid = False
        if atomObj.get_parent().get_resname() in [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "ILE", "LEU", "LYS", "MET", "PRO", "SER", "THR", "VAL", "HIS", "TRP", "PHE" , "TYR" ] :
                aminoacid = True
        
        if aminoacid:
            if atomObj.get_name() in [ "O", "OXT"] :
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"] = 1
            elif atomKey in self.number_of_pi_electrons_in_bonds_in_sidechains:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"] = self.number_of_pi_electrons_in_bonds_in_sidechains[atomKey]
            else:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"] = 0
                
            if atomObj.get_name() == "N" :
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_conjugate_bonds"] = 1
            elif atomKey in self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_conjugate_bonds"] = self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains[atomKey]
            else:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_conjugate_bonds"] = 0
                
        else:
            
            if atomObj.get_name() in self.number_of_pi_electrons_in_bonds_ligands:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"] = self.number_of_pi_electrons_in_bonds_ligands[atomObj.get_name()]
            else:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"] = 0
                
            if atomObj.get_name() in self.number_of_pi_electrons_in_conjugate_bonds_in_ligands:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_conjugate_bonds"] = self.number_of_pi_electrons_in_conjugate_bonds_in_ligands[atomObj.get_name()]
            else:
                self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_conjugate_bonds"] = 0
                
#        print("elektrony pi w wiazaniach 2 i 3: ",self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"]  )
#        print("elektrony pi w sprzezonych wiazaniach ", self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_conjugate_bonds"])

    def set_number_of_protons_to_add(self, atom):
        number_of_protons_to_add  = 8
        number_of_protons_to_add -= self.valence_electrons[self.atomList[atom].element]
        number_of_protons_to_add -= len(list(nx.neighbors(self.moleculeGraph, atom)))
        
        number_of_pi_electrons_in_double_and_triple_bonds = self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"]
        number_of_protons_to_add -= number_of_pi_electrons_in_double_and_triple_bonds
        number_of_protons_to_add += int(self.moleculeGraph.nodes[atom]["charge"])
        
        self.moleculeGraph.nodes[atom]["number_of_protons_to_add"] =  number_of_protons_to_add
#        print("liczba protonow do dodania: ",self.moleculeGraph.nodes[atom]["number_of_protons_to_add"]  )

    def set_steric_number_and_lone_pairs(self, atom):
        steric_number = 0

        steric_number += self.valence_electrons[self.atomList[atom].element]
        steric_number += len(list(nx.neighbors(self.moleculeGraph, atom)))
        steric_number += self.moleculeGraph.nodes[atom]["number_of_protons_to_add"]
        steric_number -= self.moleculeGraph.nodes[atom]["charge"]
        steric_number -= self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_bonds"]
        steric_number -= self.moleculeGraph.nodes[atom]["number_of_pi_electrons_in_conjugate_bonds"]

        self.moleculeGraph.nodes[atom]["steric_number"] = math.floor(steric_number/2.0)

        self.moleculeGraph.nodes[atom]["number_of_lone_pairs"] = steric_number - len(list(nx.neighbors(self.moleculeGraph, atom))) - self.moleculeGraph.nodes[atom]["number_of_protons_to_add"]

        self.moleculeGraph.nodes[atom]["steric_number_and_lone_pairs_set"] = True
#        print("liczba steryczna: ", self.moleculeGraph.nodes[atom]["steric_number"] )
        return


    def add_protons(self, atom):
        # decide which method to use
        
        if self.moleculeGraph.nodes[atom]["steric_number"] in list(self.protonation_methods.keys()):
            self.protonation_methods[self.moleculeGraph.nodes[atom]["steric_number"]](atom)

        return


    def trigonal(self, atom):
        number_of_protons_to_add = self.moleculeGraph.nodes[atom]["number_of_protons_to_add"]
        
        if number_of_protons_to_add == 0:
            return
        
        rot_angle = math.radians(120.0)
        bonded_atoms_ids = list(self.moleculeGraph.neighbors(atom))
        
        if len(bonded_atoms_ids) == 0:
            return

        if len(bonded_atoms_ids) == 1:
            A = self.atomList[atom].get_coord()
            B = self.atomList[ bonded_atoms_ids[0] ].get_coord()
            
            Bneighbors = list(self.moleculeGraph.neighbors(bonded_atoms_ids[0]))
            for cCandidate in Bneighbors:
                if cCandidate != atom and self.atomList[cCandidate].element in [ "N" , "C" ]:
                    C = self.atomList[cCandidate].get_coord()
                    norm_vec = -normalize(np.cross(A-B, B-C))
                    break
            else:
                norm_vec = get_ortonormal(B-A)
            
            bondDirection = B-A
            for i in range(number_of_protons_to_add):
                bondDirection = rotateVector(bondDirection, norm_vec, rot_angle)
                bondDirection = normalize(bondDirection)* self.bond_lengths[ self.atomList[atom].element ]
                
                newAtomCoords = A + bondDirection
                self.hydrogenAtomsList.append(HydrogenAtom(newAtomCoords))
                self.connected2H.append(atom)
                
        elif len(bonded_atoms_ids) == 2:
            A = self.atomList[atom].get_coord()
            B = self.atomList[ bonded_atoms_ids[0] ].get_coord()
            C = self.atomList[ bonded_atoms_ids[1] ].get_coord()
            
            AB = normalize(B-A)
            AC = normalize(C-A)
            
            bondDirection = -(AB+AC)
            bondDirection = normalize(bondDirection)* self.bond_lengths[ self.atomList[atom].element ]
            
            newAtomCoords = A + bondDirection
            
            self.hydrogenAtomsList.append(HydrogenAtom(newAtomCoords))
            self.connected2H.append(atom)

        return


    def tetrahedral(self, atom):
        number_of_protons_to_add = self.moleculeGraph.nodes[atom]["number_of_protons_to_add"]
        rot_angle = math.radians(109.5)
        
        if number_of_protons_to_add == 0:
            return
        
        bonded_atoms_ids = list(self.moleculeGraph.neighbors(atom))
        
        if len(bonded_atoms_ids) == 0:
            return

        if len(bonded_atoms_ids) == 1:
            A = self.atomList[atom].get_coord()
            B = self.atomList[ bonded_atoms_ids[0] ].get_coord()
            
            norm_vec = get_ortonormal(B-A)
            dih_rot = math.radians(120)
            bondDirection = rotateVector(B-A, norm_vec, rot_angle)
            
            for i in range(number_of_protons_to_add):
                bondDirection = rotateVector(bondDirection, B-A, dih_rot)
                bondDirection = normalize(bondDirection)* self.bond_lengths[ self.atomList[atom].element ]
                
                newAtomCoords = A + bondDirection
                self.hydrogenAtomsList.append(HydrogenAtom(newAtomCoords))
                self.connected2H.append(atom)
        # 1 bond
        
        elif len(bonded_atoms_ids) == 2:
            A = self.atomList[atom].get_coord()
            B = self.atomList[ bonded_atoms_ids[0] ].get_coord()
            C = self.atomList[ bonded_atoms_ids[1] ].get_coord()
            
            AB = normalize(B-A)
            AC = normalize(C-A)
            
            axis = AB+AC
            bondDirection = rotateVector(-AB, axis, math.radians(90))
            bondDirection = normalize(bondDirection)*self.bond_lengths[ self.atomList[atom].element ]
            
            newAtomCoords = A + bondDirection
            self.hydrogenAtomsList.append(HydrogenAtom(newAtomCoords))
            self.connected2H.append(atom)
            
            if number_of_protons_to_add > 1:
                bondDirection = rotateVector(bondDirection, axis, math.radians(180))
                newAtomCoords = A + bondDirection
                self.hydrogenAtomsList.append(HydrogenAtom(newAtomCoords))
                self.connected2H.append(atom)
                
        elif len(bonded_atoms_ids) == 3:
            A = self.atomList[atom].get_coord()
            B = self.atomList[ bonded_atoms_ids[0] ].get_coord()
            C = self.atomList[ bonded_atoms_ids[1] ].get_coord()
            D = self.atomList[ bonded_atoms_ids[2] ].get_coord()
            
            AB = normalize(B-A)
            AC = normalize(C-A)
            AD = normalize(D-A)
            
            bondDirection = -(AB+AC+AD)
            bondDirection = bondDirection*self.bond_lengths[ self.atomList[atom].element ]
            
            
            newAtomCoords = A + bondDirection
            self.hydrogenAtomsList.append(HydrogenAtom(newAtomCoords))
            self.connected2H.append(atom)
            

        return

    
def get_ortonormal(vec):
    closest2zeroIndex = -1
    dist = 100
    for i, el in enumerate(vec):
        if abs(el) < dist:
            closest2zeroIndex = i
            dist = abs(el)
            
    vecOut = np.array([0.0, 0.0, 0.0])
    
    aInd = (closest2zeroIndex+1)%3
    bInd = (closest2zeroIndex+2)%3
    
    if abs(vec[bInd]) < 0.001:
        vecOut[bInd]=1
        return vecOut
    
    vecOut[aInd] =1
    vecOut[bInd] = -vec[aInd]/vec[bInd]
    
    return normalize(vecOut)