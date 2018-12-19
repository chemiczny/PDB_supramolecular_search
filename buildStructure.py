#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 20:18:20 2018

@author: michal
"""

from Bio.PDB.StructureBuilder import StructureBuilder

def buildStructure(atomList, exceptResidue):
    sb = StructureBuilder() 
    sb.init_structure('pdb') 
    sb.init_seg(' ') 
    sb.init_model(0) 
    sb.init_chain('A') 
    
    # Atom 
#    parent_id = atomList[0].parent.parent.id 
#    sb.structure[0]['A'].id = parent_id 
    
    for resid, pdb_object in enumerate(atomList):
        if pdb_object.get_parent() != exceptResidue:
            sb.init_residue("WTF", " ", resid, " ")
            sb.structure[0]['A'].child_list[-1].add(pdb_object) 
   
              # Return structure 
    return sb.structure