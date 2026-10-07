#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 12:53:24 2018

@author: michal

In older versions of Biopython method get_fullid() returns empty tuple
"""

def createResId(residue):
    idNo = residue.get_id()[1]
    chain = residue.get_parent().get_id()
    name = residue.get_resname()
    firstCoord = list(residue.get_atoms())[0].get_coord()
    coordStr = str(firstCoord[0])+str(firstCoord[1])+str(firstCoord[2])
    
    return chain+str(idNo)+name+coordStr

def createResIdFromAtom(atom):
    return createResId( atom.get_parent() )