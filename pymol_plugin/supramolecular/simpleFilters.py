#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:37:23 2018

@author: michal
"""
cationicAA = ["LYS","ARG"]
acidicAA = ["ASP","GLU"]
restAA = ["ALA", "CYS","GLY","ILE","LEU","MET","ASN","PRO","GLN","SER","THR","VAL"]
ringAA = ["PHE", "HIS", "TRP", "TYR"]
NU = ["A","G","T","C","U","I","DA", "DC", "DG", "DT", "DI" ]

def noAAinPiAcids(actualData):
    return actualData[ ~actualData["Pi acid Code"].isin(ringAA ) ]

def noAAinPiRes(actualData):
    return actualData[ ~actualData["Pi res code"].isin(ringAA) ]

def noAAinAnions(actualData):
    return actualData[~actualData["Anion code"].isin( acidicAA) ]
   
def noNUinAnions(actualData):
    return actualData[ ~actualData["Anion code"].isin( NU) ]
   
def noNUinPiAcids(actualData):
    return actualData[~actualData["Pi acid Code"].isin(NU )  ]

def noNUinPiRes(actualData):
    return actualData[~actualData["Pi res code"].isin( NU ) ]

def noAAinHAcceptors(actualData):
    return actualData[~actualData["Acceptor code"].isin( acidicAA) ]

def noNUinHAcceptors(actualData):
    return actualData[~actualData["Acceptor code"].isin( NU ) ]

def noAAinHDonors(actualData):
    return actualData[~actualData["Donor code"].isin([ "TYR" , "PHE" , "HIS" ,"TRP" ]) ]

def noNUinHDonors(actualData):
    return actualData[~actualData["Donor code"].isin( NU ) ]

def noAAinCations(actualData):
    return actualData[~actualData["Cation code"].isin(cationicAA)]