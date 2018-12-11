#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:37:23 2018

@author: michal
"""

def noAAinPiAcids(actualData):
    return actualData[(actualData["Pi acid Code"] != "TYR") & (actualData[ "Pi acid Code" ] != "PHE" ) & ( actualData["Pi acid Code"] != "HIS" ) & ( actualData["Pi acid Code"] != "TRP" ) ]

def noAAinPiRes(actualData):
    return actualData[(actualData["Pi res code"] != "TYR") & (actualData[ "Pi res code" ] != "PHE" ) & ( actualData["Pi res code"] != "HIS" ) & ( actualData["Pi res code"] != "TRP" ) ]

def noAAinAnions(actualData):
    return actualData[(actualData["Anion code"] != "GLU") & (actualData[ "Anion code" ] != "ASP" ) & ( actualData["Anion code"] != "TYR" ) & ( actualData["Anion code"] != "CYS") ]
   
def noNUinAnions(actualData):
    return actualData[(actualData["Anion code"] != "G") & (actualData[ "Anion code" ] != "A" ) & ( actualData["Anion code"] != "C" ) & ( actualData["Anion code"] != "T") & ( actualData["Anion code"] != "U") ]
   
def noNUinPiAcids(actualData):
    return actualData[(actualData["Pi acid Code"] != "G") & (actualData[ "Pi acid Code" ] != "A" ) & ( actualData["Pi acid Code"] != "C" ) & ( actualData["Pi acid Code"] != "T" ) & ( actualData["Pi acid Code"] != "U" )  ]

def noNUinPiRes(actualData):
    return actualData[(actualData["Pi res code"] != "G") & (actualData[ "Pi res code" ] != "A" ) & ( actualData["Pi res code"] != "C" ) & ( actualData["Pi res code"] != "T" ) & ( actualData["Pi res code"] != "U" ) ]

def noAAinHAcceptors(actualData):
    return actualData[(actualData["Acceptor code"] != "GLU") & (actualData[ "Acceptor code" ] != "ASP" ) & ( actualData["Acceptor code"] != "TYR" ) & ( actualData["Acceptor code"] != "CYS") ]

def noNUinHAcceptors(actualData):
    return actualData[(actualData["Acceptor code"] != "G") & (actualData[ "Acceptor code" ] != "A" ) & ( actualData["Acceptor code"] != "C" ) & ( actualData["Acceptor code"] != "T" ) & ( actualData["Acceptor code"] != "U" ) ]

def noAAinHDonors(actualData):
    return actualData[(actualData["Donor code"] != "TYR") & (actualData[ "Donor code" ] != "PHE" ) & ( actualData["Donor code"] != "HIS" ) & ( actualData["Donor code"] != "TRP" ) ]

def noNUinHDonors(actualData):
    return actualData[(actualData["Donor code"] != "G") & (actualData[ "Donor code" ] != "A" ) & ( actualData["Donor code"] != "C" ) & ( actualData["Donor code"] != "T" ) & ( actualData["Donor code"] != "U" ) ]
