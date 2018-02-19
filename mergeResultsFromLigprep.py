# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 19:24:23 2018

@author: michal
"""

from PDB_requests import getLigandCodeFromSdf
import pandas as pd

sdfFromLigprep = "sdf/ligprep_2-out_cutted.sdf"
fullLigandAnionData = "logs/AnionPiLigandFullxhThresholds.log"

anionsCodes = getLigandCodeFromSdf(sdfFromLigprep)
anionsCodes = list(set(anionsCodes))
anionsCodes.append("CL")

logData = pd.read_table(fullLigandAnionData)
logData = logData[ logData["Residue Name"].isin(anionsCodes) ]
    
logData = logData[ [ "Ligand Code" , "PDB Code" ] ]
logData = logData.drop_duplicates()

tolerances = pd.read_table("resolu.idx" , sep = "\t;\t", engine = "python")
tolerances = tolerances.drop_duplicates(subset = "PDB Code") 
#Zduplikowane PDB nie roznia sie istotnie rozdzielczoscia
#Tylko do odsiewu

logData = logData.merge(tolerances, on = "PDB Code")


logData = logData[  logData["RESOLUTION"] < 2.5  ]

logData = logData.sort_values( by = [ "PDB Code" ] )
logData.to_csv( "logs/MergeResultsFromLigprep.log", sep = ":", columns = [ "Ligand Code" , "PDB Code" ] , header = False, index = False )