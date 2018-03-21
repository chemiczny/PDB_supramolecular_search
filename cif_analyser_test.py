# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 18:24:14 2017

@author: michal
"""
from cif_analyser import getRingsCentroids
import numpy as np
from fakeBiopython import  myResidue
import glob
        
def writeLigandAndCentorids( ligandAtoms, centroids ):
    atomNo = len(ligandAtoms) + len(centroids)*2
    
    xyz = open("xyz/test_centroids.xyz", 'a+')
    xyz.write( str(atomNo)+"\ndupa\n" )
    for atom in ligandAtoms:
        coord = atom.get_coord()
        xyz.write(atom.get_name()[0]+" "+str(coord[0])+" "+str(coord[1])+" "+str(coord[2])+"\n")
        
    for centroid in centroids:
        normVec = centroid["normVec"]
        centroid = np.array(centroid["coords"])
        
        newGhost = centroid+normVec
        print(newGhost)
        
        xyz.write("X "+str(centroid[0])+" "+str(centroid[1])+" "+str(centroid[2])+"\n")
        xyz.write("X "+str(newGhost[0])+" "+str(newGhost[1])+" "+str(newGhost[2])+"\n")
        
    xyz.close()
    
    
myRes = myResidue()
ligands = glob.glob("xyz/xyz/ligands/*.xyz")
i = 0
for ligand in ligands:
    i+=1
    myRes.read_xyz(ligand)
    centroids = getRingsCentroids( myRes )
    writeLigandAndCentorids( myRes.get_atoms(), centroids )

print("przelecialo: ", i)