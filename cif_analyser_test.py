# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 18:24:14 2017

@author: michal
"""
from cif_analyser import getRingsCentroids
import numpy as np
from fakeBiopython import  myResidue
        
def writeLigandAndCentorids( ligandAtoms, centroids ):
    atomNo = len(ligandAtoms) + len(centroids)*2
    
    xyz = open("xyz/cit_centroids.xyz", 'a+')
    xyz.write( str(atomNo)+"\n\n" )
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
#myRes.read_xyz("testosterone.xyz")
#centroids = getRingsCentroids( myRes )
#myRes.read_xyz("ligand.xyz")
#centroids = getRingsCentroids( myRes )
#writeLigandAndCentorids( myRes.get_atoms(), centroids )
#myRes.read_xyz("peptyd.xyz", False)
#centroids = getRingsCentroids( myRes )
#writeLigandAndCentorids( myRes.get_atoms(), centroids )
myRes.read_xyz("xyz/ligands/UMP.xyz", False)
centroids = getRingsCentroids( myRes )
#writeLigandAndCentorids( myRes.get_atoms(), centroids )
#myRes.read_xyz("antracen.xyz", False)
#centroids = getRingsCentroids( myRes )
#writeLigandAndCentorids( myRes.get_atoms(), centroids )
#myRes.read_xyz("cycle_test.xyz", False)
#centroids = getRingsCentroids( myRes )
#writeLigandAndCentorids( myRes.get_atoms(), centroids )
#myRes.read_xyz("C60.xyz", False)
#centroids = getRingsCentroids( myRes )
writeLigandAndCentorids( myRes.get_atoms(), centroids )