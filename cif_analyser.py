# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 13:56:41 2017

Skrypt powstał we wspolpracy z najpiekniejsza kobieta na swiecie, Zrodzona
w Olkuszu, Pania mojego serca Emilia Kuzniak.

Sluzy do analizy pojedynczego pliku cif celem znalezienia potencjalnych struktur
zawierajacych oddzialywania supramolekularne anion-pi oraz obliczenie
wielkosci geometrycznych tychze czasteczek na potrzeby dalszej analizy
"""


from Bio.PDB import FastMMCIFParser, NeighborSearch, Selection
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import time
import numpy as np
from supramolecularLogging import writeSupramolecularSearchHeader, writeSupramolecularSearchResults
from ringDetection import getRingsCentroids
from anionRecogniser import extractNeighbours
        

def findSupramolecularAnionPiLigand( ligandCode, cifFile, PDBcode, ligprepData = None ):
    """
    Przeanalizuj pojedynczy plik cif pod katem oddzialywan suporamolekularnych
    zwiazanych z konkretnym ligandem
    
    Wejscie:
    ligandCode - kod liganda
    cifFile    - plik cif pobrany z bazy PDB
    
    Wyjscie:
    Hehehe, czas pokaze...
    """
    parser = FastMMCIFParser()
    
    try:
        structure = parser.get_structure('temp', cifFile)
    except:
        print("Biopytong nie ogarnia!", cifFile)
        #Zeby zobaczyc co sie dzieje
        return True        
        
    supramolecularFound = False 
    
    for modelIndex, model in enumerate(structure):
        atoms = Selection.unfold_entities(model, 'A')  
        ns = NeighborSearch(atoms)
           
        for residue in model.get_residues():
            if ligandCode == residue.get_resname():
                supramolecularFound = supramolecularFound or analysePiacid(residue, PDBcode, modelIndex, ns, ligprepData)
            
    return supramolecularFound
    
def findSupramolecularAnionPiAllLigands( cifFile, PDBcode, ligprepData = None ):
    """
    Przeanalizuj pojedynczy plik cif pod katem oddzialywan suporamolekularnych
    zwiazanych z konkretnym ligandem
    
    Wejscie:
    ligandCode - kod liganda
    cifFile    - plik cif pobrany z bazy PDB
    
    Wyjscie:
    Hehehe, czas pokaze...
    """
    parser = FastMMCIFParser()
    
    try:
        structure = parser.get_structure('temp', cifFile)
    except:
        print("Biopytong nie ogarnia!", cifFile)
        #Zeby zobaczyc co sie dzieje
        return True        
        
    supramolecularFound = False 
    notPiacids = [ "HOH", "DOD", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "ILE", "LEU", "LYS", "MET", "PRO", "SER", "THR", "VAL" ] 
        
    resolution = readResolution(cifFile)
#    print("jade z pliku: ", cifFile)
    for modelIndex, model in enumerate(structure):
        atoms = Selection.unfold_entities(model, 'A')  
        ns = NeighborSearch(atoms)
           
        for residue in model.get_residues():
            residueName = residue.get_resname().upper()
            if not residueName in notPiacids  :
#                print("Analizuje: ", residueName)
                if analysePiacid(residue, PDBcode, modelIndex, ns, ligprepData, resolution):
                    supramolecularFound = True
            
    return supramolecularFound
    
def readResolution( cifFile ):
    mmcif_dict = MMCIF2Dict(cifFile)
    res_keys = ["_refine.ls_d_res_high" , "_reflns_shell.d_res_high" ]
    res_list = []
    for key in res_keys:
        if key in mmcif_dict:
            resolution = list(mmcif_dict[key])
            if len(resolution) > 1:
                return -2
            try:
                res_list.append(float(resolution[0]))
            except:
                print("jakis syf: ", cifFile, key, mmcif_dict[key])
                
    if not res_list:
        return -1
                
    variance = np.var(res_list)
    if variance < 0.001:
        return res_list[0]
    else:
        return -3
        
    
def analysePiacid(ligand, PDBcode, modelIndex, ns, ligprepData, resolution):
    centroids = getRingsCentroids( ligand )
    ligandCode = ligand.get_resname()
#    print("Znalazlem pierscienie w ilosci: ", len(centroids))
    
    ligandWithAnions = False
    allExtractedAtomsForLigand =[]
    
    for centroid in centroids:
        distance = 4.5
        neighbors = ns.search(np.array(centroid["coords"]), distance, 'A')
        extractedAtoms = extractNeighbours( neighbors, ligandCode, ns )
            
        if ligprepData:
            extractedAtoms = anionScreening( extractedAtoms, ligprepData )
            
        cationNear = []
        if len(extractedAtoms) > 0:
            cationNear = searchForCation( centroid["coords"], ns )
        
        extractedAtoms =  writeSupramolecularSearchResults(ligand, PDBcode, centroid, extractedAtoms, modelIndex, resolution, cationNear)
        
        if len(extractedAtoms) > 0:
            ligandWithAnions = True
            allExtractedAtomsForLigand += extractedAtoms
            saveLigandWithAnion(ligand, ligandCode, PDBcode, modelIndex, centroid, extractedAtoms)
        
    if ligandWithAnions:
        allAnions = getResiduesListFromAtomData( allExtractedAtomsForLigand ) 
        anionsNames = []
        
        for anion in allAnions:
            anionsNames.append( anion.get_resname() )
        
        saveLigand(ligand, ligandCode, PDBcode, modelIndex )
        saveLigandEnv(ligand, ligandCode, PDBcode, modelIndex, centroids, anionsNames, ns)
        
    return ligandWithAnions

def searchForCation ( point, ns  ):
    distance = 10
    neighbors = ns.search(np.array(point), distance, 'A')
    
    metalCationsFound = []
    metalCations =  [
        "Li", "Be", 
        "Na", "Mg", "Al", 
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", 
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", 
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", 
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]
    
    for atom in neighbors:
        if atom.element.upper() in [ element.upper() for element in metalCations ]:
            metalCationsFound.append( atom.element.upper() )
    
    return metalCationsFound
    
def getResiduesListFromAtomData( atomDataList ):
    atomsList = []
    for atomData in atomDataList:
        atomsList.append(atomData["Atom"])
    return Selection.unfold_entities(atomsList, 'R')  
    
def saveLigandWithAnion(ligand, ligandCode, PDBcode, modelIndex, centroid, extractedAtoms):        
    anions = getResiduesListFromAtomData( extractedAtoms )
    ligandAtoms = Selection.unfold_entities(ligand, 'A')  
    ligandId = str(ligand.get_id()[1])
    chain = ligand.get_parent().get_id()
    directory = "xyz/ligands_anions/"
    
    if ligandCode.upper() in [ "HIS" , "PHE", "TYR", "TRP" ]:
        directory = "xyz/aminoacids_anions/"
        
    for anion in anions:     
        atomsList = Selection.unfold_entities(anion, 'A') +  ligandAtoms
        
        anionName = anion.get_resname()
        anionId = str(anion.get_id()[1])
        anionChain = anion.get_parent().get_id()
        xyzName = directory+anionName+"_ANION.xyz"

        comment = "PDBCode: "+PDBcode+" Pi-acid-code: "+ligandCode+" Pi-acid-id: "+ligandId+" Pi-acid-chain: "+chain+" Anion-id: "+anionId+" Anion-chain: "+anionChain+"_ModelNo: "+str(modelIndex)
        appendXYZ( xyzName, atomsList, comment, [ centroid ] )
            
        

def saveLigand(ligand, ligandCode, PDBcode, modelIndex):
    ligandAtoms = Selection.unfold_entities(ligand, 'A')  
    ligandAtoms = Selection.unfold_entities(ligand, 'A')  
    ligandId = str(ligand.get_id()[1])
    chain = ligand.get_parent().get_id()
    
    directory = "xyz/ligands/"
    
    if ligandCode.upper() in [ "HIS" , "PHE", "TYR", "TRP" ]:
        directory = "xyz/aminoacids/"
        
    xyzName = directory+ligandCode+".xyz"
    comment = "PDBCode: "+PDBcode+" Pi-acid-id: "+ligandId+" Pi-acid-chain: "+chain+"_ModelNo: "+str(modelIndex)
    appendXYZ( xyzName, ligandAtoms, comment )

def saveLigandEnv(ligand, ligandCode, PDBcode, modelIndex, centroids, anionsNames, ns):
    ligandAtoms = Selection.unfold_entities(ligand, 'A')  
    neighbors = []
    distance = 4.5
    ligandId = str(ligand.get_id()[1])
    chain = ligand.get_parent().get_id()
    
    for atom in ligandAtoms:
        neighbors += ns.search(np.array(atom.get_coord()), distance, 'A')
        
    neighbors = Selection.unfold_entities(neighbors, 'R')
    atomsList = []
    for neighbor in neighbors:
        if not neighbor.get_resname() in [ "HOH", "DOD" ]:
            atomsList += Selection.unfold_entities(neighbor, 'A') 
         
    directory = "xyz/ligands_ENV/"
    if ligandCode.upper() in [ "HIS" , "PHE", "TYR", "TRP" ]:
        directory = "xyz/aminoacids_ENV/"
        
    xyzName = directory+ligandCode+"_ENV.xyz"
    comment = "PDBCode: "+PDBcode+" Pi-acid-id: "+ligandId+" Pi-acid-chain: "+chain+"_ModelNo: "+str(modelIndex)+"_Anions: "+",".join(anionsNames)
    
    appendXYZ( xyzName, atomsList, comment, centroids )
    
def appendXYZ( xyzName, atomsList, comment = "", centroids = [] ):
    atomNo = len(atomsList)+len(centroids)*2
    
    xyz = open(xyzName, 'a+')
    xyz.write( str(atomNo)+"\n" )
    xyz.write(comment+"\n")
    for atom in atomsList:
        coord = atom.get_coord()
        xyz.write(atom.element+" "+str(coord[0])+" "+str(coord[1])+" "+str(coord[2])+"\n")
        
    for centroid in centroids:
        normVec = centroid["normVec"]
        centroid = np.array(centroid["coords"])
        
        newGhost = centroid+normVec
        
        xyz.write("X "+str(centroid[0])+" "+str(centroid[1])+" "+str(centroid[2])+"\n")
        xyz.write("X "+str(newGhost[0])+" "+str(newGhost[1])+" "+str(newGhost[2])+"\n")
        
    xyz.close()
    
def anionScreening( atoms, ligprepData ):
    selectedAtoms = []
    anionsNames = ligprepData["anionNames"]
    
    for atomData in atoms:
        parentName = atomData["Atom"].get_parent().get_resname()        
        if parentName in anionsNames:
            selectedAtoms.append(atomData)
        
    return selectedAtoms
    

if __name__ == "__main__":
    writeSupramolecularSearchHeader( )
    timeStart = time.time()
#    findSupramolecularAnionPiLigand( "MCY", "cif/106d.cif", "106D" )
    #findSupramolecularAnionPiLigand( "NCZ", "cif/1j5i.cif", "1J5I" )
    findSupramolecularAnionPiAllLigands( "cif/1bp0.cif", "1BP0")
    findSupramolecularAnionPiAllLigands( "cif/3bdj.cif", "3BDJ")
#    findSupramolecularAnionPiLigand( "7NC", "cif/5wqk.cif", "5WQK" )
#    findSupramolecularAnionPiLigand( "HPA", "cif/3nrz.cif", "3NRZ" )
#    findSupramolecularAnionPiLigand( "LUM", "cif/1he5.cif", "1HE5" )
#    findSupramolecularAnionPiLigand( "NAP", "cif/3bcj.cif", "3BCJ" )
#    findSupramolecularAnionPiLigand( "NAP", "cif/4lbs.cif", "4LBS" )
    timeStop = time.time()
    print("Calosc: ", timeStop-timeStart)
#findSupramolecular( "LUM", 666 )
