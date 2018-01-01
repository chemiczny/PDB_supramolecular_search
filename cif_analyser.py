# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 13:56:41 2017

Skrypt powstał we wspolpracy z najpiekniejsza kobieta na swiecie, Zrodzona
w Olkuszu, Pania mojego serca Emilia Kuzniak.

Sluzy do analizy pojedynczego pliku cif celem znalezienia potencjalnych struktur
zawierajacych oddzialywania supramolekularne anion-pi oraz obliczenie
wielkosci geometrycznych tychze czasteczek na potrzeby dalszej analizy
"""


from Bio.PDB import MMCIFParser, NeighborSearch, Selection
from math import sqrt, acos, degrees
import time
import networkx as nx
import numpy as np

def normalize(v):
    """
    Funkcja pomocnicza. Normalizuje wektor v.
    
    Wejscie:
    v - numpy numeric array, wektor do normalizacji
    
    Wyjcie:
    v - znormalizowany wektor v
    """
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def getAverageCoords ( allAtomsList, atomsIndList ):
    """
    Oblicz srednie wspolrzedne x, y, z wybranych atomow
    
    Wejscie:
    allAtomsList - lista obiektow Atom (cala czasteczka)
    atomsIndList - lista indeksow atomow, ktorych usrednione polozenie
                    ma byc obliczone
                    
    Wyjscie:
    averageCoords - 3-elementowa lista zawierajaca usrednione wspolrzedne 
                    x, y, z wybranych atomow
    """
    averageCoords = [ 0. , 0. , 0. ]
    
    for atomInd in atomsIndList:
        new_coords = allAtomsList[atomInd].get_coord()
        
        for coord in range(3):
            averageCoords[coord] += new_coords[coord]
            
    atomsNo = float(len(atomsIndList))
    
    for coord in range(3):
        averageCoords[coord] /= atomsNo
        
    return averageCoords
    
def isFlat(allAtomsList, atomsIndList):
    """
    Sprawdz czy wybrane atomy leza w jednej plaszczyznie.
    Procedura: na podstawie polozen trzech pierwszych atomow wyznacza sie 
    wektor normalnych do wyznaczonej przez nich plaszczyzny. Nastepnie 
    sprawdzane jest czy kolejne wiazania tworza wektory prostopadle
    do wektora normalnego. Dopuszczalne jest odchylenie 5 stopni.
    
    TODO: Nie lepiej byloby obliczyc tensor momentu bezwladnosci, 
    zdiagonalizowac go i rozstrzygnac na podstawie jego wartosci wlasnych?
    
    Wejscie:
    allAtomsList - lista obiektow Atom (cala czasteczka)
    atomsIndList - lista indeksow atomow, ktore maja byc zweryfikowane
                    pod wzgledem lezenia w jednej plaszczyznie
                    
    Wyjscie:
    verdict - slownik, posiada klucze: isFlat (zmienna logiczna, True jesli
                struktura jest plaska), normVec (3-elementowa lista float,
                wspolrzedne wektora normalnego plaszczyzny, jesli struktura nie
                jest plaska wszystkie jego wspolrzedne sa rowne 0)
    """
        
    verdict = { 'isFlat' : False, 'normVec' : [ 0, 0, 0] }
    
    if len(atomsIndList) < 3:
        print("Znaleziono 2-wu elementowy cykl!!!")
        return verdict
    
    A = np.array( allAtomsList[ atomsIndList[0] ].get_coord() )
    B = np.array( allAtomsList[ atomsIndList[1] ].get_coord() )
    C = np.array( allAtomsList[ atomsIndList[2] ].get_coord() )
    
    vec1 = A-B
    vec2 = B-C
    
    norm_vec = normalize(np.cross(vec1, vec2))
    
    if len(atomsIndList) <= 3:
#        print("LOL strasznie krotki pierscien")
        #verdict['isFlat'] = True
        verdict['normVec'] = norm_vec
        return verdict    
    
    lastAtom = C    
    for i in range( 3, len(atomsIndList)  ):
        atomInd = atomsIndList[i]
        D = np.array( allAtomsList[ atomInd ].get_coord() )
        new_vec = normalize( lastAtom - D )
        
        if abs( np.inner( new_vec, norm_vec ) ) > 0.09:
            return verdict
            
    verdict['isFlat'] = True
    verdict['normVec'] = norm_vec
    return verdict

def getRingsCentroids( molecule ):
    """
    Znajdz pierscienie w czasteczce i wyznacz wspolrzedne ich srodkow jesli
    sa one aromatyczne.
    
    Procedura: 
    - utworz graf na podstawie danych o atomach, wszystkie atomy lezace blizej 
        niz 1.8 A sa traktowane jako wierzcholki grafu polaczone krawedzia
    - znajdz fundamentalne cykle w grafie (Stosowany algorytm:
    Paton, K. An algorithm for finding a fundamental set of cycles of a graph. 
    Comm. ACM 12, 9 (Sept 1969), 514-518.)
    - odrzuc cykle, ktore skladaja sie z wiecej niz 6 wierzcholkow
    - sprawdz czy znalezione pierscienie sa plaskie
    
    Wejscie:
    -molecule - obiekt Residue (Biopython)
    
    Wyjscie:
    -centroids - lista slownikow z danymi o znalezionych (potencjalnie) 
        aromatycznych pierscienieniach. Slownik zawiera klucze:
        coords (wspolrzedne srodka pierscienia), 
        normVec (wektor normalnych plaszczyzny pierscienia)
    """
    atoms = list(molecule.get_atoms())
    G = nx.Graph()
    
    for atom1Ind in range(len(atoms)):
        atom1 = atoms[atom1Ind]
        if atom1.get_name() == "H":
            continue
        
        for atom2Ind in range(atom1Ind+1, len(atoms)):
            atom2 = atoms[atom2Ind]            
            if atom2.get_name() == "H":
                continue
            
            distance = atom1 - atom2
            
            if distance < 1.8 :
                G.add_edge(atom1Ind, atom2Ind)
                
    cycles = list(nx.cycle_basis(G))
    centroids = []
    for cycle in cycles:
        if len(cycle) > 6:
            continue
        
        flatAnalyse = isFlat(atoms, cycle)
        if not flatAnalyse['isFlat']:
            continue
        
        centroids.append({ "coords" : getAverageCoords( atoms, cycle), "normVec" : flatAnalyse["normVec"] })
        
    return centroids
        

def findSupramolecularAnionPiLigand( ligandCode, cifFile, PDBcode ):
    """
    Przeanalizuj pojedynczy plik cif pod katem oddzialywan suporamolekularnych
    zwiazanych z konkretnym ligandem
    
    Wejscie:
    ligandCode - kod liganda
    cifFile    - plik cif pobrany z bazy PDB
    
    Wyjscie:
    Hehehe, czas pokaze...
    """
    parser = MMCIFParser()
    
    try:
        structure = parser.get_structure('temp', cifFile)
    except:
        print("Biopytong nie ogarnia!", cifFile)
        #Zeby zobaczyc co sie dzieje
        return True        
        
    atoms = Selection.unfold_entities(structure, 'A')  
    ns = NeighborSearch(atoms)
    supramolecularFound = False    
    
    ligands = []
    for residue in structure.get_residues():
        if ligandCode == residue.get_resname():
            ligands.append(residue)
    for ligand in ligands:
        centroids = getRingsCentroids( ligand )
        #print(centroids)
        
        for centroid in centroids:
            neighbors = ns.search(np.array(centroid["coords"]), 5, 'A')
            extractedAtoms = extractNeighbours( neighbors, ligandCode )
            if len(extractedAtoms) > 0:
                supramolecularFound = True
            
            writeSupramolecularSearchResults(ligandCode, PDBcode, centroid, extractedAtoms)
            
    return supramolecularFound
    
            
def writeSupramolecularSearchHeader( ):
    """
    Zapisz naglowki do pliku z wynikami:
    """
    resultsFileName = "logs/anionPiLigand.log"
    resultsFile = open(resultsFileName, "w")
    resultsFile.write("PDB Code\tLigand Code\tResidue Name\tAnion type\t")
    resultsFile.write("Atom symbol\tDistance\tAngle\n")
    resultsFile.close()
    
            
def writeSupramolecularSearchResults( ligandCode, PDBcode, centroid, extractedAtoms ):
    """
    Zapisz dane do pliku z wynikami
    """
    resultsFileName = "logs/anionPiLigand.log"
    
    resultsFile = open(resultsFileName, "a+")
    
    for atomData in extractedAtoms:
        residueName = atomData["Atom"].get_parent().get_resname()
        resultsFile.write(PDBcode+"\t")
        resultsFile.write(ligandCode+"\t")
        resultsFile.write(residueName+"\t")
        resultsFile.write(atomData["AnionType"]+"\t")
        resultsFile.write(atomData["Atom"].element+"\t")
        
        distance = atomDistanceFromCentroid( atomData["Atom"], centroid )
        angle = atomAngleNomVecCentroid( atomData["Atom"], centroid )
        
        resultsFile.write(str(distance)+"\t")
        resultsFile.write(str(angle)+"\n")
        
        
    
    resultsFile.close()
    
def atomDistanceFromCentroid( atom, centroid ):
    """
    Funkcja pomocnicza, oblicza odleglosc pomiedzy atomem a srodkiem 
    pierscienia
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    centroid - slownik, klucze: coords, normVec
    
    Wyjscie:
    odleglosc (float)
    """
    atomCoords = atom.get_coord()
    centoridCoords = centroid["coords"]
    
    dist = 0
    for atomCoord, centroidCoord in  zip( atomCoords, centoridCoords ):
        dist+= (atomCoord-centroidCoord)*(atomCoord-centroidCoord)
        
    return sqrt(dist)

def atomAngleNomVecCentroid( atom, centroid ):
    """
    Funkcja pomocnicza, oblicza kat pomiedzy kierunkiem od srodka 
    pierscienia do atomu a wektorem normalnym plaszczyzny pierscienia
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    centroid - slownik, klucze: coords, normVec
    
    Wyjscie:
    kat w stopniach
    """
    atomCoords = np.array(atom.get_coord())
    centroidCoords = np.array(centroid["coords"])
    normVec = centroid["normVec"]
    
    centrAtomVec= normalize( atomCoords - centroidCoords )
    inner_prod = np.inner( normVec, centrAtomVec )
    
    return degrees( acos(inner_prod) )
    
            
def extractNeighbours( atomList, ligandCode ):
    """
    Wydziel atomy, ktore moga byc anionami w sasiedztwie liganda
    
    Wejscie:
    atomList - lista obiektow Atom (Biopython)
    ligandCode - kod liganda
    
    Wysjcie:
    extractedAtoms - lista slownikow z informacjami o potencjalnym anionie,
        klucze: Atom - atom, AnionType - rodzaj anionu
    """
    extractedAtoms = []
    
    for atom in atomList:
        element_symbol = atom.element
        
        if not isItWorthAnalyzing(atom, ligandCode):
            continue
        
        potentiallySupramolecular = False
        anionType = ""
        
        if element_symbol == "O":
            potentiallySupramolecular, anionType = handleOxygen(atom)
        elif element_symbol in [ "F", "CL", "BR", "I" ]:
            potentiallySupramolecular, anionType = handleHalogens(atom)
        elif element_symbol in [ "S", "SE" ]:
            potentiallySupramolecular, anionType = handleChalcogens(atom)
        elif element_symbol in [ "N", "P" ]:
            potentiallySupramolecular, anionType = handlePnictogens(atom)
        else:
            potentiallySupramolecular, anionType = handleOthers(atom)            
            
        if potentiallySupramolecular:
#            print("cos watergo uwagi :D", atom.get_fullname())
            extractedAtoms.append({ "Atom" : atom, "AnionType" : anionType})
            
    return extractedAtoms

def isItWorthAnalyzing(atom, ligandCode):
    """
    Szybka, wstepna weyfikacja atomu jako potencjalnego anionu
    Odrzucane sa te atomy, ktore naleza do niekwasnych aminokwasow,
    czasteczki wody oraz ligandu
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    ligandCode - kod liganda
    
    Wyjscie:
    True/False - czy atom moze byc potencjalnym anionem?
    """
    parent_name = atom.get_parent().get_resname()
    
#    bad_parents = [ "HOH", "ALA", "ARG", "ASN", "CYS", "GLN", "GLY", "HIS",
#        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
        
    bad_parents = [ "HOH" ]
        
    if parent_name == ligandCode:
        return False
    
    if parent_name in bad_parents:
        return False
        
    return True

def handleOxygen( atom ):
    """
    Zweryfikuj atom tlenu jako potencjalny anion.
    W tej chwili jako aniony identyfikujemy wszystkie atomy tlenu, ktore:
    -nie naleza do czasteczki wody
    -jest sa zwiazane z atomem C, to atom C jest zwiazany z co najmniej dwoma
        atomami O (wykluczamy alkohole, etery, aldehydy, ketony...)
        
    TODO: a estry? a C-koniec bialka?
    
    Wejscie:
    atom - obiekt Atom (Biopython), tlen
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    graph, oxygenInd = molecule2graph(atom, atoms )
    
    oxygen_neighbors = []
    oxygen_neighbors = graph.neighbors(oxygenInd)
    
    if len(oxygen_neighbors) == 0:
        return True, "Complex-O?"
    elif len(oxygen_neighbors) > 1:
        print("Prawdopodobnie ester lub eter")
        return False, "O"
        
    oxygen_neighbor_index = oxygen_neighbors[0]
    oxygen_neighbor_symbol = atoms[ oxygen_neighbor_index ].element
    
    center_neighbors = graph.neighbors( oxygen_neighbor_index )
    
    oxygens_found = 0
    for unknown_atom in center_neighbors:
        if atoms[unknown_atom].element == "O":
            oxygens_found+=1
            
    if oxygen_neighbor_symbol == "C" and oxygens_found < 2:
        return False, oxygen_neighbor_symbol+str(oxygens_found)+"O"
    
    return True, oxygen_neighbor_symbol+"O"+str(oxygens_found)
    
    

def handleHalogens( atom):
    """
    Zweryfikuj atom fluorowca jako potencjalny anion.
    W chwili obecnej przyjmujemy, ze atom fluorowca jest anionem
    wtedy gdy jest nie zwiazany z zadna wieksza czasteczka (czyli wykrywamy
    aniony proste: F-, CL-, BR-, I-)
    
    Wejscie:
    atom - obiekt Atom (Biopython), fluorowiec
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) == 1:
        return True, atom.element
        
    graph, halogenInd = molecule2graph(atom, atoms )
    halogen_neighbors = graph.neighbors(halogenInd)
    
    if len(halogen_neighbors) == 0:
        return True, "Complex-"+atom.element+"?"
        
    return False, atom.element
    

def handleChalcogens(atom):
    """
    Zweryfikuj pozostale tlenowce (oprocz tlenu) jako potencjalne aniony.
    Podobnie jak w przypadku fluorowcow. W ten sposob ogarniamy tylko S2-
    
    Wejscie:
    atom - obiekt Atom (Biopython), tlenowiec
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) == 1:
        return True, atom.element
        
    graph, chalcogenInd = molecule2graph(atom, atoms )
    chalcogen_neighbors = graph.neighbors(chalcogenInd)
    
    if len(chalcogen_neighbors) == 0:
        return True, "Complex-"+atom.element+"?"
        
    return False, atom.element

def handlePnictogens(atom):
    """
    Zweryfikuj atom azotowca jako potencjalny anion.
    W chwili obecnej rozwazamy tylko jony CN- oraz N3-
    
    Wejscie:
    atom - obiekt Atom (Biopython), azotowiec
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    atoms = list(atom.get_parent().get_atoms())
    
    if len(atoms) == 2:
    
        is_carbon = False
        
        for unknown in atoms:
            if unknown.element == "C":
                is_carbon = True
                
        return is_carbon, "CN"
        
    elif len(atoms) ==3:
        onlyNitrogen = True
        
        for unknown in atoms:
            if unknown.element != "N":
                onlyNitrogen = False
                
        return onlyNitrogen, "N3"
    
    return False, atom.element

def handleOthers(atom):
    """
    Zweryfikuj niezakwalifikowany atom jako potencjalny anion.
    Obecnie nie wiemy jakie kryteria tu wrzucic wiec ta procedura zawsze zwraca
    False
    
    Wejscie:
    atom - obiekt Atom (Biopython)
    
    Wyjscie:
    True/False, anionType - czy atom moze byc potencjalnym anionem?, rodzaj
                            anionu (string)
    """
    return False, atom.element

def molecule2graph( atom, atoms ):
    """
    Konwersja czasteczki na graf (networkx)
    
    Wejscie:
    atom - obiekt Atom (Biopython), ktorego polozenie w grafie jest istotne
    atoms - lista wszystkich atomow
    
    Wyjscie:
    G, atomInd - graf (networkx), indeks wejsciowego atomu (wierzcholek w grafie)
    """
    atomName = atom.get_fullname()    
    thresholds = { "C" : 1.8, "O" : 1.8, "N" : 1.8, "S" : 2.1,
                  "F" : 1.6, "CL" : 2.0, "BR" : 2.1, "I" : 2.2 }
    
    G = nx.Graph()
    atoms_found = []
    for atom1Ind in range(len(atoms)):
        threshold1 = 2.2
        atom1 = atoms[atom1Ind]
        if atom1.element == "H":
            continue
        
        if atom1.element in thresholds.keys():
            threshold1 = thresholds[atom1.element]
        
        if atom1.get_fullname() == atomName:
            atoms_found.append(atom1Ind)
        
        for atom2Ind in range(atom1Ind+1, len(atoms)):
            threshold2 = 2.2
            atom2 = atoms[atom2Ind]            
            if atom2.element == "H":
                continue
            
            if atom2.element in thresholds.keys():
                threshold2 = thresholds[atom2.element]
            
            distance = atom1 - atom2
            
            threshold = max( threshold1, threshold2 )
            if distance < threshold :
                G.add_edge(atom1Ind, atom2Ind)
                
    if len(atoms_found) != 1:
        print("WTF!? ", atoms_found)
        
    if not atoms_found[0] in G.nodes():
        print("Nie znalazlem "+ atom.element+" w grafie!")
        G.add_node( atoms_found[0] )
        
    return G, atoms_found[0]

if __name__ == "__main__":
    writeSupramolecularSearchHeader( )
    timeStart = time.time()
    findSupramolecularAnionPiLigand( "HPA", "cif/3nrz.cif", "3NRZ" )
#    findSupramolecularAnionPiLigand( "LUM", "cif/1he5.cif", "1HE5" )
#    findSupramolecularAnionPiLigand( "NAP", "cif/3bcj.cif", "3BCJ" )
#    findSupramolecularAnionPiLigand( "NAP", "cif/4lbs.cif", "4LBS" )
    timeStop = time.time()
    print("Calosc: ", timeStop-timeStart)
#findSupramolecular( "LUM", 666 )