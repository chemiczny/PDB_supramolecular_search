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
    if len(atomsIndList) <= 3:
        print("LOL strasznie krotki pierscien")
        return True
        
    verdict = { 'isFlat' : False, 'normVec' : [ 0, 0, 0] }
    A = np.array( allAtomsList[ atomsIndList[0] ].get_coord() )
    B = np.array( allAtomsList[ atomsIndList[1] ].get_coord() )
    C = np.array( allAtomsList[ atomsIndList[2] ].get_coord() )
    
    vec1 = A-B
    vec2 = B-C
    
    norm_vec = normalize(np.cross(vec1, vec2))
    
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
        

def findSupramolecularAnionPiLigand( ligandCode, cifFile ):
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
    structure = parser.get_structure('temp', cifFile)
    atoms = Selection.unfold_entities(structure, 'A')    
    
    ligands = []
    for residue in structure.get_residues():
        if ligandCode == residue.get_resname():
            print("Znalazlem!", ligandCode)
            ligands.append(residue)
    
    for ligand in ligands:
        centroids = getRingsCentroids( ligand )
        #print(centroids)
        
        ns = NeighborSearch(atoms)
        for centroid in centroids:
            neighbors = ns.search(np.array(centroid["coords"]), 5, 'A')
            extractedAtoms = extractNeighbours( neighbors, ligandCode )
            
def extractNeighbours( atomList, ligandCode ):
    """
    Wydziel atomy, ktore moga byc anionami w sasiedztwie liganda
    
    Wejscie:
    atomList - lista obiektow Atom (Biopython)
    ligandCode - kod liganda
    
    Wysjcie:
    jeszcze nieokreslone :)
    """
    extractedAtoms = []
    
    for atom in atomList:
        element_symbol = atom.element
        
        if not isItWorthAnalyzing(atom, ligandCode):
            continue
        
        potentiallySupramolecular = False
        
        if element_symbol == "O":
            potentiallySupramolecular = handleOxygen(atom)
        elif element_symbol in [ "F", "CL", "BR", "I" ]:
            potentiallySupramolecular = handleHalogens(atom)
        elif element_symbol in [ "S", "SE" ]:
            potentiallySupramolecular = handleChalcogens(atom)
        elif element_symbol in [ "N", "P" ]:
            potentiallySupramolecular = handlePnictogens(atom)
        else:
            potentiallySupramolecular = handleOthers(atom)            
            
        if potentiallySupramolecular:
            print("cos watergo uwagi :D", atom.get_fullname())
            extractedAtoms.append(atom)
            
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
    
    bad_parents = [ "HOH", "ALA", "ARG", "ASN", "CYS", "GLN", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]
        
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
    
    if len(oxygen_neighbors) > 1:
        print("Prawdopodobnie ester lub eter")
        return
        
    oxygen_neighbor_index = oxygen_neighbors[0]
    oxygen_neighbor_symbol = atoms[ oxygen_neighbor_index ].element
    
    print("Moj somsiad: ", oxygen_neighbor_symbol)
    
    center_neighbors = graph.neighbors( oxygen_neighbor_index )
    
    oxygens_found = 0
    for unknown_atom in center_neighbors:
        if atoms[unknown_atom].element == "O":
            oxygens_found+=1
            
    if oxygen_neighbor_symbol == "C" and oxygens_found < 2:
        return False, oxygen_neighbor_symbol+str(oxygens_found)+"O"
    
    return True, oxygen_neighbor_symbol+str(oxygens_found)+"O"
    
    

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
    
    if len(atoms) > 1:
        return False, atom.element
    
    return True, atom.element

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
    
    if len(atoms) > 1:
        return False, atom.element
    
    return True, atom.element

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
    
    G = nx.Graph()
    atoms_found = []
    for atom1Ind in range(len(atoms)):
        atom1 = atoms[atom1Ind]
        if atom1.element == "H":
            continue
        
        if atom1.get_fullname() == atomName:
            atoms_found.append(atom1Ind)
        
        for atom2Ind in range(atom1Ind+1, len(atoms)):
            atom2 = atoms[atom2Ind]            
            if atom2.element == "H":
                continue
            
            distance = atom1 - atom2
            
            if distance < 1.8 :
                G.add_edge(atom1Ind, atom2Ind)
                
    if len(atoms_found) != 1:
        print("WTF!? ", atoms_found)
        
    return G, atoms_found[0]

if __name__ == "__main__":
    findSupramolecularAnionPiLigand( "LUM", "cif/1he5.cif" )
    findSupramolecularAnionPiLigand( "NAP", "cif/3bcj.cif" )
    findSupramolecularAnionPiLigand( "NAP", "cif/4lbs.cif" )
#findSupramolecular( "LUM", 666 )