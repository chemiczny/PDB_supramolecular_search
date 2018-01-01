"""
Skrypt powstał we wspolpracy z najpiekniejsza kobieta na swiecie, Zrodzona
w Olkuszu, Pania mojego serca Emilia Kuzniak.

Sluzy do wyszukiwania kodow PDB na podstawie kodow ligandow.
"""
import requests

ses = requests.Session()

"""
Klasa pomocnicza do przechowywania informacji o konkretnym rekordzie PDB
"""
class PDBdata:
    def __init__(self, PDBcode, method):
        self.PDBcode = PDBcode
        self.method = method


def getPdbCode( ligandCode ):
    """
    Funkcja do znajdowania unikalnych rekordow PDB powiazanych z wprowadzonym 
    kodem liganda. Z znalezionych wynikow usuwane sa duplikaty
    
    Wejscie:
    ligandcode - string, kod liganda
    
    Wyjcie:
    ligand2pdbData - slownik, kluczem jest kod liganda, wartoscia lista obiektow
                    PDBdata
    ligandStack    - lista kodow ligandow, kolejnosc jest zgodna z kolejnoscia
                    znajdowania kodow w bazie. Wszystkie kody dotycza tej samej
                    struktury
    """
    ligand2PDBcodes, ligandStack = getPdbCodeAndReplacedLigandCode( ligandCode)
    
    ligand2pdbData = {}
    
    for ligand in ligandStack:
        pdbsProcessed = []
        pdbData = getAllPDBcodes( ligand )
        
        #Jesli dla obecnego kodu liganda znaleziono cokolwiek
        # to tworzymy d;a niego tablice w slowniku z wynikami
        if len( pdbData ) > 0:
            ligand2pdbData[ligand] = []
        elif ligand in ligand2PDBcodes:
            if len(ligand2PDBcodes[ligand]) > 0 :
                ligand2pdbData[ligand] = []
            
        # Do slownika z wynikami zapisujemy znalezione dane, dbamy
        # o nie zapisywanie duplikatow
        for pdb in pdbData:
            if not pdb.PDBcode in pdbsProcessed:
                ligand2pdbData[ligand].append( pdb )
                pdbsProcessed.append(pdb.PDBcode)
                
        if ligand in ligand2PDBcodes:
            for pdb in ligand2PDBcodes[ligand]:
                if not pdb in pdbsProcessed:
                    ligand2pdbData[ligand].append( PDBdata( pdb, "Unknown"  ) )
                    pdbsProcessed.append(pdb)
                    print("Nieznana metoda! ", ligand, pdb)
                
    return ligand2pdbData, ligandStack

def getAllPDBcodes( ligandCode ):
    """
    Funkcja sluzy do znajdowania rekordow PDB powiazanych z konretnym ligandem.
    Ze znalezionych rekordow nie sa usuwane duplikaty.
    Wyszukiwania sa dokonywane przez ... (nie wiem jak nazwac ta druga wyszukiwarke)
    
    Wejscie:
    ligandcode - string, kod liganda
    
    Wyjcie:
    ligand2pdbData - slownik, kluczem jest kod liganda, wartoscia lista obiektow
                    PDBdata
    """
    ligandSearchAdress = "http://ligand-expo.rcsb.org/pyapps/ldHandler.py"
    
    payload = {'formid': 'cc-db-inst-search' ,'targetId': ligandCode, 'operation':'idsearch'}
    
    r = ses.post(ligandSearchAdress, params=payload)
    
    htmlText = r.text
    htmlTextSpl = htmlText.split('\n')
    
    allPdbCodes = []
    lineInd = 0    
    PdbNo = 0
    PdbFound = 0
    
    for line in htmlTextSpl:
        if "No results found for this query" in line:
            break
        elif "Count in released entries" in line:
            PdbNo = int( line.split("td>")[3].split("<")[0] )
        elif "rs1-" in line:
            newPdbCode = htmlTextSpl[lineInd+1].split(">")[2].split("<")[0]
            newPdbCode = newPdbCode.upper().strip()
            newMethod  = htmlTextSpl[lineInd+3].split(">")[1].split("<")[0]
            
            allPdbCodes.append( PDBdata( newPdbCode, newMethod ) )
            PdbFound+=1
        
    
        lineInd+= 1
        
    if PdbFound != PdbNo:
        print("PdbFound != PdbNo")
        print("Pdbfound: ", PdbFound)
        print("PdbNo: ", PdbNo)
        print("Ligand code: ", ligandCode)
        
    return allPdbCodes

def getPdbCodeAndReplacedLigandCode( ligandCode ):
    """
    Funkcja sluzy do znajdowania rekordow PDB powiazanych z konretnym ligandem.
    Wyszukiwania sa dokonywane przez ... (nie wiem jak nazwac ta druga wyszukiwarke)
    
    Wejscie:
    ligandcode - string, kod liganda
    
    Wyjcie:
    ligand2pdbCodes - slownik, kluczem jest kod liganda, wartoscia lista kodow
                    PDB
    
    ligandStack    - lista kodow ligandow, kolejnosc jest zgodna z kolejnoscia
                    znajdowania kodow w bazie. Wszystkie kody dotycza tej samej
                    struktury
    """
    ligandSearchAdress = "http://ligand-expo.rcsb.org/pyapps/ldHandler.py"

    payload = {'formid': 'cc-index-search' ,'target': ligandCode, 'operation':'ccid'}
    r = ses.get(ligandSearchAdress, params=payload)

    htmlText = r.text
    htmlTextSpl = htmlText.split('\n')

    foundInd = 0
    PDBCode= None
    
    ligand2PDBcodes = {}
    ligandStack = [ ligandCode ]
    for line in htmlTextSpl:
        if "Model PDB code" in line:
            lineWithPDBCode = htmlTextSpl[foundInd+1]
            PDBCode = lineWithPDBCode.split(">")[1].split("<")[0].strip()
            PDBCode = PDBCode.upper().strip()
            
            if ligandCode in ligand2PDBcodes:
                print("Dodaje nowy PDB code do istniejacego rekordu")
                ligand2PDBcodes[ligandCode].append( PDBCode )
            else:
                ligand2PDBcodes[ligandCode] = [ PDBCode ]

        elif "Replaced by" in line:
            print("Znaleziono 'Replaced by' dla "+ligandCode)
            lineWithLigandCode = htmlTextSpl[foundInd+1]
            newLigandCode = lineWithLigandCode.split(">")[1].split("<")[0]
            newLigand2PDBcode, newLigandStack = getPdbCodeAndReplacedLigandCode( newLigandCode )
            ligandStack += newLigandStack
            
            for newLigand in newLigand2PDBcode:
                if newLigand in ligand2PDBcodes:
                    print( "Powtarzajacy sie ligand code z zamiany kodu" )
                    ligand2PDBcodes[newLigand] +=  newLigand2PDBcode[ newLigand ] 
                else:
                    ligand2PDBcodes[newLigand] =  newLigand2PDBcode[ newLigand ] 

	    
        foundInd+=1
        
    return ligand2PDBcodes, ligandStack

def getLigandCodeFromSdf ( sdfFileName  ):
    """
    Funkcja sluzy do pobierania kodow ligandow z pliku .sdf
    
    Wejscie:
    sdfFileName - nazwa pliku sdf
    
    Wyjscie:
    ligandCodes - lista znalezionych kodow ligandow
    """
    sdfFile= open(sdfFileName, 'r' )

    line = sdfFile.readline()
    ligandCodes = []
    while line:
        if "field_0" in line:
            lineWithData = sdfFile.readline()
            ligandCodes.append( lineWithData.strip() )
            
        line = sdfFile.readline()

    sdfFile.close()

    return ligandCodes

def addDataToOutput( ligands2PDBdata, ligandsStack):
    """
    Funkcja sluzy do wpisania znalezionych danych dotyczacych pojedynczej struktury
    liganda (a moze ligandu?). Dane wspisywane sa do pliku wynikowego, jego nazwa
    jest okreslona przez wartosc zmiennej globalnej sdfOutput
    
    Wejscie:
    ligands2PDBdata - slownik, kluczem jest kod liganda, wartoscia lista obiektow
                    PDBdata
    ligandStack     - lista znalezionych kodow odpowiadajacych pojedynczemu 
                    ligandowi.
    """
    outputName = sdfOutput
    outputFile = open( outputName, "a+" )


    firstLigand = True
    for ligand in ligandsStack:
        if not firstLigand:
            outputFile.write(" Replaced by: ")
        outputFile.write(" ligand code: "+ligand+": " )
        if ligand in ligands2PDBdata:
            for pdbData in ligands2PDBdata[ ligand ]:
                outputFile.write(" pdb code: "+ pdbData.PDBcode+" method: "+pdbData.method )
            
        firstLigand = False
            
    outputFile.write("\n\n")

    outputFile.close()

def addTextToOutput( text  ):
    """
    Funckja zapisuje dane na temat ligandow, dla ktorych nie znaleziono ani jednego
    rekordu w bazie PDB. Zapis nastepuje do pliku PDBwrong.log
    """
    outputName = "PDBwrong.log"
    outputFile = open( outputName, "a+" )

    outputFile.write(text+"\n\n")

    outputFile.close()

"""
Wlasciwa czesc kodu:
1. Pobierz kody ligandow z pliku sdf
2. Dla kazdego z ligandow znajdz kody PDB i zapisz je do pliku
"""
#sdfInput =  "aromaty_wiecej_niz_1_pierscien_podst_elektrofilowe_2.sdf"

sdfInput = "wiecej_niz_1_pierscien_obecny_aromat_i_metal.sdf"
sdfOutput = sdfInput[0:-3]+"log"

ligandyEmilki = getLigandCodeFromSdf( sdfInput )
ligandsNo = len(ligandyEmilki)
print("Znaleziono: "+str(ligandsNo)+" kodow ligandow")

ligandInd = 0

for ligand in ligandyEmilki:
    ligands2PDBdata, ligandsStack = getPdbCode( ligand )
    ligandInd+= 1
    
    if len(ligandsStack)>2:
        print("Wiecej niz dwa ligandy na stosie! "+str(ligandsStack))
    
    if ligands2PDBdata:
        addDataToOutput( ligands2PDBdata, ligandsStack)
    else:
        addTextToOutput("Nie mozna znalezc pdb code dla: "+ligand)

    if ligandInd % 20 == 0:
        print("Postep: "+str(ligandInd)+"/"+str(ligandsNo))
