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
    def __init__(self, pdb_code, method):
        self.pdb_code = pdb_code
        self.method = method


def get_pdb_code(ligand_code):
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
    ligand2_pdb_codes, ligand_stack = get_pdb_code_and_replaced_ligand_code(ligand_code)

    ligand2pdb_data = {}

    for ligand in ligand_stack:
        pdbs_processed = []
        pdb_data = get_all_pdb_codes(ligand)

        # Jesli dla obecnego kodu liganda znaleziono cokolwiek
        # to tworzymy d;a niego tablice w slowniku z wynikami
        if len(pdb_data) > 0:
            ligand2pdb_data[ligand] = []
        elif ligand in ligand2_pdb_codes:
            if len(ligand2_pdb_codes[ligand]) > 0:
                ligand2pdb_data[ligand] = []

        # Do slownika z wynikami zapisujemy znalezione dane, dbamy
        # o nie zapisywanie duplikatow
        for pdb in pdb_data:
            if pdb.pdb_code not in pdbs_processed:
                ligand2pdb_data[ligand].append(pdb)
                pdbs_processed.append(pdb.pdb_code)

        if ligand in ligand2_pdb_codes:
            for pdb in ligand2_pdb_codes[ligand]:
                if pdb not in pdbs_processed:
                    ligand2pdb_data[ligand].append(PDBdata(pdb, "Unknown"))
                    pdbs_processed.append(pdb)
                    print("Nieznana metoda! ", ligand, pdb)

    return ligand2pdb_data, ligand_stack


def get_all_pdb_codes(ligand_code):
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
    ligand_search_adress = "http://ligand-expo.rcsb.org/pyapps/ldHandler.py"

    payload = {
        "formid": "cc-db-inst-search",
        "targetId": ligand_code,
        "operation": "idsearch",
    }

    r = ses.post(ligand_search_adress, params=payload)

    html_text = r.text
    html_text_spl = html_text.split("\n")

    all_pdb_codes = []
    line_ind = 0
    pdb_no = 0
    pdb_found = 0

    for line in html_text_spl:
        if "No results found for this query" in line:
            break
        elif "Count in released entries" in line:
            pdb_no = int(line.split("td>")[3].split("<")[0])
        elif "rs1-" in line:
            new_pdb_code = html_text_spl[line_ind + 1].split(">")[2].split("<")[0]
            new_pdb_code = new_pdb_code.upper().strip()
            new_method = html_text_spl[line_ind + 3].split(">")[1].split("<")[0]

            all_pdb_codes.append(PDBdata(new_pdb_code, new_method))
            pdb_found += 1

        line_ind += 1

    if pdb_found != pdb_no:
        print("PdbFound != PdbNo")
        print("Pdbfound: ", pdb_found)
        print("PdbNo: ", pdb_no)
        print("Ligand code: ", ligand_code)

    return all_pdb_codes


def get_pdb_code_and_replaced_ligand_code(ligand_code):
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
    ligand_search_adress = "http://ligand-expo.rcsb.org/pyapps/ldHandler.py"

    payload = {"formid": "cc-index-search", "target": ligand_code, "operation": "ccid"}
    r = ses.get(ligand_search_adress, params=payload)

    html_text = r.text
    html_text_spl = html_text.split("\n")

    found_ind = 0
    pdb_code = None

    ligand2_pdb_codes = {}
    ligand_stack = [ligand_code]
    for line in html_text_spl:
        if "Model PDB code" in line:
            line_with_pdb_code = html_text_spl[found_ind + 1]
            pdb_code = line_with_pdb_code.split(">")[1].split("<")[0].strip()
            pdb_code = pdb_code.upper().strip()

            if ligand_code in ligand2_pdb_codes:
                print("Dodaje nowy PDB code do istniejacego rekordu")
                ligand2_pdb_codes[ligand_code].append(pdb_code)
            else:
                ligand2_pdb_codes[ligand_code] = [pdb_code]

        elif "Replaced by" in line:
            print("Znaleziono 'Replaced by' dla " + ligand_code)
            line_with_ligand_code = html_text_spl[found_ind + 1]
            new_ligand_code = line_with_ligand_code.split(">")[1].split("<")[0]
            new_ligand2_pdb_code, new_ligand_stack = (
                get_pdb_code_and_replaced_ligand_code(new_ligand_code)
            )
            ligand_stack += new_ligand_stack

            for new_ligand in new_ligand2_pdb_code:
                if new_ligand in ligand2_pdb_codes:
                    print("Powtarzajacy sie ligand code z zamiany kodu")
                    ligand2_pdb_codes[new_ligand] += new_ligand2_pdb_code[new_ligand]
                else:
                    ligand2_pdb_codes[new_ligand] = new_ligand2_pdb_code[new_ligand]

        found_ind += 1

    return ligand2_pdb_codes, ligand_stack


def get_ligand_code_from_sdf(sdf_file_name):
    """
    Funkcja sluzy do pobierania kodow ligandow z pliku .sdf

    Wejscie:
    sdfFileName - nazwa pliku sdf

    Wyjscie:
    ligandCodes - lista znalezionych kodow ligandow
    """
    sdf_file = open(sdf_file_name, "r")

    line = sdf_file.readline()
    ligand_codes = []
    while line:
        if "field_0" in line:
            line_with_data = sdf_file.readline()
            ligand_codes.append(line_with_data.strip())

        line = sdf_file.readline()

    sdf_file.close()

    return ligand_codes


def add_data_to_output(ligands2_pdb_data, ligands_stack):
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
    output_name = sdf_output
    output_file = open(output_name, "a+")

    first_ligand = True
    for ligand in ligands_stack:
        if not first_ligand:
            output_file.write(" Replaced by: ")
        output_file.write(" ligand code: " + ligand + ": ")
        if ligand in ligands2_pdb_data:
            for pdb_data in ligands2_pdb_data[ligand]:
                output_file.write(
                    " pdb code: " + pdb_data.pdb_code + " method: " + pdb_data.method
                )

        first_ligand = False

    output_file.write("\n\n")

    output_file.close()


def add_text_to_output(text):
    """
    Funckja zapisuje dane na temat ligandow, dla ktorych nie znaleziono ani jednego
    rekordu w bazie PDB. Zapis nastepuje do pliku PDBwrong.log
    """
    output_name = "PDBwrong.log"
    output_file = open(output_name, "a+")

    output_file.write(text + "\n\n")

    output_file.close()


if __name__ == "__main__":
    """
    Wlasciwa czesc kodu:
    1. Pobierz kody ligandow z pliku sdf
    2. Dla kazdego z ligandow znajdz kody PDB i zapisz je do pliku
    """
    # sdfInput =  "sdf/aromaty_wiecej_niz_1_pierscien_podst_elektrofilowe_2.sdf"

    sdf_input = "sdf/wiecej_niz_1_pierscien_obecny_aromat_i_metal.sdf"
    sdf_output = sdf_input[0:-3] + "log"

    ligandy_emilki = get_ligand_code_from_sdf(sdf_input)
    ligands_no = len(ligandy_emilki)
    print("Znaleziono: " + str(ligands_no) + " kodow ligandow")

    ligand_ind = 0

    for ligand in ligandy_emilki:
        ligands2_pdb_data, ligands_stack = get_pdb_code(ligand)
        ligand_ind += 1

        if len(ligands_stack) > 2:
            print("Wiecej niz dwa ligandy na stosie! " + str(ligands_stack))

        if ligands2_pdb_data:
            add_data_to_output(ligands2_pdb_data, ligands_stack)
        else:
            add_text_to_output("Nie mozna znalezc pdb code dla: " + ligand)

        if ligand_ind % 20 == 0:
            print("Postep: " + str(ligand_ind) + "/" + str(ligands_no))
