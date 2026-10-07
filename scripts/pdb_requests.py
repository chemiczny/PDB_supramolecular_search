"""
Written in collaboration with the most beautiful woman in the world, born
in Olkusz, the lady of my heart, Emilia Kuzniak.

Finds PDB codes of structures containing ligands listed in an SDF file.
"""

import requests

ses = requests.Session()

"""Helper class holding information about a single PDB entry."""


class PDBdata:
    def __init__(self, pdb_code, method):
        self.pdb_code = pdb_code
        self.method = method


def get_pdb_code(ligand_code):
    """Find unique PDB entries containing the given ligand.

    Args:
        ligand_code: ligand code.

    Returns:
        Tuple (ligand2pdb_data, ligand_stack): a dict mapping ligand codes to
        lists of PDBdata objects (duplicates removed), and the list of ligand
        codes in the order they were found. All codes refer to the same
        ligand (codes may have been replaced over time).
    """
    ligand2_pdb_codes, ligand_stack = get_pdb_code_and_replaced_ligand_code(ligand_code)

    ligand2pdb_data = {}

    for ligand in ligand_stack:
        pdbs_processed = []
        pdb_data = get_all_pdb_codes(ligand)

        # If anything was found for the current ligand code,
        # create a list for it in the results dict
        if len(pdb_data) > 0:
            ligand2pdb_data[ligand] = []
        elif ligand in ligand2_pdb_codes:
            if len(ligand2_pdb_codes[ligand]) > 0:
                ligand2pdb_data[ligand] = []

        # Store the data found in the results dict, skipping duplicates
        for pdb in pdb_data:
            if pdb.pdb_code not in pdbs_processed:
                ligand2pdb_data[ligand].append(pdb)
                pdbs_processed.append(pdb.pdb_code)

        if ligand in ligand2_pdb_codes:
            for pdb in ligand2_pdb_codes[ligand]:
                if pdb not in pdbs_processed:
                    ligand2pdb_data[ligand].append(PDBdata(pdb, "Unknown"))
                    pdbs_processed.append(pdb)
                    print("Unknown method! ", ligand, pdb)

    return ligand2pdb_data, ligand_stack


def get_all_pdb_codes(ligand_code):
    """Find all PDB entries containing the given ligand, using Ligand Expo.

    Duplicates are not removed.

    Args:
        ligand_code: ligand code.

    Returns:
        List of PDBdata objects.
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
    """Find PDB codes for a ligand, following "Replaced by" links in Ligand Expo.

    Args:
        ligand_code: ligand code.

    Returns:
        Tuple (ligand2_pdb_codes, ligand_stack): a dict mapping ligand codes to
        lists of PDB codes, and the list of ligand codes in the order they were
        found. All codes refer to the same ligand.
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
                print("Adding a new PDB code to an existing record")
                ligand2_pdb_codes[ligand_code].append(pdb_code)
            else:
                ligand2_pdb_codes[ligand_code] = [pdb_code]

        elif "Replaced by" in line:
            print("Found 'Replaced by' for " + ligand_code)
            line_with_ligand_code = html_text_spl[found_ind + 1]
            new_ligand_code = line_with_ligand_code.split(">")[1].split("<")[0]
            new_ligand2_pdb_code, new_ligand_stack = (
                get_pdb_code_and_replaced_ligand_code(new_ligand_code)
            )
            ligand_stack += new_ligand_stack

            for new_ligand in new_ligand2_pdb_code:
                if new_ligand in ligand2_pdb_codes:
                    print("Duplicate ligand code after code replacement")
                    ligand2_pdb_codes[new_ligand] += new_ligand2_pdb_code[new_ligand]
                else:
                    ligand2_pdb_codes[new_ligand] = new_ligand2_pdb_code[new_ligand]

        found_ind += 1

    return ligand2_pdb_codes, ligand_stack


def get_ligand_code_from_sdf(sdf_file_name):
    """Read ligand codes from an SDF file.

    Args:
        sdf_file_name: path to the SDF file.

    Returns:
        List of ligand codes.
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
    """Append the data found for a single ligand to the output file.

    The output file name is given by the global variable sdf_output.

    Args:
        ligands2_pdb_data: dict mapping ligand codes to lists of PDBdata objects.
        ligands_stack: list of ligand codes referring to the same ligand.
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
    """Append a message about a ligand with no PDB entries to PDBwrong.log."""
    output_name = "PDBwrong.log"
    output_file = open(output_name, "a+")

    output_file.write(text + "\n\n")

    output_file.close()


if __name__ == "__main__":
    """
    Main part:
    1. Read ligand codes from the SDF file.
    2. For each ligand, find its PDB codes and write them to the output file.
    """
    # sdfInput =  "sdf/aromaty_wiecej_niz_1_pierscien_podst_elektrofilowe_2.sdf"

    sdf_input = "sdf/wiecej_niz_1_pierscien_obecny_aromat_i_metal.sdf"
    sdf_output = sdf_input[0:-3] + "log"

    sdf_ligand_codes = get_ligand_code_from_sdf(sdf_input)
    ligands_no = len(sdf_ligand_codes)
    print("Found: " + str(ligands_no) + " ligand codes")

    ligand_ind = 0

    for ligand in sdf_ligand_codes:
        ligands2_pdb_data, ligands_stack = get_pdb_code(ligand)
        ligand_ind += 1

        if len(ligands_stack) > 2:
            print("More than two ligands on the stack! " + str(ligands_stack))

        if ligands2_pdb_data:
            add_data_to_output(ligands2_pdb_data, ligands_stack)
        else:
            add_text_to_output("Cannot find a PDB code for: " + ligand)

        if ligand_ind % 20 == 0:
            print("Progress: " + str(ligand_ind) + "/" + str(ligands_no))
