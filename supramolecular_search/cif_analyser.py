"""
Created on Sun Dec 17 13:56:41 2017

Written in collaboration with the most beautiful woman in the world, born
in Olkusz, the lady of my heart, Emilia Kuzniak.

Analyses a single CIF file to find supramolecular interactions involving
aromatic rings (anion-pi, cation-pi, pi-pi, ...) and computes their geometric
parameters for further analysis.
"""

import logging
from os.path import getsize

from supramolecular_search.config import configure

config = configure()

from Bio.PDB import FastMMCIFParser, NeighborSearch, Selection
from supramolecular_search.primitive_cif2dict import PrimitiveCif2Dict
import numpy as np
from supramolecular_search.supramolecular_logging import SupramolecularLogger
from supramolecular_search.ring_detection import (
    get_rings_centroids,
    find_in_graph,
    is_flat_primitive,
    normalize,
    molecule2graph,
)
from supramolecular_search.protonate import Protonate
from supramolecular_search.anion_recogniser import AnionRecogniser, create_res_id
from multiprocessing import current_process
import networkx as nx
from collections import defaultdict
from time import time


logger = logging.getLogger(__name__)


def find_supramolecular(cif_data):
    """Analyse a single CIF file and write the interactions found to the logs.

    Args:
        cif_data: tuple (cif_file, pdb_code, log_id) passed to CifAnalyser.
    """

    cif_analyser = CifAnalyser(*cif_data)
    cif_analyser.analyse_cif()


class CifAnalyser:
    def __init__(self, cif_file, pdb_code, log_id="default"):
        global config
        self.cif_file = cif_file
        self.pdb_code = pdb_code

        if log_id == "default":
            self.file_id = current_process().name
        else:
            self.file_id = log_id

        self.ns = None
        self.resolution = None
        self.method = None
        self.h_atoms_present = False

        self.structure_type = "Unknown"

        self.anion_recogniser = AnionRecogniser()

        self.aa_cation_radius = 5.0
        self.metal_cation_radius = 10
        self.h_bonds_radius = 3.5

        self.small_cutting_radius = 5.0
        self.big_cutting_radius = 12

        self.aromatic_aa_counter = {}
        self.aromatic_aa = ["PHE", "HIS", "TRP", "TYR"]

        self.supra_logger = SupramolecularLogger(
            pdb_code, self.file_id, config["scratch"]
        )

    def init_aromatic_aa_counter(self):
        self.aromatic_aa_counter = {}

        for aa_code in self.aromatic_aa:
            self.aromatic_aa_counter[aa_code] = 0

    def analyse_cif(self):
        self.supra_logger.write_additional_info("Starting analysis: " + self.pdb_code)
        parser = FastMMCIFParser(QUIET=True)

        time_start = time()
        try:
            time1 = time()
            structure = parser.get_structure("temp", self.cif_file)
            time_bio_parsing = time() - time1
        except Exception:
            error_message = "FastMMCIFParser cannot handle with: " + self.cif_file
            self.supra_logger.increment_partial_progress()
            self.supra_logger.write_additional_info(error_message)
            return True

        self.supra_logger.write_additional_info(
            "File size: " + str(getsize(self.cif_file))
        )

        not_piacids = [
            "HOH",
            "DOD",
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PRO",
            "SER",
            "THR",
            "VAL",
        ]

        time2 = time()
        self.resolution, self.method = self.read_resolution_and_method()
        time_resolution_reading = time() - time2

        self.h_atoms_present = hydrogens_present(structure)
        self.structure_type = self.determine_structure_type(structure)

        times = {
            "getCentroidTime": 0,
            "extractAnionsAroundRingTime": 0,
            "extractCationsTime": 0,
            "cationAnalysisTime": 0,
            "extractPiPiTime": 0,
            "extractHBondsTime": 0,
            "timeBioParsing": time_bio_parsing,
            "timeResolutionReading": time_resolution_reading,
        }

        for model_index, model in enumerate(structure):
            self.init_aromatic_aa_counter()
            residue2counts = defaultdict(int)

            atoms = Selection.unfold_entities(model, "A")
            not_disordered_atoms = []
            for atom in atoms:
                if not atom.is_disordered() or atom.get_altloc() == "A":
                    not_disordered_atoms.append(atom)

            if len(not_disordered_atoms) == 0:
                continue
            self.ns = NeighborSearch(not_disordered_atoms)
            for residue in model.get_residues():
                residue_name = residue.get_resname().upper()
                first_atom = list(residue.get_atoms())[0]
                not_asolution = (
                    first_atom.is_disordered() and first_atom.get_altloc() != "A"
                )

                if not not_asolution:
                    residue2counts[residue_name] += 1

                if residue_name not in not_piacids:
                    supramolecular_found, new_times = self.analyse_piacid(
                        residue, model_index, structure
                    )
                    for key in new_times:
                        times[key] += new_times[key]

            self.supra_logger.write_additional_info("model: " + str(model_index))
            for rc in residue2counts:
                self.supra_logger.write_additional_info(
                    rc + " " + str(residue2counts[rc])
                )

            self.supra_logger.write_additional_info("Recognised aromatic AA")
            for aa_code in self.aromatic_aa_counter:
                self.supra_logger.write_additional_info(
                    aa_code + " " + str(self.aromatic_aa_counter[aa_code])
                )

        self.supra_logger.increment_partial_progress()
        time_stop = time()
        time_taken = time_stop - time_start

        for key in times:
            self.supra_logger.write_additional_info(
                "time " + key + " " + str(times[key])
            )

        self.supra_logger.write_additional_info(
            "Analysis finished " + self.pdb_code + " time: " + str(time_taken)
        )

    def determine_structure_type(self, structure):
        all_aa = [
            "ALA",
            "CYS",
            "GLY",
            "ILE",
            "LEU",
            "MET",
            "ASN",
            "PRO",
            "GLN",
            "SER",
            "THR",
            "VAL",
            "ASP",
            "GLU",
            "PHE",
            "HIS",
            "TRP",
            "TYR",
            "LYS",
            "ARG",
        ]
        dna_nucleotides = ["DA", "DC", "DG", "DT", "DI"]
        rna_nucleotides = ["A", "G", "T", "C", "U", "I"]

        aa_counter = 0
        dna_counter = 0
        rna_counter = 0

        for res in structure.get_residues():
            res_name = res.get_resname().upper()

            if res_name in all_aa:
                aa_counter += 1
            elif res_name in dna_nucleotides:
                dna_counter += 1
            elif res_name in rna_nucleotides:
                rna_counter += 1

        structure_type = []

        if aa_counter > 0:
            structure_type.append("protein")

        if dna_counter > 0:
            structure_type.append("DNA")

        if rna_counter > 0:
            structure_type.append("RNA")

        if not structure_type:
            return "other"

        return "-".join(structure_type)

    def read_resolution_and_method(self):
        try:
            mmcif_dict_parser = PrimitiveCif2Dict(
                self.cif_file,
                ["_refine.ls_d_res_high", "_reflns_shell.d_res_high", "_exptl.method"],
            )
            mmcif_dict = mmcif_dict_parser.result
        except Exception:
            error_message = "primitiveCif2Dict cannot handle with: " + self.cif_file
            self.supra_logger.write_additional_info(error_message)
            return -666, -666, "Unknown"

        method = "Unknown"
        if "_exptl.method" in mmcif_dict:
            method = mmcif_dict["_exptl.method"]
            if len(method) == 1:
                method = method[0]
            else:
                method = sorted(method)

        res_key = "_refine.ls_d_res_high"
        resolution = []
        if res_key in mmcif_dict:
            resolution = mmcif_dict[res_key]

        if not resolution:
            return -1, method

        if len(resolution) == 1:
            if isfloat(resolution[0]):
                return resolution[0], method
            else:
                return -1, method
        else:
            res_floats = []
            for res in resolution:
                if isfloat(res):
                    res_floats.append(float(res))

            if not res_floats:
                return -1, method

            res_floats = sorted(res_floats)
            return res_floats[0], method

    def analyse_piacid(self, ligand, model_index, structure):
        first_atom = list(ligand.get_atoms())[0]
        if first_atom.is_disordered() and first_atom.get_altloc() != "A":
            return False, {}

        get_centroid_time = 0
        extract_anions_around_ring_time = 0
        extract_cations_time = 0
        cation_analysis_time = 0
        extract_pi_pi_time = 0
        extract_h_bonds_time = 0

        time0 = time()
        centroids, ligand_graph = get_rings_centroids(ligand, True)
        get_centroid_time += time() - time0

        ligand_with_anions = False
        ligand_code = ligand.get_resname().upper()

        for centroid in centroids:
            self.anion_recogniser.clean_properties_pack()

            if ligand_code in self.aromatic_aa:
                self.aromatic_aa_counter[ligand_code] += 1

            big_distance = self.big_cutting_radius
            distance = self.small_cutting_radius

            atoms = list(
                self.ns.search(np.array(centroid["coords"]), big_distance, "A")
            )

            ns_small = NeighborSearch(atoms)
            neighbors = list(
                ns_small.search(np.array(centroid["coords"]), distance, "A")
            )

            time1 = time()
            extracted_anion_atoms = self.anion_recogniser.extract_anion_atoms(
                neighbors, ligand, ns_small
            )
            extract_anions_around_ring_time += time() - time1
            methyls = extract_aa_methyls(neighbors)

            if len(extracted_anion_atoms) > 0 or True:
                self.write_geometric_properties(ligand, centroid, model_index)

                time2 = time()
                extracted_metal_cations = extract_metal_cations(
                    centroid["coords"], ns_small, self.metal_cation_radius
                )
                extracted_aa_cations = extract_aa_cations(
                    centroid["coords"], ns_small, self.aa_cation_radius
                )
                extract_cations_time += time() - time2

                cation_ring_len_chains = []
                cation_complex_data = []

                time3 = time()
                for cat in extracted_metal_cations:
                    cation_ring_len_chains.append(
                        find_chain_len_cation_ring(
                            cat, ligand, centroid, self.ns, ligand_graph, self.file_id
                        )
                    )
                    cation_complex_data.append(
                        self.find_cation_complex(cat, self.ns, ligand)
                    )
                cation_analysis_time += time() - time3

                extracted_cations = extracted_metal_cations + extracted_aa_cations
                cation_ring_len_chains_full = cation_ring_len_chains + len(
                    extracted_aa_cations
                ) * [(0, False, -1)]

                self.supra_logger.write_cation_pi_results(
                    ligand,
                    centroid,
                    extracted_cations,
                    cation_ring_len_chains_full,
                    model_index,
                )
                self.supra_logger.write_metal_ligand_results(
                    extracted_metal_cations, cation_complex_data, model_index
                )

                time4 = time()
                extracted_centroids, ring_molecules = extract_ring_centroids(
                    centroid["coords"], ligand, ns_small
                )
                extract_pi_pi_time += time() - time4

                self.supra_logger.write_pi_pi_results(
                    ligand, centroid, ring_molecules, extracted_centroids, model_index
                )

                for atom in extracted_anion_atoms:
                    self.supra_logger.write_anion_cation_results(
                        atom["Atom"], ligand, centroid, extracted_cations, model_index
                    )

                    if atom["Atom"].anionData.h_bonds_analyzed:
                        continue

                    time5 = time()
                    h_donors = extract_hbonds(
                        atom,
                        ns_small,
                        self.h_bonds_radius,
                        self.h_atoms_present,
                        self.file_id,
                        structure,
                    )
                    extract_h_bonds_time += time() - time5

                    self.supra_logger.write_hbonds_results(h_donors, atom, model_index)
                    atom["Atom"].anionData.h_bonds_analyzed = True

            extracted_atoms = self.supra_logger.write_anion_pi_results(
                ligand,
                centroid,
                extracted_anion_atoms,
                model_index,
                self.resolution,
                self.method,
                self.structure_type,
            )

            if methyls:
                self.supra_logger.write_methyl_pi_results(
                    ligand,
                    centroid,
                    methyls,
                    model_index,
                    self.resolution,
                    self.method,
                    self.structure_type,
                )

            if len(extracted_atoms) > 0:
                ligand_with_anions = True

        times = {
            "getCentroidTime": get_centroid_time,
            "extractAnionsAroundRingTime": extract_anions_around_ring_time,
            "extractCationsTime": extract_cations_time,
            "cationAnalysisTime": cation_analysis_time,
            "extractPiPiTime": extract_pi_pi_time,
            "extractHBondsTime": extract_h_bonds_time,
        }

        return ligand_with_anions, times

    def write_geometric_properties(self, ligand, centroid, model_index):

        for geometric_property in self.anion_recogniser.properties2calculate_pack:
            if geometric_property.kind == "plane":
                self.supra_logger.write_anion_pi_planar_results(
                    ligand,
                    centroid,
                    geometric_property,
                    model_index,
                    geometric_property.anion_group_id,
                )
            elif geometric_property.kind == "line":
                self.supra_logger.write_anion_pi_linear_results(
                    ligand,
                    centroid,
                    geometric_property,
                    model_index,
                    geometric_property.anion_group_id,
                )
            elif geometric_property.kind == "lineSymmetric":
                self.supra_logger.write_anion_pi_linear_results(
                    ligand,
                    centroid,
                    geometric_property,
                    model_index,
                    geometric_property.anion_group_id,
                    True,
                )

    def find_cation_complex(self, cation, ns, ligand):
        if hasattr(cation, "analysedAsComplex"):
            return {"complex": False, "coordNo": 0, "ligands": []}

        cation.analysedAsComplex = True
        potential_ligands = ns.search(cation.get_coord(), 2.85, "A")
        anion_space_with_cation = ns.search(cation.get_coord(), 4.5, "A")
        metals = [
            "LI",
            "BE",
            "NA",
            "MG",
            "AL",
            "K",
            "CA",
            "SC",
            "TI",
            "V",
            "CR",
            "MN",
            "FE",
            "CO",
            "NI",
            "CU",
            "ZN",
            "GA",
            "GE",
            "AS",
            "RB",
            "SR",
            "Y",
            "ZR",
            "NB",
            "MO",
            "TC",
            "RU",
            "RH",
            "PD",
            "AG",
            "CD",
            "IN",
            "SN",
            "SB",
            "TE",
            "CS",
            "BA",
            "LA",
            "CE",
            "PR",
            "ND",
            "PM",
            "SM",
            "EU",
            "GD",
            "TB",
            "DY",
            "HO",
            "ER",
            "TM",
            "YB",
            "LU",
            "HF",
            "TA",
            "W",
            "RE",
            "OS",
            "IR",
            "PT",
            "AU",
            "HG",
            "TL",
            "PB",
            "BI",
            "PO",
            "AT",
            "FR",
            "RA",
            "AC",
            "TH",
            "PA",
            "U",
            "NP",
            "PU",
            "AM",
            "CM",
            "BK",
            "CF",
            "ES",
            "FM",
            "MD",
            "NO",
            "LR",
            "RF",
            "DB",
            "SG",
            "BH",
            "HS",
            "MT",
            "DS",
            "RG",
        ]

        anion_space = []
        for a in anion_space_with_cation:
            if a.element.upper() not in metals:
                anion_space.append(a)

        if len(anion_space) == 0:
            return {"complex": False, "coordNo": 0, "ligands": []}

        class FakeRes:
            def __init__(self):
                pass

            def get_resname(self):
                return "test_resname"

            def get_id(self):
                return "###"

            def get_atoms(self):
                return [FakeRes()]

            def get_parent(self):
                return FakeRes()

            def get_coord(self):
                return [0, 0, 0]

        resultant_vector = np.array([0.0, 0.0, 0.0])
        coord_no = 0
        ligands = []

        summary = {}

        self.anion_recogniser.extract_anion_atoms(potential_ligands, FakeRes(), ns)

        for atom in potential_ligands:
            if atom.element == "H":
                continue
            if atom != cation:
                new_vector = normalize(atom.get_coord() - cation.get_coord())
                resultant_vector += new_vector
                coord_no += 1
                if hasattr(atom, "anionData"):
                    is_anion = atom.anionData.charged
                    anion_type = atom.anionData.anion_type
                else:
                    is_anion = False
                    anion_type = atom.element

                anion_id = -1
                if is_anion:
                    anion_id = atom.anionData.anion_id

                potental_ligand_name = atom.get_parent().get_resname()
                if potental_ligand_name in summary:
                    summary[potental_ligand_name] += 1
                else:
                    summary[potental_ligand_name] = 1

                ligands.append(
                    {
                        "isAnion": is_anion,
                        "anionType": anion_type,
                        "atom": atom,
                        "AnionId": anion_id,
                    }
                )

        str_summary = []
        for key in sorted(list(summary.keys())):
            if summary[key] > 1:
                str_summary.append(key + "-" + str(summary[key]))
            else:
                str_summary.append(key)

        str_summary = cation.element + "_" + "_".join(str_summary)

        if coord_no == 0:
            return {
                "complex": False,
                "coordNo": 0,
                "ligands": [],
                "summary": str_summary,
            }

        vector_len = np.linalg.norm(resultant_vector)

        if vector_len < 0.2:
            return {
                "complex": True,
                "coordNo": coord_no,
                "ligands": ligands,
                "summary": str_summary,
            }
        else:
            return {
                "complex": False,
                "coordNo": coord_no,
                "ligands": ligands,
                "summary": str_summary,
            }


def extract_metal_cations(point, ns, distance):
    neighbors = ns.search(np.array(point), distance, "A")

    metal_cations_found = []

    metal_cations = [
        "Li",
        "Be",
        "Na",
        "Mg",
        "Al",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U",
        "Np",
        "Pu",
        "Am",
        "Cm",
        "Bk",
        "Cf",
        "Es",
        "Fm",
        "Md",
        "No",
        "Lr",
        "Rf",
        "Db",
        "Sg",
        "Bh",
        "Hs",
        "Mt",
        "Ds",
        "Rg",
    ]

    for atom in neighbors:
        if atom.element.upper() in [element.upper() for element in metal_cations]:
            metal_cations_found.append(atom)

    return metal_cations_found


def extract_aa_cations(point, ns, distance):
    neighbors = ns.search(np.array(point), distance, "A")

    aa_cations_found = []

    for atom in neighbors:
        if (
            atom.get_parent().get_resname().upper() in ["ARG", "LYS"]
            and atom.get_name() != "N"
            and atom.element.upper() == "N"
        ):
            aa_cations_found.append(atom)

    return aa_cations_found


def extract_aa_methyls(neighbors):
    methyls = []

    aa2methyls = {
        "ALA": ["CB"],
        "ILE": ["CG2", "CD1"],
        "LEU": ["CD1", "CD2"],
        "MET": ["CE"],
        "THR": ["CG2"],
        "VAL": ["CG1", "CG2"],
    }
    for atom in neighbors:
        resname = atom.get_parent().get_resname().upper()
        atom_name = atom.get_name()
        if resname in aa2methyls:
            if atom_name in aa2methyls[resname]:
                methyls.append({"Atom": atom, "AnionType": atom_name, "AnionId": 0})

    return methyls


def extract_hbonds(atom, ns_small, distance, h_atoms_present, file_id, structure):

    neighbors = ns_small.search(np.array(atom["Atom"].get_coord()), distance + 2, "A")
    acceptor_residue = atom["Atom"].get_parent()

    if not h_atoms_present:
        protonation_worker = Protonate()
        protonation_worker.protonate(neighbors, atom["Atom"])
        graph = protonation_worker.molecule_graph
        neighbors = protonation_worker.atom_list
        anion_atom_ind = protonation_worker.anion_id
    else:
        graph, anion_atom_ind = molecule2graph(
            neighbors, atom["Atom"], False, False, omit_metals=True
        )

    donors = []

    for potential_donor_ind in graph.nodes():
        element = neighbors[potential_donor_ind].element

        if element in ["O", "N"]:
            connected = list(graph.neighbors(potential_donor_ind))

            if atom["Atom"] - neighbors[potential_donor_ind] > distance:
                continue

            if neighbors[potential_donor_ind].get_parent() == acceptor_residue:
                continue

            connected_elements = [neighbors[c].element for c in connected]
            connected_elements = list(set(connected_elements))

            for c in connected:
                element = neighbors[c].element
                if element != "H":
                    continue
                else:
                    donors.append(
                        {
                            "donor": neighbors[potential_donor_ind],
                            "hydrogen": neighbors[c],
                            "HFromExp": h_atoms_present,
                        }
                    )

    return donors


def find_chain_len_cation_ring(
    cation, pi_acid, centroid_data, ns, ligand_graph, file_id
):
    pi_acid_atoms = list(pi_acid.get_atoms())

    cation_neighbours = ns.search(cation.get_coord(), 2.85, "A")
    first_atom_in_ring = centroid_data["cycleAtoms"][0]
    shortest_path = []
    distance_cat_chain = -1
    something_found = False
    for cat_n in cation_neighbours:
        if cat_n.element == "H":
            continue
        if cat_n == cation:
            continue
        if cat_n.get_parent() == pi_acid:
            temp_graph, cat_n_index = find_in_graph(ligand_graph, cat_n, pi_acid_atoms)
            if cat_n_index is None:
                logger.warning(
                    "Cation neighbour not found in ligand graph of %s %s chain %s",
                    pi_acid.get_resname(),
                    pi_acid.get_id(),
                    pi_acid.get_parent().get_id(),
                )
                continue

            if not nx.has_path(ligand_graph, first_atom_in_ring, cat_n_index):
                continue

            new_path = nx.shortest_path(ligand_graph, first_atom_in_ring, cat_n_index)
            new_path = [
                node for node in new_path if node not in centroid_data["cycleAtoms"]
            ]
            if len(new_path) < len(shortest_path) or not something_found:
                distance_cat_chain = cation - cat_n
                shortest_path = new_path

            something_found = True

    if not something_found:
        return 0, False, -1

    flat = True
    for node in shortest_path:
        neighbors = list(nx.neighbors(ligand_graph, node))
        if len(neighbors) > 2:
            verdict = is_flat_primitive(pi_acid_atoms, neighbors + [node], 0.25)
            if not verdict["isFlat"]:
                flat = False

    return len(shortest_path) + 1, flat, distance_cat_chain


def extract_ring_centroids(point, residue, ns):
    dist_atom = 5.3
    dist_cent = 5.0
    neighbors = ns.search(np.array(point), dist_atom, "A")
    residues = Selection.unfold_entities(neighbors, "R")
    point = np.array(point)

    not_aromatic = [
        "HOH",
        "DOD",
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PRO",
        "SER",
        "THR",
        "VAL",
    ]

    res1id = create_res_id(residue)
    centroids_found = []
    ring_molecules = []

    for res in residues:
        if res.get_resname() in not_aromatic:
            continue

        if create_res_id(res) == res1id:
            continue

        centroids = get_rings_centroids(res)
        for centroid in centroids:
            dist = np.linalg.norm(point - np.array(centroid["coords"]))
            if dist < dist_cent:
                centroid["distance"] = dist
                centroids_found.append(centroid)
                ring_molecules.append(res)

    return centroids_found, ring_molecules


def get_residues_list_from_atom_data(atom_data_list):
    atoms_list = []
    for atom_data in atom_data_list:
        atoms_list.append(atom_data["Atom"])
    return Selection.unfold_entities(atoms_list, "R")


def anion_screening(atoms, ligprep_data):
    selected_atoms = []
    anions_names = ligprep_data["anionNames"]

    for atom_data in atoms:
        parent_name = atom_data["Atom"].get_parent().get_resname()
        if parent_name in anions_names:
            selected_atoms.append(atom_data)

    return selected_atoms


def hydrogens_present(structure):
    for a in structure.get_atoms():
        if a.element == "H":
            return True

    return False


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    pass
