"""
Created on Sat Apr 21 13:58:24 2018

@author: michal
"""

import logging
import numpy as np
from supramolecular_search.numpy_utilities import normalize
import networkx as nx


logger = logging.getLogger(__name__)


def get_average_coords(all_atoms_list, atoms_ind_list):
    """
    Oblicz srednie wspolrzedne x, y, z wybranych atomow

    Wejscie:
    all_atoms_list - lista obiektow Atom (cala czasteczka)
    atoms_ind_list - lista indeksow atomow, ktorych usrednione polozenie
                    ma byc obliczone

    Wyjscie:
    average_coords - 3-elementowa lista zawierajaca usrednione wspolrzedne
                    x, y, z wybranych atomow
    """
    average_coords = [0.0, 0.0, 0.0]

    for atom_ind in atoms_ind_list:
        new_coords = all_atoms_list[atom_ind].get_coord()

        for coord in range(3):
            average_coords[coord] += new_coords[coord]

    atoms_no = float(len(atoms_ind_list))

    for coord in range(3):
        average_coords[coord] /= atoms_no

    return average_coords


def is_flat(all_atoms_list, atoms_ind_list, substituents):
    """
    Sprawdz czy wybrane atomy leza w jednej plaszczyznie.
    Procedura: na podstawie polozen trzech pierwszych atomow wyznacza sie
    wektor normalnych do wyznaczonej przez nich plaszczyzny. Nastepnie
    sprawdzane jest czy kolejne wiazania tworza wektory prostopadle
    do wektora normalnego. Dopuszczalne jest odchylenie 5 stopni.

    TODO: Nie lepiej byloby obliczyc tensor momentu bezwladnosci,
    zdiagonalizowac go i rozstrzygnac na podstawie jego wartosci wlasnych?

    Wejscie:
    all_atoms_list - lista obiektow Atom (cala czasteczka)
    atoms_ind_list - lista indeksow atomow, ktore maja byc zweryfikowane
                    pod wzgledem lezenia w jednej plaszczyznie

    Wyjscie:
    verdict - slownik, posiada klucze: is_flat (zmienna logiczna, True jesli
                struktura jest plaska), norm_vec (3-elementowa lista float,
                wspolrzedne wektora normalnego plaszczyzny, jesli struktura nie
                jest plaska wszystkie jego wspolrzedne sa rowne 0)
    """

    verdict = {"isFlat": False, "normVec": [0, 0, 0]}

    if len(atoms_ind_list) < 3:
        logger.warning("Found a cycle with fewer than 3 atoms: %s", atoms_ind_list)
        return verdict

    norm_vec = get_norm_vec(all_atoms_list, atoms_ind_list)

    if len(atoms_ind_list) <= 3:
        verdict["normVec"] = norm_vec
        return verdict

    centroid = get_average_coords(all_atoms_list, atoms_ind_list)
    for i in range(1, len(atoms_ind_list)):
        atom_ind = atoms_ind_list[i]
        point_d = np.array(all_atoms_list[atom_ind].get_coord())
        new_vec = normalize(centroid - point_d)

        if abs(np.inner(new_vec, norm_vec)) > 0.09:
            return verdict

    for substituent_key in substituents:
        point_d = np.array(all_atoms_list[substituent_key].get_coord())
        point_e = np.array(all_atoms_list[substituents[substituent_key]].get_coord())
        new_vec = normalize(point_e - point_d)

        if abs(np.inner(new_vec, norm_vec)) > 0.20:
            return verdict

    verdict["isFlat"] = True
    verdict["normVec"] = norm_vec
    verdict["coords"] = centroid
    return verdict


def is_flat_primitive(all_atoms_list, atoms_ind_list, max_dist=0.15):
    """
    Sprawdz czy wybrane atomy leza w jednej plaszczyznie.
    Procedura: na podstawie polozen trzech pierwszych atomow wyznacza sie
    wektor normalnych do wyznaczonej przez nich plaszczyzny. Nastepnie
    sprawdzane jest czy kolejne wiazania tworza wektory prostopadle
    do wektora normalnego. Dopuszczalne jest odchylenie 5 stopni.

    TODO: Nie lepiej byloby obliczyc tensor momentu bezwladnosci,
    zdiagonalizowac go i rozstrzygnac na podstawie jego wartosci wlasnych?

    Wejscie:
    all_atoms_list - lista obiektow Atom (cala czasteczka)
    atoms_ind_list - lista indeksow atomow, ktore maja byc zweryfikowane
                    pod wzgledem lezenia w jednej plaszczyznie

    Wyjscie:
    verdict - slownik, posiada klucze: is_flat (zmienna logiczna, True jesli
                struktura jest plaska), norm_vec (3-elementowa lista float,
                wspolrzedne wektora normalnego plaszczyzny, jesli struktura nie
                jest plaska wszystkie jego wspolrzedne sa rowne 0)
    """

    verdict = {"isFlat": False, "normVec": [0, 0, 0]}
    if len(atoms_ind_list) < 3:
        logger.warning("Found a cycle with fewer than 3 atoms: %s", atoms_ind_list)
        return verdict

    norm_vec = get_norm_vec(all_atoms_list, atoms_ind_list)

    if len(atoms_ind_list) <= 3:
        verdict["normVec"] = norm_vec
        return verdict

    centroid = get_average_coords(all_atoms_list, atoms_ind_list)
    plane_offset = -np.inner(centroid, norm_vec)

    for atom_ind in atoms_ind_list:
        atom_coord = np.array(all_atoms_list[atom_ind].get_coord())
        atom_dist = abs(np.inner(norm_vec, atom_coord) + plane_offset)
        if atom_dist > max_dist:
            return verdict

    verdict["isFlat"] = True
    verdict["normVec"] = norm_vec
    verdict["coords"] = centroid
    return verdict


def get_norm_vec(all_atoms_list, atoms_ind_list):
    norm_vec = np.array([0.0, 0.0, 0.0])
    expanded_list = atoms_ind_list + atoms_ind_list[:2]
    for i in range(len(expanded_list) - 2):
        point_a = np.array(all_atoms_list[expanded_list[i]].get_coord())
        point_b = np.array(all_atoms_list[expanded_list[i + 1]].get_coord())
        point_c = np.array(all_atoms_list[expanded_list[i + 2]].get_coord())

        vec1 = point_a - point_b
        vec2 = point_b - point_c

        norm_vec += normalize(np.cross(vec1, vec2))

    return normalize(norm_vec)


def get_ring_elements(cycle, atoms):
    element_dict = {}

    for ai in cycle:
        a = atoms[ai]
        new_element = a.element
        if new_element in element_dict:
            element_dict[new_element] += 1
        else:
            element_dict[new_element] = 1

    ring_elements = ""
    elements = sorted(list(element_dict.keys()))

    for e in elements:
        el_no = element_dict[e]
        if el_no == 1:
            ring_elements += e
        else:
            ring_elements += e + str(el_no)

    return ring_elements


def get_rings_centroids(molecule, return_graph=False):
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
        norm_vec (wektor normalnych plaszczyzny pierscienia)
    """
    atoms = list(molecule.get_atoms())
    graph = molecule2graph(atoms, None, True, True, True)

    cycles = list(nx.cycle_basis(graph))
    centroids = []

    for cycle_id, cycle in enumerate(cycles):
        if not only_ligh_atoms_in_cycle(cycle, atoms):
            continue

        if len(cycle) > 6 or len(cycle) < 5:
            continue

        substituents = get_substituents(graph, cycle)
        flat_analyse = is_flat(atoms, cycle, substituents)
        if not flat_analyse["isFlat"]:
            continue

        ring_elements = get_ring_elements(cycle, atoms)
        centroids.append(
            {
                "coords": flat_analyse["coords"],
                "normVec": flat_analyse["normVec"],
                "ringSize": len(cycle),
                "cycleAtoms": cycle,
                "cycleId": cycle_id,
                "ringElements": ring_elements,
            }
        )

    if return_graph:
        return centroids, graph
    else:
        return centroids


def only_ligh_atoms_in_cycle(cycle, atoms):
    for atom_ind in cycle:
        if atoms[atom_ind].element not in ["C", "O", "N"]:
            return False
    return True


def get_substituents(graph_molecule, cycle):
    substituents = {}

    for atom in cycle:
        candidates = graph_molecule.neighbors(atom)
        for candidate in candidates:
            if candidate not in cycle:
                substituents[atom] = candidate

    return substituents


def molecule2graph(
    atoms, atom=None, return_subgraph=True, omit_hydrogens=True, omit_metals=False
):
    """
    Konwersja czasteczki na graf (networkx)

    Wejscie:
    atom - obiekt Atom (Biopython), ktorego polozenie w grafie jest istotne
    atoms - lista wszystkich atomow

    Wyjscie:
    G, atom_ind - graf (networkx), indeks wejsciowego atomu (wierzcholek w grafie)
    """

    radius = {
        "H": 0.32,
        "D": 0.32,
        "HE": 0.46,
        "LI": 1.33,
        "BE": 1.02,
        "B": 0.85,
        "C": 0.75,
        "N": 0.71,
        "O": 0.63,
        "F": 0.64,
        "NE": 0.67,
        "NA": 1.55,
        "MG": 1.39,
        "AL": 1.26,
        "SI": 1.16,
        "P": 1.11,
        "S": 1.03,
        "CL": 0.99,
        "AR": 0.96,
        "K": 1.96,
        "CA": 1.71,
        "SC": 1.48,
        "TI": 1.36,
        "V": 1.34,
        "CR": 1.22,
        "MN": 1.19,
        "FE": 1.16,
        "CO": 1.11,
        "NI": 1.10,
        "CU": 1.12,
        "ZN": 1.18,
        "GA": 1.24,
        "GE": 1.21,
        "AS": 1.21,
        "SE": 1.16,
        "BR": 1.14,
        "KR": 1.17,
        "RB": 2.10,
        "SR": 1.85,
        "Y": 1.63,
        "ZR": 1.54,
        "NB": 1.47,
        "MO": 1.38,
        "TC": 1.28,
        "RU": 1.25,
        "RH": 1.25,
        "PD": 1.20,
        "AG": 1.28,
        "CD": 1.36,
        "IN": 1.42,
        "SN": 1.40,
        "SB": 1.40,
        "TE": 1.36,
        "I": 1.33,
        "XE": 1.31,
        "CS": 2.32,
        "BA": 1.96,
        "HF": 1.52,
        "TA": 1.46,
        "W": 1.37,
        "RE": 1.31,
        "OS": 1.29,
        "IR": 1.22,
        "PT": 1.23,
        "AU": 1.24,
        "HG": 1.33,
        "TL": 1.44,
        "PB": 1.44,
        "BI": 1.51,
        "PO": 1.45,
        "AT": 1.47,
        "RN": 1.42,
        "FR": 2.23,
        "RA": 2.01,
        "RF": 1.57,
        "DB": 1.49,
        "SG": 1.43,
        "BH": 1.41,
        "HS": 1.34,
        "MT": 1.29,
        "DS": 1.28,
        "RG": 1.21,
        "LA": 1.80,
        "CE": 1.63,
        "PR": 1.76,
        "ND": 1.74,
        "PM": 1.73,
        "SM": 1.72,
        "EU": 1.68,
        "GD": 1.69,
        "TB": 1.68,
        "DY": 1.67,
        "HO": 1.66,
        "ER": 1.65,
        "TM": 1.64,
        "YB": 1.70,
        "LU": 1.62,
        "AC": 1.86,
        "TH": 1.75,
        "PA": 1.69,
        "U": 1.70,
        "NP": 1.71,
        "PU": 1.72,
        "AM": 1.66,
        "CM": 1.66,
        "BK": 1.68,
        "CF": 1.68,
        "ES": 1.65,
        "FM": 1.67,
        "MD": 1.73,
        "NO": 1.76,
        "LR": 1.61,
    }
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

    graph = nx.Graph()
    atoms_found = []
    for atom1_ind, atom1 in enumerate(atoms):
        if atom1.element == "H" and omit_hydrogens:
            continue

        if atom1.element.upper() in metals and omit_metals:
            continue

        if atom1.element not in radius:
            continue

        graph.add_node(atom1_ind, element=atom1.element)

        radius1 = radius[atom1.element]

        if atom is not None:
            if atom == atom1:
                atoms_found.append(atom1_ind)

        for atom2_ind, atom2 in enumerate(atoms[atom1_ind + 1 :], atom1_ind + 1):
            if atom2.element == "H" and omit_hydrogens:
                continue

            if atom2.element.upper() in metals and omit_metals:
                continue

            if atom2.element not in radius:
                continue

            radius2 = radius[atom2.element]

            distance = atom1 - atom2

            threshold = 1.2 * (radius1 + radius2)
            if distance < threshold:
                graph.add_edge(atom1_ind, atom2_ind)

    if atom is not None:
        if len(atoms_found) != 1:
            logger.warning(
                "Expected exactly one matching atom in graph, found: %s", atoms_found
            )

        if atoms_found[0] not in graph.nodes():
            logger.warning("Atom %s not found in graph", atom.element)
            graph.add_node(atoms_found[0])

        if return_subgraph:
            nodes2stay = []
            for node in graph.nodes():
                if nx.has_path(graph, node, atoms_found[0]):
                    nodes2stay.append(node)

            return graph.subgraph(nodes2stay), atoms_found[0]
        else:
            return graph, atoms_found[0]

    return graph


def molecule_fragment2graph(atoms, atom, max_dist, return_subgraph=True):
    """
    Konwersja czasteczki na graf (networkx)

    Wejscie:
    atom - obiekt Atom (Biopython), ktorego polozenie w grafie jest istotne
    atoms - lista wszystkich atomow

    Wyjscie:
    G, atom_ind - graf (networkx), indeks wejsciowego atomu (wierzcholek w grafie)
    """
    radius = {
        "H": 0.32,
        "D": 0.32,
        "HE": 0.46,
        "LI": 1.33,
        "BE": 1.02,
        "B": 0.85,
        "C": 0.75,
        "N": 0.71,
        "O": 0.63,
        "F": 0.64,
        "NE": 0.67,
        "NA": 1.55,
        "MG": 1.39,
        "AL": 1.26,
        "SI": 1.16,
        "P": 1.11,
        "S": 1.03,
        "CL": 0.99,
        "AR": 0.96,
        "K": 1.96,
        "CA": 1.71,
        "SC": 1.48,
        "TI": 1.36,
        "V": 1.34,
        "CR": 1.22,
        "MN": 1.19,
        "FE": 1.16,
        "CO": 1.11,
        "NI": 1.10,
        "CU": 1.12,
        "ZN": 1.18,
        "GA": 1.24,
        "GE": 1.21,
        "AS": 1.21,
        "SE": 1.16,
        "BR": 1.14,
        "KR": 1.17,
        "RB": 2.10,
        "SR": 1.85,
        "Y": 1.63,
        "ZR": 1.54,
        "NB": 1.47,
        "MO": 1.38,
        "TC": 1.28,
        "RU": 1.25,
        "RH": 1.25,
        "PD": 1.20,
        "AG": 1.28,
        "CD": 1.36,
        "IN": 1.42,
        "SN": 1.40,
        "SB": 1.40,
        "TE": 1.36,
        "I": 1.33,
        "XE": 1.31,
        "CS": 2.32,
        "BA": 1.96,
        "HF": 1.52,
        "TA": 1.46,
        "W": 1.37,
        "RE": 1.31,
        "OS": 1.29,
        "IR": 1.22,
        "PT": 1.23,
        "AU": 1.24,
        "HG": 1.33,
        "TL": 1.44,
        "PB": 1.44,
        "BI": 1.51,
        "PO": 1.45,
        "AT": 1.47,
        "RN": 1.42,
        "FR": 2.23,
        "RA": 2.01,
        "RF": 1.57,
        "DB": 1.49,
        "SG": 1.43,
        "BH": 1.41,
        "HS": 1.34,
        "MT": 1.29,
        "DS": 1.28,
        "RG": 1.21,
        "LA": 1.80,
        "CE": 1.63,
        "PR": 1.76,
        "ND": 1.74,
        "PM": 1.73,
        "SM": 1.72,
        "EU": 1.68,
        "GD": 1.69,
        "TB": 1.68,
        "DY": 1.67,
        "HO": 1.66,
        "ER": 1.65,
        "TM": 1.64,
        "YB": 1.70,
        "LU": 1.62,
        "AC": 1.86,
        "TH": 1.75,
        "PA": 1.69,
        "U": 1.70,
        "NP": 1.71,
        "PU": 1.72,
        "AM": 1.66,
        "CM": 1.66,
        "BK": 1.68,
        "CF": 1.68,
        "ES": 1.65,
        "FM": 1.67,
        "MD": 1.73,
        "NO": 1.76,
        "LR": 1.61,
    }

    graph = nx.Graph()
    atoms_found = []
    for atom1_ind, atom1 in enumerate(atoms):
        if atom1.element == "H":
            continue

        if atom1.element not in radius:
            continue

        if atom1 - atom > max_dist:
            continue

        graph.add_node(atom1_ind, element=atom1.element)

        radius1 = radius[atom1.element]

        if atom is not None:
            if atom == atom1:
                atoms_found.append(atom1_ind)

        for atom2_ind, atom2 in enumerate(atoms[atom1_ind + 1 :], atom1_ind + 1):
            if atom2.element == "H":
                continue

            if atom2.element not in radius:
                continue

            if atom2 - atom > max_dist:
                continue

            radius2 = radius[atom2.element]

            distance = atom1 - atom2

            threshold = 1.2 * (radius1 + radius2)
            if distance < threshold:
                graph.add_edge(atom1_ind, atom2_ind)

    if atom is not None:
        if len(atoms_found) != 1:
            logger.warning(
                "Expected exactly one matching atom in graph, found: %s", atoms_found
            )

        if atoms_found[0] not in graph.nodes():
            logger.warning("Atom %s not found in graph", atom.element)
            graph.add_node(atoms_found[0])

        if return_subgraph:
            nodes2stay = []
            for node in graph.nodes():
                if nx.has_path(graph, node, atoms_found[0]):
                    nodes2stay.append(node)

            return graph.subgraph(nodes2stay), atoms_found[0]
        else:
            return graph, atoms_found[0]

    return graph


def find_in_graph(graph, atom, atom_list):
    nodes2stay = []
    atoms_found = []

    for node in graph.nodes:
        if atom_list[node] == atom:
            atoms_found.append(node)

    if len(atoms_found) != 1:
        logger.warning(
            "Expected exactly one matching atom in graph, found: %s", atoms_found
        )

    if atoms_found:
        for node in graph.nodes():
            if nx.has_path(graph, node, atoms_found[0]):
                nodes2stay.append(node)

        return graph.subgraph(nodes2stay), atoms_found[0]
    else:
        return graph.subgraph(nodes2stay), None
