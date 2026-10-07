"""
Created on Sat Apr 21 14:04:23 2018

@author: michal
"""

from supramolecular_search.config import configure

configure()

from Bio.PDB import Selection, NeighborSearch
from supramolecular_search.ring_detection import (
    get_substituents,
    is_flat,
    is_flat_primitive,
    molecule2graph,
)
from supramolecular_search.anion_template_creator import (
    AnionMatcher,
    ANION_TEMPLATES_DIR,
)
import json
from os.path import join
from glob import glob
from networkx.readwrite.json_graph import node_link_graph
from copy import copy
from supramolecular_search.biopython_utilities import (
    create_res_id,
    create_res_id_from_atom,
)


class AnionData:
    def __init__(self, anion_type, charged, anion_id, properties):
        self.anion_type = anion_type
        self.charged = charged
        self.anion_id = anion_id
        self.h_bonds_analyzed = False
        self.properties = properties


class Property:
    def __init__(
        self, prop_id, kind, atoms_involved, anion_group_id, directional_vector=[]
    ):
        self.unique_id = prop_id
        self.kind = kind
        self.atoms_involved = atoms_involved
        self.directional_vector = directional_vector
        self.anion_group_id = anion_group_id


class AnionRecogniser:
    def __init__(self):
        self.templates = get_all_templates()
        self.property_id = 0

        self.properties2calculate_pack = []

    def clean_properties_pack(self):
        self.properties2calculate_pack = []

    def extract_anion_atoms(self, atom_list_orig, ligand, ns):
        """Extract atoms that may belong to anions near the given ligand.

        Args:
            atom_list_orig: list of Biopython Atom objects to search.
            ligand: Biopython Residue whose own atoms are skipped.
            ns: NeighborSearch over the structure.

        Returns:
            List of dicts describing potential anion atoms, with keys "Atom",
            "AnionType" and "AnionId".
        """
        extracted_atoms = []
        atom_list = []

        present_properites_id = set([])

        for atom in atom_list_orig:
            if atom.get_parent() == ligand:
                continue
            if hasattr(atom, "anionData"):
                if atom.anionData.charged:
                    extracted_atoms.append(
                        {
                            "Atom": atom,
                            "AnionType": atom.anionData.anion_type,
                            "AnionId": atom.anionData.anion_id,
                        }
                    )
                    for prop in atom.anionData.properties:
                        if prop.unique_id not in present_properites_id:
                            present_properites_id.add(prop.unique_id)
                            self.properties2calculate_pack.append(prop)
            else:
                atom_list.append(atom)

        residues_neighbor = Selection.unfold_entities(atom_list, "R")
        res_id2atoms, res_id2res = create_res_dicts(
            atom_list, residues_neighbor, ligand
        )
        aminoacids = [
            "ALA",
            "ARG",
            "ASN",
            "CYS",
            "GLN",
            "GLY",
            "HIS",
            "GLU",
            "ASP",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ]

        for res_id in res_id2atoms:
            res_atoms = res_id2atoms[res_id]
            elements = atoms_list2_elements_list(res_atoms)
            res_name = res_id2res[res_id].get_resname()

            if res_name.upper() in aminoacids and not (
                "O" in elements or "S" in elements
            ):
                continue
            atoms = get_residue_with_connections(res_atoms, ns)

            ns_res = NeighborSearch(atoms)
            for atom in res_atoms:
                potentially_supramolecular = False
                anion_type = ""
                potentially_supramolecular, anion_type = self.search_in_anion_templates(
                    atom, atoms, ns_res
                )

                if potentially_supramolecular:
                    extracted_atoms.append(
                        {
                            "Atom": atom,
                            "AnionType": atom.anionData.anion_type,
                            "AnionId": atom.anionData.anion_id,
                        }
                    )

        return extracted_atoms

    def search_in_anion_templates(self, atom, atoms, ns):
        if hasattr(atom, "anionData"):
            if atom.anionData.charged:
                return True, atom.anionData.anion_type
            else:
                return False, atom.anionData.anion_type

        element = atom.element.upper()

        if element not in self.templates:
            return False, element

        atoms5 = list(ns.search(atom.get_coord(), 5.0, "A"))

        graph, atom_ind = molecule2graph(atoms5, atom, omit_metals=True)
        composition = graph2_composition(graph)

        priority2template = self.templates[element]
        for priority in sorted(priority2template.keys()):
            for guess in priority2template[priority]:
                if not dummy_compare(copy(composition), guess):
                    continue

                match_result, anion_group = self.try2match_template(
                    graph, atom_ind, guess, atoms5
                )
                if match_result:
                    return match_result, anion_group

        return False, element

    def try2match_template(self, molecule_graph, atom_id, graph_template, atoms):
        molecule_graph.nodes[atom_id]["charged"] = True
        an_matcher = AnionMatcher(molecule_graph, graph_template)

        if graph_template.graph["fullIsomorphism"]:
            result = an_matcher.is_isomorphic()
        else:
            result = an_matcher.subgraph_is_isomorphic()
        if not result:
            molecule_graph.nodes[atom_id]["charged"] = False
            return result, molecule_graph.nodes[atom_id]["element"]

        matching = an_matcher.mapping
        reverse_mapping = {}
        for key in matching:
            reverse_mapping[matching[key]] = key

        if graph_template.graph["geometry"] == "planarWithSubstituents":
            flat_analysis = is_flat(
                atoms,
                list(matching.keys()),
                get_substituents(molecule_graph, list(matching.keys())),
            )
            if not flat_analysis["isFlat"]:
                molecule_graph.nodes[atom_id]["charged"] = False
                return False, molecule_graph.nodes[atom_id]["element"]
        elif graph_template.graph["geometry"] == "planar":
            flat_analysis = is_flat_primitive(atoms, list(matching.keys()), 0.5)
            if not flat_analysis["isFlat"]:
                molecule_graph.nodes[atom_id]["charged"] = False
                return False, molecule_graph.nodes[atom_id]["element"]

        anion_group = graph_template.graph["name"]

        if "X" in graph_template.graph["name"] and graph_template.graph["nameMapping"]:
            matching = an_matcher.mapping
            for node in graph_template.graph["nameMapping"]:
                element = get_element_from_match(matching, int(node), molecule_graph)
                anion_group = anion_group.replace("X", element)
                break

        molecule_graph.nodes[atom_id]["charged"] = False

        anion_group_id = sorted(list([atoms[aid].get_name() for aid in matching]))[0]

        properties2measure = []
        if graph_template.graph["properties2measure"]:
            self.property_id += 1
            current_property_id = self.property_id

            for property2calc in graph_template.graph["properties2measure"]:
                kind = property2calc["kind"]

                atoms_involved_indexes = property2calc["atoms"]
                atoms_involved = []

                for a_ind in atoms_involved_indexes:
                    atoms_involved.append(atoms[reverse_mapping[a_ind]])

                directional_vector = []
                if kind == "plane":
                    for point_data in property2calc["directionalVector"]:
                        key = list(point_data.keys())[0]
                        if key == "atom":
                            directional_vector.append(
                                {"atom": atoms[reverse_mapping[point_data[key]]]}
                            )
                        else:
                            atoms_required = [
                                atoms[reverse_mapping[atom_ind]]
                                for atom_ind in point_data[key]
                            ]
                            directional_vector.append({key: atoms_required})

                new_property = Property(
                    current_property_id,
                    kind,
                    atoms_involved,
                    anion_group_id,
                    directional_vector,
                )
                properties2measure.append(new_property)
                self.properties2calculate_pack.append(new_property)

        for a_id in matching:
            template_atom_id = matching[a_id]

            if template_atom_id in graph_template.graph["otherCharges"]:
                atoms[a_id].anionData = AnionData(
                    anion_group, True, anion_group_id, properties2measure
                )
            else:
                atoms[a_id].anionData = AnionData(
                    anion_group, False, anion_group_id, []
                )

        atoms[atom_id].anionData = AnionData(
            anion_group, True, anion_group_id, properties2measure
        )

        return True, anion_group


def get_residue_with_connections(atoms, ns):
    atom_parent = atoms[0].get_parent()

    parent_atoms = list(atom_parent.get_atoms())
    neighbors = []
    for atom in atoms:
        neighbors += ns.search(atom.get_coord(), 2, "A")
    neighbors = list(set(neighbors))
    new_neighbors = []
    for atom in neighbors:
        new_neighbors += ns.search(atom.get_coord(), 2, "A")

    neighbors = list(set(new_neighbors))

    new_neighbors = []
    for atom in neighbors:
        if atom.element != "C":
            new_neighbors += ns.search(atom.get_coord(), 2, "A")

    neighbors = list(set(new_neighbors))
    parent = atom_parent.get_parent()
    for neighbor in neighbors:
        if parent != neighbor.get_parent():
            parent_atoms.append(neighbor)

    return list(set(parent_atoms))


def atoms_list2_elements_list(atom_list):
    elements = set()
    for a in atom_list:
        elements.add(a.element)

    return elements


def create_res_dicts(atoms, residues, ligand):
    ligand_id = create_res_id(ligand)
    res_id2res = {}
    for res in residues:
        res_id2res[create_res_id(res)] = res

    bad_parents = ["HOH", "DOD", "OXY"]
    res_id2atoms = {}

    for a in atoms:
        res_id = create_res_id_from_atom(a)
        if res_id == ligand_id:
            continue
        resname = a.get_parent().get_resname()

        if resname in bad_parents:
            continue
        if res_id not in res_id2atoms:
            res_id2atoms[res_id] = [a]
        else:
            res_id2atoms[res_id].append(a)

    return res_id2atoms, res_id2res


def graph2_composition(graph):
    composition = {}
    for node in list(graph):
        element = graph.nodes[node]["element"]
        if element not in composition:
            composition[element] = 1
        else:
            composition[element] += 1

    return composition


def dummy_compare(composition, template):
    for node in list(template):
        if "aliases" not in template.nodes[node]:
            element = template.nodes[node]["element"]
            if element == "X":
                continue

            if element not in composition:
                return False
            else:
                composition[element] -= 1
                if composition[element] < 0:
                    return False
        elif not template.nodes[node]["aliases"]:
            element = template.nodes[node]["element"]
            if element == "X":
                continue
            if element not in composition:
                return False
            else:
                composition[element] -= 1
                if composition[element] < 0:
                    return False
    return True


def get_all_templates():
    templates = join(ANION_TEMPLATES_DIR, "*", "*.json")
    anions_templates = glob(templates)

    graph_templates = {}

    for template in anions_templates:
        json_f = open(template)
        graph_template = node_link_graph(json.load(json_f), edges="links")
        json_f.close()

        priority = graph_template.graph["priority"]
        element = graph_template.nodes[graph_template.graph["charged"]]["element"]

        if element not in graph_templates:
            graph_templates[element] = {priority: [graph_template]}
        elif priority not in graph_templates[element]:
            graph_templates[element][priority] = [graph_template]
        else:
            graph_templates[element][priority].append(graph_template)

    return graph_templates


def get_element_from_match(matching, node, graph):
    for match in matching:
        if matching[match] == node:
            return graph.nodes[match]["element"]


if __name__ == "__main__":
    pass
