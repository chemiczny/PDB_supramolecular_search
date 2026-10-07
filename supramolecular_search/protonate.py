"""
Created on Fri Jan 18 15:16:17 2019

@author: michal
"""

from supramolecular_search.ring_detection import molecule2graph
import networkx as nx
import math
import numpy as np
from supramolecular_search.numpy_utilities import normalize, rotate_vector


class HydrogenAtom(object):
    def __init__(self, coord):
        self.coord = coord
        self.element = "H"

    def get_coord(self):
        return self.coord

    def __sub__(self, atom):
        r2 = 0
        for x1, x2 in zip(self.coord, atom.coord):
            r2 += (x1 - x2) * (x1 - x2)

        return math.sqrt(r2)


class Protonate:
    """Protonates atoms using VSEPR theory"""

    def __init__(self, verbose=False):
        self.verbose = verbose

        self.valence_electrons = {
            "N": 5,
            "O": 6,
        }

        self.standard_charges = {
            "ARG-NH1": 1.0,
            "ASP-OD2": -1.0,
            "GLU-OE2": -1.0,
            "HIS-ND1": 1.0,
            "LYS-NZ": 1.0,
            "ARG-NE": 1,
            "N+": 1.0,
            "C-": -1.0,
        }

        self.sybyl_charges = {
            "N.pl3": +1,
            "N.3": +1,
            "N.4": +1,
            "N.ar": +1,
            "O.co2-": -1,
        }

        self.bond_lengths = {
            "C": 1.09,
            "N": 1.01,
            "O": 0.96,
            "F": 0.92,
            "Cl": 1.27,
            "Br": 1.41,
            "I": 1.61,
            "S": 1.35,
        }

        self.number_of_pi_electrons_in_bonds_in_backbone = {"O": 1}

        self.number_of_pi_electrons_in_conjugate_bonds_in_backbone = {"N": 1}

        self.number_of_pi_electrons_in_bonds_in_sidechains = {
            "ARG-NH1": 1,
            "ASN-OD1": 1,
            "ASP-OD1": 1,
            "GLU-OE1": 1,
            "GLN-OE1": 1,
            "HIS-ND1": 1,
        }

        self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains = {
            "ARG-NH2": 1,
            "ASN-ND2": 1,
            "GLN-NE2": 1,
            "HIS-NE2": 1,
            "TRP-NE1": 1,
            "GLU-OE2": 1,
            "ASP-OD2": 1,
            "ARG-NE": 1,
        }

        self.number_of_pi_electrons_in_bonds_ligands = {
            "N.pl3": 0,
            "O.2": 1,
            "O.co2": 1,
            "N.ar": 1,
            "N.1": 2,
        }

        self.number_of_pi_electrons_in_conjugate_bonds_in_ligands = {
            "N.am": 1,
            "N.pl3": 1,
        }

        self.protonation_methods = {4: self.tetrahedral, 3: self.trigonal}

        self.molecule_graph = {}

        self.anion_id = -1
        self.atom_list = []
        self.hydrogen_atoms_list = []
        self.connected2_h = []

        self.aa_codes = set(
            [
                "ALA",
                "ARG",
                "ASN",
                "ASP",
                "CYS",
                "GLN",
                "GLU",
                "GLY",
                "HIS",
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
        )

        return

    def protonate(self, atom_list, anion_atom):
        self.atom_list = atom_list
        self.molecule_graph, self.anion_id = molecule2graph(
            atom_list, anion_atom, False, omit_metals=True
        )
        self.molecule_graph = self.molecule_graph.copy()

        self.anion_coords = anion_atom.get_coord()
        # protonate all atoms
        for atom_id in self.molecule_graph.nodes():
            atom = self.atom_list[atom_id]
            element = atom.element
            res_name = atom.get_parent().get_resname()

            if (
                element in ["O", "N"]
                and atom != anion_atom
                and res_name in self.aa_codes
            ):
                self.protonate_atom(atom_id)

        h_ind = len(self.atom_list)
        for conn_ind in self.connected2_h:
            self.molecule_graph.add_edge(h_ind, conn_ind)
            h_ind += 1
        self.atom_list += self.hydrogen_atoms_list

        return

    def set_charge(self, atom_id):
        atom = self.atom_list[atom_id]
        # atom is a protein atom

        key = "%3s-%s" % (atom.get_parent().get_resname(), atom.get_name())
        self.molecule_graph.nodes[atom_id]["key"] = key
        if key in list(self.standard_charges.keys()):
            self.molecule_graph.nodes[atom_id]["charge"] = self.standard_charges[key]
        else:
            self.molecule_graph.nodes[atom_id]["charge"] = 0

        return

    def protonate_atom(self, atom):

        self.set_charge(atom)
        self.set_number_of_pi_electrons(atom)
        self.set_number_of_protons_to_add(atom)
        self.set_steric_number_and_lone_pairs(atom)
        self.add_protons(atom)
        return

    def set_number_of_pi_electrons(self, atom):
        atom_key = self.molecule_graph.nodes[atom]["key"]

        atom_obj = self.atom_list[atom]
        aminoacid = False
        if atom_obj.get_parent().get_resname() in [
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
            "HIS",
            "TRP",
            "PHE",
            "TYR",
        ]:
            aminoacid = True

        if aminoacid:
            if atom_obj.get_name() in ["O", "OXT"]:
                self.molecule_graph.nodes[atom]["number_of_pi_electrons_in_bonds"] = 1
            elif atom_key in self.number_of_pi_electrons_in_bonds_in_sidechains:
                self.molecule_graph.nodes[atom]["number_of_pi_electrons_in_bonds"] = (
                    self.number_of_pi_electrons_in_bonds_in_sidechains[atom_key]
                )
            else:
                self.molecule_graph.nodes[atom]["number_of_pi_electrons_in_bonds"] = 0

            if atom_obj.get_name() == "N":
                self.molecule_graph.nodes[atom][
                    "number_of_pi_electrons_in_conjugate_bonds"
                ] = 1
            elif (
                atom_key in self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains
            ):
                self.molecule_graph.nodes[atom][
                    "number_of_pi_electrons_in_conjugate_bonds"
                ] = self.number_of_pi_electrons_in_conjugate_bonds_in_sidechains[
                    atom_key
                ]
            else:
                self.molecule_graph.nodes[atom][
                    "number_of_pi_electrons_in_conjugate_bonds"
                ] = 0

        else:
            if atom_obj.get_name() in self.number_of_pi_electrons_in_bonds_ligands:
                self.molecule_graph.nodes[atom]["number_of_pi_electrons_in_bonds"] = (
                    self.number_of_pi_electrons_in_bonds_ligands[atom_obj.get_name()]
                )
            else:
                self.molecule_graph.nodes[atom]["number_of_pi_electrons_in_bonds"] = 0

            if (
                atom_obj.get_name()
                in self.number_of_pi_electrons_in_conjugate_bonds_in_ligands
            ):
                self.molecule_graph.nodes[atom][
                    "number_of_pi_electrons_in_conjugate_bonds"
                ] = self.number_of_pi_electrons_in_conjugate_bonds_in_ligands[
                    atom_obj.get_name()
                ]
            else:
                self.molecule_graph.nodes[atom][
                    "number_of_pi_electrons_in_conjugate_bonds"
                ] = 0

    def set_number_of_protons_to_add(self, atom):
        number_of_protons_to_add = 8
        number_of_protons_to_add -= self.valence_electrons[self.atom_list[atom].element]
        number_of_protons_to_add -= len(list(nx.neighbors(self.molecule_graph, atom)))

        number_of_pi_electrons_in_double_and_triple_bonds = self.molecule_graph.nodes[
            atom
        ]["number_of_pi_electrons_in_bonds"]
        number_of_protons_to_add -= number_of_pi_electrons_in_double_and_triple_bonds
        number_of_protons_to_add += int(self.molecule_graph.nodes[atom]["charge"])

        self.molecule_graph.nodes[atom]["number_of_protons_to_add"] = (
            number_of_protons_to_add
        )

    def set_steric_number_and_lone_pairs(self, atom):
        steric_number = 0

        steric_number += self.valence_electrons[self.atom_list[atom].element]
        steric_number += len(list(nx.neighbors(self.molecule_graph, atom)))
        steric_number += self.molecule_graph.nodes[atom]["number_of_protons_to_add"]
        steric_number -= self.molecule_graph.nodes[atom]["charge"]
        steric_number -= self.molecule_graph.nodes[atom][
            "number_of_pi_electrons_in_bonds"
        ]
        steric_number -= self.molecule_graph.nodes[atom][
            "number_of_pi_electrons_in_conjugate_bonds"
        ]

        self.molecule_graph.nodes[atom]["steric_number"] = math.floor(
            steric_number / 2.0
        )

        self.molecule_graph.nodes[atom]["number_of_lone_pairs"] = (
            steric_number
            - len(list(nx.neighbors(self.molecule_graph, atom)))
            - self.molecule_graph.nodes[atom]["number_of_protons_to_add"]
        )

        self.molecule_graph.nodes[atom]["steric_number_and_lone_pairs_set"] = True
        return

    def add_protons(self, atom):
        # decide which method to use

        if self.molecule_graph.nodes[atom]["steric_number"] in list(
            self.protonation_methods.keys()
        ):
            self.protonation_methods[self.molecule_graph.nodes[atom]["steric_number"]](
                atom
            )

        return

    def trigonal(self, atom):
        number_of_protons_to_add = self.molecule_graph.nodes[atom][
            "number_of_protons_to_add"
        ]

        if number_of_protons_to_add == 0:
            return

        rot_angle = math.radians(120.0)
        bonded_atoms_ids = list(self.molecule_graph.neighbors(atom))

        if len(bonded_atoms_ids) == 0:
            return

        if len(bonded_atoms_ids) == 1:
            point_a = self.atom_list[atom].get_coord()
            point_b = self.atom_list[bonded_atoms_ids[0]].get_coord()

            b_neighbors = list(self.molecule_graph.neighbors(bonded_atoms_ids[0]))
            for c_candidate in b_neighbors:
                if c_candidate != atom and self.atom_list[c_candidate].element in [
                    "N",
                    "C",
                ]:
                    point_c = self.atom_list[c_candidate].get_coord()
                    norm_vec = -normalize(
                        np.cross(point_a - point_b, point_b - point_c)
                    )
                    break
            else:
                norm_vec = normalize(
                    np.cross(point_b - point_a, self.anion_coords - point_a)
                )

            bond_direction = point_b - point_a
            for i in range(number_of_protons_to_add):
                bond_direction = rotate_vector(bond_direction, norm_vec, rot_angle)
                bond_direction = (
                    normalize(bond_direction)
                    * self.bond_lengths[self.atom_list[atom].element]
                )

                new_atom_coords = point_a + bond_direction
                self.hydrogen_atoms_list.append(HydrogenAtom(new_atom_coords))
                self.connected2_h.append(atom)

        elif len(bonded_atoms_ids) == 2:
            point_a = self.atom_list[atom].get_coord()
            point_b = self.atom_list[bonded_atoms_ids[0]].get_coord()
            point_c = self.atom_list[bonded_atoms_ids[1]].get_coord()

            vec_ab = normalize(point_b - point_a)
            vec_ac = normalize(point_c - point_a)

            bond_direction = -(vec_ab + vec_ac)
            bond_direction = (
                normalize(bond_direction)
                * self.bond_lengths[self.atom_list[atom].element]
            )

            new_atom_coords = point_a + bond_direction

            self.hydrogen_atoms_list.append(HydrogenAtom(new_atom_coords))
            self.connected2_h.append(atom)

        return

    def tetrahedral(self, atom):
        number_of_protons_to_add = self.molecule_graph.nodes[atom][
            "number_of_protons_to_add"
        ]
        rot_angle = math.radians(109.5)

        if number_of_protons_to_add == 0:
            return

        bonded_atoms_ids = list(self.molecule_graph.neighbors(atom))

        if len(bonded_atoms_ids) == 0:
            return

        if len(bonded_atoms_ids) == 1:
            point_a = self.atom_list[atom].get_coord()
            point_b = self.atom_list[bonded_atoms_ids[0]].get_coord()

            norm_vec = normalize(
                np.cross(point_b - point_a, self.anion_coords - point_a)
            )
            dih_rot = math.radians(120)
            bond_direction = rotate_vector(point_b - point_a, norm_vec, rot_angle)
            bond_direction = (
                normalize(bond_direction)
                * self.bond_lengths[self.atom_list[atom].element]
            )

            for i in range(number_of_protons_to_add):
                new_atom_coords = point_a + bond_direction
                self.hydrogen_atoms_list.append(HydrogenAtom(new_atom_coords))
                self.connected2_h.append(atom)

                bond_direction = rotate_vector(
                    bond_direction, point_b - point_a, dih_rot
                )

        # 1 bond

        elif len(bonded_atoms_ids) == 2:
            point_a = self.atom_list[atom].get_coord()
            point_b = self.atom_list[bonded_atoms_ids[0]].get_coord()
            point_c = self.atom_list[bonded_atoms_ids[1]].get_coord()

            vec_ab = normalize(point_b - point_a)
            vec_ac = normalize(point_c - point_a)

            axis = vec_ab + vec_ac
            bond_direction = rotate_vector(-vec_ab, axis, math.radians(90))
            bond_direction = (
                normalize(bond_direction)
                * self.bond_lengths[self.atom_list[atom].element]
            )

            new_atom_coords = point_a + bond_direction
            self.hydrogen_atoms_list.append(HydrogenAtom(new_atom_coords))
            self.connected2_h.append(atom)

            if number_of_protons_to_add > 1:
                bond_direction = rotate_vector(bond_direction, axis, math.radians(180))
                new_atom_coords = point_a + bond_direction
                self.hydrogen_atoms_list.append(HydrogenAtom(new_atom_coords))
                self.connected2_h.append(atom)

        elif len(bonded_atoms_ids) == 3:
            point_a = self.atom_list[atom].get_coord()
            point_b = self.atom_list[bonded_atoms_ids[0]].get_coord()
            point_c = self.atom_list[bonded_atoms_ids[1]].get_coord()
            point_d = self.atom_list[bonded_atoms_ids[2]].get_coord()

            vec_ab = normalize(point_b - point_a)
            vec_ac = normalize(point_c - point_a)
            vec_ad = normalize(point_d - point_a)

            bond_direction = -(vec_ab + vec_ac + vec_ad)
            bond_direction = (
                bond_direction * self.bond_lengths[self.atom_list[atom].element]
            )

            new_atom_coords = point_a + bond_direction
            self.hydrogen_atoms_list.append(HydrogenAtom(new_atom_coords))
            self.connected2_h.append(atom)

        return


def get_ortonormal(vec):
    closest2zero_index = -1
    dist = 100
    for i, el in enumerate(vec):
        if abs(el) < dist:
            closest2zero_index = i
            dist = abs(el)

    vec_out = np.array([0.0, 0.0, 0.0])

    a_ind = (closest2zero_index + 1) % 3
    b_ind = (closest2zero_index + 2) % 3

    if abs(vec[b_ind]) < 0.001:
        vec_out[b_ind] = 1
        return vec_out

    vec_out[a_ind] = 1
    vec_out[b_ind] = -vec[a_ind] / vec[b_ind]

    return normalize(vec_out)
