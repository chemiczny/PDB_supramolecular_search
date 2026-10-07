"""
Created on Sat Apr 21 13:44:59 2018

@author: michal
"""

from math import sin, cos, radians, sqrt, acos, degrees
import numpy as np
from supramolecular_search.numpy_utilities import normalize
from os.path import isfile, join
from supramolecular_search.ring_detection import get_norm_vec, get_average_coords


def write_anion_pi_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/anionPi.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    results_file.write(
        "Anion code\tAnion chain\tAnion id\tAnion type\tAnion group id\t"
    )
    results_file.write("Atom symbol\tDistance\tAngle\t")
    results_file.write("x\th\t")
    results_file.write("CentroidId\t")
    results_file.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    results_file.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    results_file.write("Model No\tDisordered\t")
    results_file.write("Ring size\tRing elements\t")
    results_file.write("Resolution\t")
    results_file.write("Method\tStructure type\n")
    results_file.close()


def write_cation_pi_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/cationPi.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    results_file.write("Cation code\tCation chain\tCation id\t")
    results_file.write("Atom symbol\tDistance\tAngle\t")
    results_file.write("x\th\t")
    results_file.write("RingChain\t")
    results_file.write("ChainFlat\t")
    results_file.write("Cation-Chain Distance\t")
    results_file.write("CentroidId\t")
    results_file.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    results_file.write("Cation x coord\tCation y coord\tCation z coord\t")
    results_file.write("Model No\n")
    results_file.close()


def write_metal_ligand_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/metalLigand.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\t")
    results_file.write("Cation code\tCation chain\tCation id\t")
    results_file.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    results_file.write("Cation element\t")
    results_file.write("Ligand element\t")
    results_file.write("isAnion\tanionType\tDistance\t")
    results_file.write("Cation x coord\tCation y coord\tCation z coord\t")
    results_file.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    results_file.write("Complex\tSummary\tCoordNo\t")
    results_file.write("Model No\n")
    results_file.close()


def write_pi_pi_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/piPi.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    results_file.write("Pi res code\tPi res chain\tPi res id\t")
    results_file.write("Distance\tAngle\t")
    results_file.write("x\th\t")
    results_file.write("theta\t")
    results_file.write("omega\t")
    results_file.write("CentroidId\t")
    results_file.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    results_file.write("Centroid 2 x coord\tCentroid 2 y coord\tCentroid 2 z coord\t")
    results_file.write("Model No\t")
    results_file.write("Ring size 2\n")
    results_file.close()


def write_anion_cation_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/anionCation.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tCation code\tCation chain\tCation id\t")
    results_file.write("Pi acid Code\tPi acid chain\tPiacid id\t")
    results_file.write("CentroidId\t")
    results_file.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    results_file.write("Anion symbol\tCation symbol\tDistance\t")
    results_file.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    results_file.write("Cation x coord\tCation y coord\tCation z coord\t")
    results_file.write("Same semisphere\t")
    results_file.write("Latitude diff\t")
    results_file.write("Model No\n")
    results_file.close()


def write_hbonds_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/hBonds.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tAnion code\tAnion chain\tAnion id\t")
    results_file.write("Donor code\tDonor chain\tDonor id\t")
    results_file.write("Acceptor group\tAcceptor atom\tAnion group id\t")
    results_file.write("Acceptor x coord\tAcceptor y coord\tAcceptor z coord\t")
    results_file.write("Donor group\tDonor atom\t")
    results_file.write("Donor x coord\tDonor y coord\tDonor z coord\t")
    results_file.write("Hydrogen x coord\tHydrogen y coord\tHydrogen z coord\t")
    results_file.write("H from Experm\tAngle\tDistance H Acc\t")
    results_file.write("Distance Don Acc\tModel No\n")
    results_file.close()


def write_anion_pi_planar_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/planarAnionPi.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    results_file.write("CentroidId\t")
    results_file.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    results_file.write("Angle\t")
    results_file.write("DirectionalAngle\t")
    results_file.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    results_file.write(
        "Anion group x coord\tAnion group y coord\tAnion group z coord\t"
    )
    results_file.write("Model No\n")
    results_file.close()


def write_anion_pi_linear_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/linearAnionPi.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    results_file.write("CentroidId\t")
    results_file.write("Anion code\tAnion chain\tAnion id\tAnion group id\t")
    results_file.write("Angle\t")
    results_file.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    results_file.write(
        "Anion group x coord\tAnion group y coord\tAnion group z coord\t"
    )
    results_file.write("Model No\n")
    results_file.close()


def write_methyl_pi_header():
    """
    Write the column headers of the results file.
    """
    results_file_name = "logs/methylPi.log"
    results_file = open(results_file_name, "w")
    results_file.write("PDB Code\tPi acid Code\tPi acid chain\tPiacid id\t")
    results_file.write(
        "Anion code\tAnion chain\tAnion id\tAnion type\tAnion group id\t"
    )
    results_file.write("Atom symbol\tDistance\tAngle\t")
    results_file.write("x\th\t")
    results_file.write("CentroidId\t")
    results_file.write("Centroid x coord\tCentroid y coord\tCentroid z coord\t")
    results_file.write("Anion x coord\tAnion y coord\tAnion z coord\t")
    results_file.write("Model No\tDisordered\t")
    results_file.write("Ring size\tRing elements\t")
    results_file.write("Resolution\t")
    results_file.write("Method\tStructure type\n")
    results_file.close()


class SupramolecularLogger:
    def __init__(self, pdb_code, file_id=None, scratch_dir=None):
        self.pdb_code = pdb_code
        self.file_id = file_id
        self.scratch_dir = scratch_dir
        self.final_log_dir = "logs"

        if file_id:
            self.additional_info_log = join(
                self.scratch_dir, "additionalInfo{}.log".format(self.file_id)
            )
            self.anion_pi_log = join(
                self.scratch_dir, "anionPi{}.log".format(self.file_id)
            )
            self.planar_anion_pi_log = join(
                self.scratch_dir, "planarAnionPi{}.log".format(self.file_id)
            )
            self.linear_anion_pi_log = join(
                self.scratch_dir, "linearAnionPi{}.log".format(self.file_id)
            )
            self.cation_pi_log = join(
                self.scratch_dir, "cationPi{}.log".format(self.file_id)
            )
            self.metal_ligand_log = join(
                self.scratch_dir, "metalLigand{}.log".format(self.file_id)
            )
            self.pi_pi_log = join(self.scratch_dir, "piPi{}.log".format(self.file_id))
            self.anion_cation_log = join(
                self.scratch_dir, "anionCation{}.log".format(self.file_id)
            )
            self.h_bonds_log = join(
                self.scratch_dir, "hBonds{}.log".format(self.file_id)
            )
            self.methyl_pi_log = join(
                self.scratch_dir, "methylPi{}.log".format(self.file_id)
            )

        else:
            self.additional_info_log = join(self.final_log_dir, "additionalInfo.log")
            self.anion_pi_log = join(self.final_log_dir, "anionPi.log")
            self.planar_anion_pi_log = join(self.final_log_dir, "planarAnionPi.log")
            self.linear_anion_pi_log = join(self.final_log_dir, "linearAnionPi.log")
            self.cation_pi_log = join(self.final_log_dir, "cationPi.log")
            self.metal_ligand_log = join(self.final_log_dir, "metalLigand.log")
            self.pi_pi_log = join(self.final_log_dir, "piPi.log")
            self.anion_cation_log = join(self.final_log_dir, "anionCation.log")
            self.h_bonds_log = join(self.final_log_dir, "hBonds.log")
            self.methyl_log = join(self.scratch_dir, "methylPi.log")

        self.partial_progress_log = join(
            self.scratch_dir, "partialProgress{}.log".format(self.file_id)
        )

    def write_additional_info(self, message):
        results_file_name = self.additional_info_log

        results = open(results_file_name, "a+")
        results.write(message + "\n")
        results.close()

    def write_anion_pi_results(
        self,
        ligand,
        centroid,
        extracted_atoms,
        model_index,
        resolution,
        method,
        structure_type,
    ):
        """
        Write the results to the results file.
        """
        results_file_name = self.anion_pi_log
        ligand_code = ligand.get_resname()
        ligand_id = str(ligand.get_id()[1])
        ligand_chain = ligand.get_parent().get_id()
        results_file = open(results_file_name, "a+")
        new_atoms = []
        for atom_data in extracted_atoms:
            distance = atom_distance_from_centroid(atom_data["Atom"], centroid)
            angle = atom_angle_nom_vec_centroid(atom_data["Atom"], centroid)

            h = abs(cos(radians(angle)) * distance)
            x = sin(radians(angle)) * distance
            if angle > 90.0:
                angle = 180 - angle
            new_atoms.append(atom_data)

            atom_coords = atom_data["Atom"].get_coord()
            centroid_coords = centroid["coords"]

            anion = atom_data["Atom"].get_parent()
            residue_name = anion.get_resname()
            anion_chain = anion.get_parent().get_id()
            anion_id = str(anion.get_id()[1])
            results_file.write(self.pdb_code + "\t")
            results_file.write(ligand_code + "\t")
            results_file.write(ligand_chain + "\t")
            results_file.write(ligand_id + "\t")
            results_file.write(residue_name + "\t")
            results_file.write(anion_chain + "\t")
            results_file.write(anion_id + "\t")
            results_file.write(atom_data["AnionType"] + "\t")
            results_file.write(str(atom_data["AnionId"]) + "\t")
            results_file.write(atom_data["Atom"].element + "\t")

            results_file.write(str(distance) + "\t")
            results_file.write(str(angle) + "\t")

            results_file.write(str(x) + "\t")
            results_file.write(str(h) + "\t")

            results_file.write(str(centroid["cycleId"]) + "\t")
            results_file.write(str(centroid_coords[0]) + "\t")
            results_file.write(str(centroid_coords[1]) + "\t")
            results_file.write(str(centroid_coords[2]) + "\t")

            results_file.write(str(atom_coords[0]) + "\t")
            results_file.write(str(atom_coords[1]) + "\t")
            results_file.write(str(atom_coords[2]) + "\t")

            results_file.write(str(model_index) + "\t")
            results_file.write(
                str(atom_data["Atom"].get_parent().is_disordered()) + "\t"
            )
            results_file.write(str(centroid["ringSize"]) + "\t")
            results_file.write(str(centroid["ringElements"]) + "\t")
            results_file.write(str(resolution) + "\t")
            results_file.write(str(method) + "\t")
            results_file.write(str(structure_type) + "\n")

        results_file.close()

        return new_atoms

    def write_anion_pi_planar_results(
        self, ligand, centroid, plane_data, model_index, anion_group_id
    ):
        """
        Write the column headers of the results file.
        """
        results_file_name = self.planar_anion_pi_log

        ligand_code = ligand.get_resname()
        ligand_id = str(ligand.get_id()[1])
        ligand_chain = ligand.get_parent().get_id()
        results_file = open(results_file_name, "a+")

        norm_vec = get_norm_vec(
            plane_data.atoms_involved, list(range(len(plane_data.atoms_involved)))
        )

        inner_prod = np.inner(norm_vec, centroid["normVec"])
        if abs(inner_prod) > 1.0:
            if abs(inner_prod) < 1.1:
                angle = 0
            else:
                angle = 666
        else:
            angle = degrees(acos(inner_prod))

        if angle > 90.0:
            angle = 180 - angle

        atom = plane_data.atoms_involved[0]
        anion = atom.get_parent()
        residue_name = anion.get_resname()
        anion_chain = anion.get_parent().get_id()
        anion_id = str(anion.get_id()[1])
        centroid_coords = centroid["coords"]
        anion_group_coords = get_average_coords(
            plane_data.atoms_involved, list(range(len(plane_data.atoms_involved)))
        )

        directional_angle = calc_directional_vector(plane_data, centroid)

        results_file.write(self.pdb_code + "\t")
        results_file.write(ligand_code + "\t")
        results_file.write(ligand_chain + "\t")
        results_file.write(ligand_id + "\t")
        results_file.write(str(centroid["cycleId"]) + "\t")
        results_file.write(residue_name + "\t")
        results_file.write(anion_chain + "\t")
        results_file.write(anion_id + "\t")
        results_file.write(str(anion_group_id) + "\t")
        results_file.write(str(angle) + "\t")
        results_file.write(str(directional_angle) + "\t")

        results_file.write(str(centroid_coords[0]) + "\t")
        results_file.write(str(centroid_coords[1]) + "\t")
        results_file.write(str(centroid_coords[2]) + "\t")

        results_file.write(str(anion_group_coords[0]) + "\t")
        results_file.write(str(anion_group_coords[1]) + "\t")
        results_file.write(str(anion_group_coords[2]) + "\t")

        results_file.write(str(model_index) + "\n")
        results_file.close()

    def write_anion_pi_linear_results(
        self,
        ligand,
        centroid,
        line_data,
        model_index,
        anion_group_id,
        symmetrize_alpha=False,
    ):
        results_file_name = self.linear_anion_pi_log

        ligand_code = ligand.get_resname()
        ligand_id = str(ligand.get_id()[1])
        ligand_chain = ligand.get_parent().get_id()
        results_file = open(results_file_name, "a+")

        vector = (
            line_data.atoms_involved[1].get_coord()
            - line_data.atoms_involved[0].get_coord()
        )
        vector = normalize(vector)

        inner_prod = np.inner(vector, centroid["normVec"])
        if abs(inner_prod) > 1.0:
            if abs(inner_prod) < 1.1:
                angle = 0
            else:
                angle = 666
        else:
            angle = degrees(acos(inner_prod))

        if symmetrize_alpha and angle > 90.0:
            angle = 180 - angle

        atom = line_data.atoms_involved[0]
        anion = atom.get_parent()
        residue_name = anion.get_resname()
        anion_chain = anion.get_parent().get_id()
        anion_id = str(anion.get_id()[1])
        centroid_coords = centroid["coords"]
        anion_group_coords = get_average_coords(
            line_data.atoms_involved, list(range(len(line_data.atoms_involved)))
        )

        results_file.write(self.pdb_code + "\t")
        results_file.write(ligand_code + "\t")
        results_file.write(ligand_chain + "\t")
        results_file.write(ligand_id + "\t")
        results_file.write(str(centroid["cycleId"]) + "\t")
        results_file.write(residue_name + "\t")
        results_file.write(anion_chain + "\t")
        results_file.write(anion_id + "\t")
        results_file.write(str(anion_group_id) + "\t")
        results_file.write(str(angle) + "\t")

        results_file.write(str(centroid_coords[0]) + "\t")
        results_file.write(str(centroid_coords[1]) + "\t")
        results_file.write(str(centroid_coords[2]) + "\t")

        results_file.write(str(anion_group_coords[0]) + "\t")
        results_file.write(str(anion_group_coords[1]) + "\t")
        results_file.write(str(anion_group_coords[2]) + "\t")

        results_file.write(str(model_index) + "\n")
        results_file.close()

    def write_cation_pi_results(
        self, ligand, centroid, extracted_atoms, cation_ring_chain_lens, model_index
    ):
        """
        Write the results to the results file.
        """
        if not extracted_atoms:
            return

        results_file_name = self.cation_pi_log
        ligand_code = ligand.get_resname()
        ligand_id = str(ligand.get_id()[1])
        ligand_chain = ligand.get_parent().get_id()
        results_file = open(results_file_name, "a+")

        for atom, chain_len in zip(extracted_atoms, cation_ring_chain_lens):
            distance = atom_distance_from_centroid(atom, centroid)
            angle = atom_angle_nom_vec_centroid(atom, centroid)

            h = abs(cos(radians(angle)) * distance)
            x = sin(radians(angle)) * distance
            if angle > 90.0:
                angle = 180 - angle

            atom_coords = atom.get_coord()
            centroid_coords = centroid["coords"]

            cation = atom.get_parent()
            residue_name = cation.get_resname()
            cation_chain = cation.get_parent().get_id()
            cation_id = str(cation.get_id()[1])
            results_file.write(self.pdb_code + "\t")
            results_file.write(ligand_code + "\t")
            results_file.write(ligand_chain + "\t")
            results_file.write(ligand_id + "\t")
            results_file.write(residue_name + "\t")
            results_file.write(cation_chain + "\t")
            results_file.write(cation_id + "\t")
            results_file.write(atom.element + "\t")

            results_file.write(str(distance) + "\t")
            results_file.write(str(angle) + "\t")

            results_file.write(str(x) + "\t")
            results_file.write(str(h) + "\t")

            results_file.write(str(chain_len[0]) + "\t")
            results_file.write(str(chain_len[1]) + "\t")
            results_file.write(str(chain_len[2]) + "\t")

            results_file.write(str(centroid["cycleId"]) + "\t")
            results_file.write(str(centroid_coords[0]) + "\t")
            results_file.write(str(centroid_coords[1]) + "\t")
            results_file.write(str(centroid_coords[2]) + "\t")

            results_file.write(str(atom_coords[0]) + "\t")
            results_file.write(str(atom_coords[1]) + "\t")
            results_file.write(str(atom_coords[2]) + "\t")

            results_file.write(str(model_index) + "\n")

        results_file.close()

    def write_metal_ligand_results(self, extracted_atoms, complex_data, model_index):
        """
        Write the results to the results file.
        """
        if not extracted_atoms:
            return

        results_file_name = self.metal_ligand_log
        results_file = open(results_file_name, "a+")

        for atom, comp_data in zip(extracted_atoms, complex_data):
            if not comp_data["ligands"]:
                continue

            cation = atom.get_parent()
            cation_residue_name = cation.get_resname()
            cation_chain = cation.get_parent().get_id()
            cation_id = str(cation.get_id()[1])
            cation_coords = atom.get_coord()

            for ligand_data in comp_data["ligands"]:
                distance = atom - ligand_data["atom"]

                ligand = ligand_data["atom"].get_parent()
                ligand_residue_name = ligand.get_resname()
                ligand_chain = ligand.get_parent().get_id()
                ligand_id = str(ligand.get_id()[1])

                results_file.write(self.pdb_code + "\t")

                results_file.write(cation_residue_name + "\t")
                results_file.write(cation_chain + "\t")
                results_file.write(cation_id + "\t")

                results_file.write(ligand_residue_name + "\t")
                results_file.write(ligand_chain + "\t")
                results_file.write(ligand_id + "\t")
                results_file.write(str(ligand_data["AnionId"]) + "\t")

                results_file.write(atom.element + "\t")
                results_file.write(ligand_data["atom"].element + "\t")

                results_file.write(str(ligand_data["isAnion"]) + "\t")
                results_file.write(ligand_data["anionType"] + "\t")

                results_file.write(str(distance) + "\t")

                ligand_coords = ligand_data["atom"].get_coord()

                results_file.write(str(cation_coords[0]) + "\t")
                results_file.write(str(cation_coords[1]) + "\t")
                results_file.write(str(cation_coords[2]) + "\t")

                results_file.write(str(ligand_coords[0]) + "\t")
                results_file.write(str(ligand_coords[1]) + "\t")
                results_file.write(str(ligand_coords[2]) + "\t")

                results_file.write(str(comp_data["complex"]) + "\t")
                results_file.write(str(comp_data["summary"]) + "\t")
                results_file.write(str(comp_data["coordNo"]) + "\t")

                results_file.write(str(model_index) + "\n")

        results_file.close()

    def write_pi_pi_results(
        self, ligand, centroid, extracted_res, extracted_centroids, model_index
    ):
        """
        Write the results to the results file.
        """
        results_file_name = self.pi_pi_log
        ligand_code = ligand.get_resname()
        ligand_id = str(ligand.get_id()[1])
        ligand_chain = ligand.get_parent().get_id()
        results_file = open(results_file_name, "a+")
        new_atoms = []
        for res, cent in zip(extracted_res, extracted_centroids):
            distance = cent["distance"]
            theta = angle_between_norm_vec(centroid, cent)
            angle = angle_norm_vec_point(centroid, cent["coords"])
            omega = calc_omega(centroid, cent)

            h = abs(cos(radians(angle)) * distance)
            x = sin(radians(angle)) * distance

            if angle > 90.0:
                angle = 180 - angle

            if theta > 90.0:
                theta = 180 - theta

            if omega > 90.0:
                omega = 180 - omega

            centroid2_coords = cent["coords"]
            centroid_coords = centroid["coords"]

            residue_name = res.get_resname()
            res_chain = res.get_parent().get_id()
            res_id = str(res.get_id()[1])
            results_file.write(self.pdb_code + "\t")
            results_file.write(ligand_code + "\t")
            results_file.write(ligand_chain + "\t")
            results_file.write(ligand_id + "\t")
            results_file.write(residue_name + "\t")
            results_file.write(res_chain + "\t")
            results_file.write(res_id + "\t")

            results_file.write(str(distance) + "\t")
            results_file.write(str(angle) + "\t")

            results_file.write(str(x) + "\t")
            results_file.write(str(h) + "\t")

            results_file.write(str(theta) + "\t")
            results_file.write(str(omega) + "\t")

            results_file.write(str(centroid["cycleId"]) + "\t")
            results_file.write(str(centroid_coords[0]) + "\t")
            results_file.write(str(centroid_coords[1]) + "\t")
            results_file.write(str(centroid_coords[2]) + "\t")

            results_file.write(str(centroid2_coords[0]) + "\t")
            results_file.write(str(centroid2_coords[1]) + "\t")
            results_file.write(str(centroid2_coords[2]) + "\t")

            results_file.write(str(model_index) + "\t")
            results_file.write(str(cent["ringSize"]) + "\n")

        results_file.close()

        return new_atoms

    def write_anion_cation_results(
        self, anion_atom, ligand, centroid, extracted_cations, model_index
    ):
        """
        Write the results to the results file.
        """
        results_file_name = self.anion_cation_log
        anion = anion_atom.get_parent()
        anion_code = anion.get_resname()
        anion_id = str(anion.get_id()[1])
        anion_chain = anion.get_parent().get_id()
        anion_coord = anion_atom.get_coord()

        ligand_code = ligand.get_resname()
        ligand_id = str(ligand.get_id()[1])
        ligand_chain = ligand.get_parent().get_id()

        results_file = open(results_file_name, "a+")
        new_atoms = []

        anion_angle = atom_angle_nom_vec_centroid(anion_atom, centroid)
        first_semisphere_anion = anion_angle < 90

        for cat in extracted_cations:
            distance = anion_atom - cat

            cation_angle = atom_angle_nom_vec_centroid(cat, centroid)
            first_semisphere_cation = cation_angle < 90
            same_semisphere = first_semisphere_anion == first_semisphere_cation

            cat_coord = cat.get_coord()
            cat_res = cat.get_parent()

            residue_name = cat_res.get_resname()
            res_chain = cat_res.get_parent().get_id()
            res_id = str(cat_res.get_id()[1])
            results_file.write(self.pdb_code + "\t")
            results_file.write(residue_name + "\t")
            results_file.write(res_chain + "\t")
            results_file.write(res_id + "\t")

            results_file.write(ligand_code + "\t")
            results_file.write(ligand_chain + "\t")
            results_file.write(ligand_id + "\t")
            results_file.write(str(centroid["cycleId"]) + "\t")

            results_file.write(anion_code + "\t")
            results_file.write(anion_chain + "\t")
            results_file.write(anion_id + "\t")
            results_file.write(str(anion_atom.anionData.anion_id) + "\t")

            results_file.write(anion_atom.element + "\t")
            results_file.write(cat.element + "\t")

            results_file.write(str(distance) + "\t")

            results_file.write(str(anion_coord[0]) + "\t")
            results_file.write(str(anion_coord[1]) + "\t")
            results_file.write(str(anion_coord[2]) + "\t")

            results_file.write(str(cat_coord[0]) + "\t")
            results_file.write(str(cat_coord[1]) + "\t")
            results_file.write(str(cat_coord[2]) + "\t")

            results_file.write(str(same_semisphere) + "\t")
            results_file.write(str(abs(anion_angle - cation_angle)) + "\t")
            results_file.write(str(model_index) + "\n")

        results_file.close()

        return new_atoms

    def write_hbonds_results(self, h_donors, atom, model_index):
        results_file_name = self.h_bonds_log

        anion_atom = atom["Atom"]
        anion = anion_atom.get_parent()
        anion_code = anion.get_resname()
        anion_id = str(anion.get_id()[1])
        anion_chain = anion.get_parent().get_id()
        anion_coord = anion_atom.get_coord()

        results_file = open(results_file_name, "a+")

        for h_don_data in h_donors:
            h_don = h_don_data["donor"]
            distance = anion_atom - h_don

            h_don_coords = h_don.get_coord()
            h_don_res = h_don.get_parent()

            residue_name = h_don_res.get_resname()
            res_chain = h_don_res.get_parent().get_id()
            res_id = str(h_don_res.get_id()[1])

            results_file.write(self.pdb_code + "\t")

            results_file.write(anion_code + "\t")
            results_file.write(anion_chain + "\t")
            results_file.write(anion_id + "\t")

            results_file.write(residue_name + "\t")
            results_file.write(res_chain + "\t")
            results_file.write(res_id + "\t")

            results_file.write(atom["AnionType"] + "\t")
            results_file.write(anion_atom.element + "\t")
            results_file.write(str(anion_atom.anionData.anion_id) + "\t")

            results_file.write(str(anion_coord[0]) + "\t")
            results_file.write(str(anion_coord[1]) + "\t")
            results_file.write(str(anion_coord[2]) + "\t")

            results_file.write(h_don.element + "\t" + h_don.element + "\t")

            results_file.write(str(h_don_coords[0]) + "\t")
            results_file.write(str(h_don_coords[1]) + "\t")
            results_file.write(str(h_don_coords[2]) + "\t")

            hydrogen_atom = h_don_data["hydrogen"]
            hydrogen_coords = hydrogen_atom.get_coord()

            results_file.write(str(hydrogen_coords[0]) + "\t")
            results_file.write(str(hydrogen_coords[1]) + "\t")
            results_file.write(str(hydrogen_coords[2]) + "\t")

            results_file.write(str(h_don_data["HFromExp"]) + "\t")

            vec1 = normalize(anion_coord - hydrogen_coords)
            vec2 = normalize(h_don_coords - hydrogen_coords)
            angle = degrees(acos(np.inner(vec1, vec2)))
            h_distance = anion_atom - hydrogen_atom

            results_file.write(str(angle) + "\t")
            results_file.write(str(h_distance) + "\t")

            results_file.write(str(distance) + "\t")

            results_file.write(str(model_index) + "\n")

        results_file.close()

    def write_methyl_pi_results(
        self,
        ligand,
        centroid,
        extracted_atoms,
        model_index,
        resolution,
        method,
        structure_type,
    ):
        """
        Write the results to the results file.
        """
        results_file_name = self.methyl_pi_log
        ligand_code = ligand.get_resname()
        ligand_id = str(ligand.get_id()[1])
        ligand_chain = ligand.get_parent().get_id()
        results_file = open(results_file_name, "a+")
        new_atoms = []
        for atom_data in extracted_atoms:
            distance = atom_distance_from_centroid(atom_data["Atom"], centroid)
            angle = atom_angle_nom_vec_centroid(atom_data["Atom"], centroid)

            h = abs(cos(radians(angle)) * distance)
            x = sin(radians(angle)) * distance
            if angle > 90.0:
                angle = 180 - angle
            new_atoms.append(atom_data)

            atom_coords = atom_data["Atom"].get_coord()
            centroid_coords = centroid["coords"]

            anion = atom_data["Atom"].get_parent()
            residue_name = anion.get_resname()
            anion_chain = anion.get_parent().get_id()
            anion_id = str(anion.get_id()[1])
            results_file.write(self.pdb_code + "\t")
            results_file.write(ligand_code + "\t")
            results_file.write(ligand_chain + "\t")
            results_file.write(ligand_id + "\t")
            results_file.write(residue_name + "\t")
            results_file.write(anion_chain + "\t")
            results_file.write(anion_id + "\t")
            results_file.write(atom_data["AnionType"] + "\t")
            results_file.write(str(atom_data["AnionId"]) + "\t")
            results_file.write(atom_data["Atom"].element + "\t")

            results_file.write(str(distance) + "\t")
            results_file.write(str(angle) + "\t")

            results_file.write(str(x) + "\t")
            results_file.write(str(h) + "\t")

            results_file.write(str(centroid["cycleId"]) + "\t")
            results_file.write(str(centroid_coords[0]) + "\t")
            results_file.write(str(centroid_coords[1]) + "\t")
            results_file.write(str(centroid_coords[2]) + "\t")

            results_file.write(str(atom_coords[0]) + "\t")
            results_file.write(str(atom_coords[1]) + "\t")
            results_file.write(str(atom_coords[2]) + "\t")

            results_file.write(str(model_index) + "\t")
            results_file.write(
                str(atom_data["Atom"].get_parent().is_disordered()) + "\t"
            )
            results_file.write(str(centroid["ringSize"]) + "\t")
            results_file.write(str(centroid["ringElements"]) + "\t")
            results_file.write(str(resolution) + "\t")
            results_file.write(str(method) + "\t")
            results_file.write(str(structure_type) + "\n")

        results_file.close()

        return new_atoms

    def increment_partial_progress(self):
        file_name = self.partial_progress_log
        if not isfile(file_name):
            partial_progress_file = open(file_name, "w")
            partial_progress_file.write("1")
            partial_progress_file.close()

        else:
            partial_progress_file = open(file_name, "r")
            actual_no = int(partial_progress_file.readline())
            partial_progress_file.close()

            actual_no += 1
            partial_progress_file = open(file_name, "w")
            partial_progress_file.write(str(actual_no))
            partial_progress_file.close()


def atom_distance_from_centroid(atom, centroid):
    """Return the distance between an atom and a ring centroid.

    Args:
        atom: Biopython Atom.
        centroid: dict with keys "coords" and "normVec".

    Returns:
        Distance (float).
    """
    atom_coords = atom.get_coord()
    centorid_coords = centroid["coords"]

    dist = 0
    for atom_coord, centroid_coord in zip(atom_coords, centorid_coords):
        dist += (atom_coord - centroid_coord) * (atom_coord - centroid_coord)

    return sqrt(dist)


def atom_angle_nom_vec_centroid(atom, centroid):
    """Return the angle between the centroid-atom direction and the ring normal.

    Args:
        atom: Biopython Atom.
        centroid: dict with keys "coords" and "normVec".

    Returns:
        Angle in degrees.
    """
    atom_coords = np.array(atom.get_coord())
    centroid_coords = np.array(centroid["coords"])
    norm_vec = centroid["normVec"]

    centr_atom_vec = normalize(atom_coords - centroid_coords)
    inner_prod = np.inner(norm_vec, centr_atom_vec)

    return degrees(acos(inner_prod))


def calc_omega(pi_acid_centroid, pi_res_centroid):
    cent1cent2_vec = np.array(pi_res_centroid["coords"]) - np.array(
        pi_acid_centroid["coords"]
    )
    cent1cent2_vec = normalize(cent1cent2_vec)

    plane_norm_vec = np.cross(cent1cent2_vec, pi_acid_centroid["normVec"])
    plane_norm_vec = normalize(plane_norm_vec)

    return degrees(acos(np.inner(plane_norm_vec, pi_res_centroid["normVec"])))


def angle_between_norm_vec(centroid1, centroid2):
    norm_vec1 = centroid1["normVec"]
    norm_vec2 = centroid2["normVec"]

    inner_prod = np.inner(norm_vec1, norm_vec2)
    if abs(inner_prod) > 1.0:
        if abs(inner_prod) < 1.1:
            return 0.0
        else:
            return 666.0

    return degrees(acos(inner_prod))


def angle_norm_vec_point(centroid, point):
    coords = np.array(point)
    centroid_coords = np.array(centroid["coords"])
    norm_vec = centroid["normVec"]

    centr_atom_vec = normalize(coords - centroid_coords)
    inner_prod = np.inner(norm_vec, centr_atom_vec)

    return degrees(acos(inner_prod))


def calc_directional_vector(plane_data, centroid):
    vec_coords = []
    for point_data in plane_data.directional_vector:
        kind = list(point_data.keys())[0]
        if kind == "atom":
            vec_coords.append(point_data[kind].get_coord())
        elif kind == "center":
            atoms_list = point_data[kind]
            vec_coords.append(
                get_average_coords(atoms_list, list(range(len(atoms_list))))
            )
        elif kind == "closest":
            atoms_list = point_data[kind]

            min_dist = 1000
            closest_atom = None

            for atom in atoms_list:
                dist = atom_distance_from_centroid(atom, centroid)

                if dist < min_dist:
                    min_dist = dist
                    closest_atom = atom

            vec_coords.append(closest_atom.get_coord())

    vec = normalize(vec_coords[1] - vec_coords[0])
    vec2 = normalize(np.array(centroid["coords"]) - vec_coords[1])

    inner_prod = np.inner(vec, vec2)
    if abs(inner_prod) > 1.0:
        if abs(inner_prod) < 1.1:
            return 0.0
        else:
            return 666.0

    return degrees(acos(inner_prod))
