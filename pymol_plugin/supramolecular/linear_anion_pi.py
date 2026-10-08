"""
Created on Fri Aug  3 10:32:43 2018

@author: michal
"""

try:
    from pymol import cmd
except ImportError:
    pass
from supramolecular_gui import SupramolecularGUI
from simple_filters import (
    no_aa_in_pi_acids,
    no_aa_in_anions,
    no_nu_in_pi_acids,
    no_nu_in_anions,
)


class LinearAnionPiGUI(SupramolecularGUI):
    def __init__(self, page, parallel_selection_function, name):
        SupramolecularGUI.__init__(self, page, parallel_selection_function, name)
        self.set_numerical_parameters({"alpha": {"header": "Angle"}})

        self.set_list_parameters(
            {
                "Pi acid": {"header": "Pi acid Code"},
                "Anions": {"header": "Anion code"},
                "PDB": {"header": "PDB Code"},
            }
        )

        self.set_sorting_parameters(
            {"Angle": "Angle", "Pi acid": "Pi acid Code", "Anion": "Anion code"},
            ["Angle", "Pi acid", "Anion"],
        )

        self.set_tree_data(
            [
                "ID",
                "PDB",
                "Pi acid",
                "Pi acid id",
                "Anion",
                "Anion id",
                "Anion gr. id",
                "alpha",
            ]
        )

        self.set_additional_checkboxes(
            [
                {"label": "No AA in Pi acids", "func": no_aa_in_pi_acids},
                {"label": "No AA in anions", "func": no_aa_in_anions},
                {"label": "No NU in Pi acids", "func": no_nu_in_pi_acids},
                {"label": "No NU in anions", "func": no_nu_in_anions},
            ]
        )

        self.set_unique_parameters(
            {
                "PDB": ["PDB Code"],
                "Pi acid": ["Pi acid Code", "Pi acid chain", "Piacid id"],
                "Anion": ["Anion code", "Anion chain", "Anion id"],
                "Anion id": ["Anion group id"],
                "Ring id": ["CentroidId"],
                "Model": ["Model No"],
            },
            ["PDB", "Pi acid", "Ring id", "Anion", "Anion id", "Model"],
        )

        self.interaction_buttons = [
            {
                "name": "All anion int.",
                "headers": [
                    "PDB Code",
                    "Anion code",
                    "Anion chain",
                    "Anion id",
                    "Model No",
                ],
            },
            {
                "name": "All pi acid int.",
                "headers": [
                    "PDB Code",
                    "Pi acid Code",
                    "Pi acid chain",
                    "Piacid id",
                    "Model No",
                ],
            },
        ]

        self.arrow_name = "linearAnionPiArrow"
        self.arrow_color = "blue red"

    def get_values(self, row_id, row):
        return (
            row_id,
            row["PDB Code"],
            row["Pi acid Code"],
            row["Pi acid chain"] + str(row["Piacid id"]),
            row["Anion code"],
            row["Anion chain"] + str(row["Anion id"]),
            row["Anion group id"],
            str(row["Angle"])[:4],
        )

    def get_selection(self, data):
        res1_id = data["Piacid id"].values[0]
        res1_chain = data["Pi acid chain"].values[0]

        res2_id = data["Anion id"].values[0]
        res2_chain = data["Anion chain"].values[0]
        cmd.hide("everything")
        selection = (
            "( "
            + "chain "
            + res1_chain
            + " and resi "
            + str(res1_id)
            + " ) or ( "
            + " chain "
            + res2_chain
            + " and resi "
            + str(res2_id)
            + ")"
        )
        return selection

    def get_selection_from_row(self, data):
        res1_id = data["Piacid id"]
        res1_chain = data["Pi acid chain"]

        res2_id = data["Anion id"]
        res2_chain = data["Anion chain"]
        selection = (
            "( "
            + "chain "
            + res1_chain
            + " and resi "
            + str(res1_id)
            + " ) or ( "
            + " chain "
            + res2_chain
            + " and resi "
            + str(res2_id)
            + ")"
        )
        return selection

    def get_arrow(self, data):
        centroid_coords = [
            data["Centroid x coord"].values[0],
            data["Centroid y coord"].values[0],
            data["Centroid z coord"].values[0],
        ]
        anion_atom_coords = [
            data["Anion group x coord"].values[0],
            data["Anion group y coord"].values[0],
            data["Anion group z coord"].values[0],
        ]

        return centroid_coords, anion_atom_coords

    def get_arrow_from_row(self, data):
        centroid_coords = [
            data["Centroid x coord"],
            data["Centroid y coord"],
            data["Centroid z coord"],
        ]
        anion_atom_coords = [
            data["Anion group x coord"],
            data["Anion group y coord"],
            data["Anion group z coord"],
        ]

        return centroid_coords, anion_atom_coords
