"""
Created on Fri Aug  3 10:32:57 2018

@author: michal
"""

try:
    from pymol import cmd
except ImportError:
    pass
from supramolecular_gui import SupramolecularGUI
from simple_filters import no_aa_in_pi_acids, no_nu_in_pi_acids, no_aa_in_cations


class CationPiGUI(SupramolecularGUI):
    def __init__(self, page, parallel_selection_function, name):
        SupramolecularGUI.__init__(self, page, parallel_selection_function, name)
        self.set_numerical_parameters(
            {
                "R": {"header": "Distance"},
                "h": {"header": "h"},
                "x": {"header": "x"},
                "alpha": {"header": "Angle"},
                "Chain size": {"header": "RingChain"},
            }
        )

        self.set_list_parameters(
            {
                "Pi acid": {"header": "Pi acid Code"},
                "Cation": {"header": "Cation code"},
                "Element": {"header": "Atom symbol"},
                "Chain": {"header": "RingChain"},
                "PDB": {"header": "PDB Code"},
            }
        )

        self.set_sorting_parameters(
            {
                "R": "Distance",
                "Angle": "Angle",
                "x": "x",
                "h": "h",
                "Pi acid": "Pi acid Code",
                "Cation": "Cation code",
                "Cat. el.": "Atom symbol",
            },
            ["R", "Angle", "x", "h", "Pi acid", "Cation", "Cat. el."],
        )

        self.set_tree_data(
            [
                "ID",
                "PDB",
                "Pi acid",
                "Pi acid id",
                "Cation",
                "Cation id",
                "Cat. el.",
                "R",
                "alpha",
                "x",
                "h",
                "chain",
                "chain Flat",
            ]
        )

        self.set_additional_checkboxes(
            [
                {"label": "No AA in Pi acids", "func": no_aa_in_pi_acids},
                {"label": "No NU in Pi acids", "func": no_nu_in_pi_acids},
                {"label": "No AA in cations", "func": no_aa_in_cations},
            ]
        )

        self.set_unique_parameters(
            {
                "PDB": ["PDB Code"],
                "Pi acid": ["Pi acid Code", "Pi acid chain", "Piacid id"],
                "Cation": ["Cation code", "Cation chain", "Cation id"],
                "Model": ["Model No"],
            },
            ["PDB", "Pi acid", "Cation", "Model"],
        )

        self.interaction_buttons = [
            {
                "name": "All cat. int.",
                "headers": [
                    "PDB Code",
                    "Cation code",
                    "Cation chain",
                    "Cation id",
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

        self.arrow_name = "cationPiArrow"
        self.arrow_color = "blue orange"

    def get_values(self, row_id, row):
        return (
            row_id,
            row["PDB Code"],
            row["Pi acid Code"],
            row["Pi acid chain"] + str(row["Piacid id"]),
            row["Cation code"],
            row["Cation chain"] + str(row["Cation id"]),
            row["Atom symbol"],
            str(row["Distance"])[:3],
            str(row["Angle"])[:4],
            str(row["x"])[:3],
            str(row["h"])[:3],
            str(row["RingChain"]),
            str(row["ChainFlat"]),
        )

    def get_selection(self, data):
        res1_id = data["Piacid id"].values[0]
        res1_chain = data["Pi acid chain"].values[0]

        res2_id = data["Cation id"].values[0]
        res2_chain = data["Cation chain"].values[0]
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

        res2_id = data["Cation id"]
        res2_chain = data["Cation chain"]
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
        point1_coords = [
            data["Centroid x coord"].values[0],
            data["Centroid y coord"].values[0],
            data["Centroid z coord"].values[0],
        ]
        point2_coords = [
            data["Cation x coord"].values[0],
            data["Cation y coord"].values[0],
            data["Cation z coord"].values[0],
        ]

        return point1_coords, point2_coords

    def get_arrow_from_row(self, data):
        point1_coords = [
            data["Centroid x coord"],
            data["Centroid y coord"],
            data["Centroid z coord"],
        ]
        point2_coords = [
            data["Cation x coord"],
            data["Cation y coord"],
            data["Cation z coord"],
        ]

        return point1_coords, point2_coords
