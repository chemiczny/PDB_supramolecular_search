"""
Created on Fri Aug  3 10:33:21 2018

@author: michal
"""

try:
    from pymol import cmd
except ImportError:
    pass
from supramolecular_gui import SupramolecularGUI
from simple_filters import (
    no_aa_in_pi_acids,
    no_aa_in_pi_res,
    no_nu_in_pi_acids,
    no_nu_in_pi_res,
)


class PiPiGUI(SupramolecularGUI):
    def __init__(self, page, parallel_selection_function, name):
        SupramolecularGUI.__init__(self, page, parallel_selection_function, name)
        self.set_numerical_parameters(
            {
                "R": {"header": "Distance"},
                "h": {"header": "h"},
                "x": {"header": "x"},
                "alpha": {"header": "Angle"},
                "theta": {"header": "theta"},
                "omega": {"header": "omega"},
            }
        )

        self.set_list_parameters(
            {
                "Pi 1": {"header": "Pi acid Code"},
                "Pi 2": {"header": "Pi res code"},
                "PDB": {"header": "PDB Code"},
            }
        )

        self.set_sorting_parameters(
            {
                "R": "Distance",
                "Angle": "Angle",
                "x": "x",
                "h": "h",
                "theta": "theta",
                "Pi 1": "Pi acid Code",
                "Pi 2": "Pi res code",
                "omega": "omega",
            },
            ["R", "Angle", "x", "h", "theta", "Pi 1", "Pi 2", "omega"],
        )

        self.set_tree_data(
            [
                "ID",
                "PDB",
                "Pi 1",
                "Pi 1 id",
                "Pi 2",
                "Pi 2 id",
                "R",
                "alpha",
                "x",
                "h",
                "theta",
                "omega",
            ]
        )

        self.set_additional_checkboxes(
            [
                {"label": "No AA in Pi 1", "func": no_aa_in_pi_acids},
                {"label": "No AA in Pi 2", "func": no_aa_in_pi_res},
                {"label": "No NU in Pi 1", "func": no_nu_in_pi_acids},
                {"label": "No NU in Pi 2", "func": no_nu_in_pi_res},
            ]
        )

        self.set_unique_parameters(
            {
                "PDB": ["PDB Code"],
                "Pi 1": ["Pi acid Code", "Pi acid chain", "Piacid id"],
                "Pi 2": ["Pi res code", "Pi res chain", "Pi res id"],
                "Model": ["Model No"],
            },
            ["PDB", "Pi 1", "Pi 2", "Model"],
        )

        self.interaction_buttons = [
            {
                "name": "All pi acid 1 int.",
                "headers": [
                    "PDB Code",
                    "Pi acid Code",
                    "Pi acid chain",
                    "Piacid id",
                    "Model No",
                ],
            },
            {
                "name": "All pi acid 2 int.",
                "headers": [
                    "PDB Code",
                    "Pi res code",
                    "Pi res chain",
                    "Pi res id",
                    "Model No",
                ],
            },
        ]

        self.arrow_name = "piPiArrow"
        self.arrow_color = "blue green"

    def get_values(self, row_id, row):
        return (
            row_id,
            row["PDB Code"],
            row["Pi acid Code"],
            row["Pi acid chain"] + str(row["Piacid id"]),
            row["Pi res code"],
            row["Pi res chain"] + str(row["Pi res id"]),
            str(row["Distance"])[:3],
            str(row["Angle"])[:4],
            str(row["x"])[:3],
            str(row["h"])[:3],
            str(row["theta"])[:4],
            str(row["omega"])[:4],
        )

    def get_selection(self, data):
        res1_id = data["Piacid id"].values[0]
        res1_chain = data["Pi acid chain"].values[0]

        res2_id = data["Pi res id"].values[0]
        res2_chain = data["Pi res chain"].values[0]
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

        res2_id = data["Pi res id"]
        res2_chain = data["Pi res chain"]
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
            data["Centroid 2 x coord"].values[0],
            data["Centroid 2 y coord"].values[0],
            data["Centroid 2 z coord"].values[0],
        ]

        return point1_coords, point2_coords

    def get_arrow_from_row(self, data):
        point1_coords = [
            data["Centroid x coord"],
            data["Centroid y coord"],
            data["Centroid z coord"],
        ]
        point2_coords = [
            data["Centroid 2 x coord"],
            data["Centroid 2 y coord"],
            data["Centroid 2 z coord"],
        ]

        return point1_coords, point2_coords
