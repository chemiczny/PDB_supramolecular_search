"""
Created on Fri Aug  3 10:33:21 2018

@author: michal
"""

try:
    from pymol import cmd
except ImportError:
    pass

from supramolecular_gui import SupramolecularGUI
from simple_filters import no_aa_in_anions, no_nu_in_anions, no_aa_in_cations


class AnionCationGUI(SupramolecularGUI):
    def __init__(self, page, parallel_selection_function, name):
        SupramolecularGUI.__init__(self, page, parallel_selection_function, name)
        self.set_numerical_parameters(
            {"R": {"header": "Distance"}, "Lat dif": {"header": "Latitude diff"}}
        )

        self.set_list_parameters(
            {
                "Cation": {"header": "Cation code"},
                "Cat. el.": {"header": "Cation symbol"},
                "Anion": {"header": "Anion code"},
                "An. at.": {"header": "Anion symbol"},
                "Pi acid": {"header": "Pi acid Code"},
                "Same hem.": {"header": "Same semisphere"},
                "PDB": {"header": "PDB Code"},
            }
        )

        self.set_sorting_parameters(
            {
                "R": "Distance",
                "Cation": "Cation code",
                "Cat. el.": "Cation symbol",
                "Anion": "Anion code",
                "An. el.": "Anion symbol",
                "Lat dif": "Latitude diff",
            },
            ["R", "Cation", "Cat. el.", "Anion", "An. el.", "Lat dif"],
        )

        self.set_tree_data(
            [
                "ID",
                "PDB",
                "Cation",
                "Cation id",
                "Cat. el.",
                "Anion",
                "Anion id",
                "Anion el.",
                "Anion gr. id",
                "Pi acid",
                "Pi acid id",
                "R",
                "Same semisphere",
                "Lat dif",
            ]
        )

        self.set_additional_checkboxes(
            [
                {"label": "No AA in anions", "func": no_aa_in_anions},
                {"label": "No NU in anions", "func": no_nu_in_anions},
                {"label": "No AA in cations", "func": no_aa_in_cations},
            ]
        )

        self.set_unique_parameters(
            {
                "PDB": ["PDB Code"],
                "Cation": ["Cation code", "Cation chain", "Cation id"],
                "Anion": ["Anion code", "Anion chain", "Anion id"],
                "Model": ["Model No"],
            },
            ["PDB", "Cation", "Anion", "Model"],
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
                "name": "All anion int.",
                "headers": [
                    "PDB Code",
                    "Anion code",
                    "Anion chain",
                    "Anion id",
                    "Model No",
                ],
            },
        ]

        self.arrow_name = "AnionCationArrow"
        self.arrow_color = "red orange"

    def get_values(self, row_id, row):
        return (
            row_id,
            row["PDB Code"],
            row["Cation code"],
            row["Cation chain"] + str(row["Cation id"]),
            row["Cation symbol"],
            row["Anion code"],
            row["Anion chain"] + str(row["Anion id"]),
            row["Anion symbol"],
            row["Anion group id"],
            row["Pi acid Code"],
            row["Pi acid chain"] + str(row["Piacid id"]),
            str(row["Distance"])[:3],
            str(row["Same semisphere"]),
            str(row["Latitude diff"]),
        )

    def get_selection(self, data):
        res1_id = data["Anion id"].values[0]
        res1_chain = data["Anion chain"].values[0]

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
        res1_id = data["Anion id"]
        res1_chain = data["Anion chain"]

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
            data["Anion x coord"].values[0],
            data["Anion y coord"].values[0],
            data["Anion z coord"].values[0],
        ]
        point2_coords = [
            data["Cation x coord"].values[0],
            data["Cation y coord"].values[0],
            data["Cation z coord"].values[0],
        ]

        return point1_coords, point2_coords

    def get_arrow_from_row(self, data):
        point1_coords = [
            data["Anion x coord"],
            data["Anion y coord"],
            data["Anion z coord"],
        ]
        point2_coords = [
            data["Cation x coord"],
            data["Cation y coord"],
            data["Cation z coord"],
        ]

        return point1_coords, point2_coords
