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
    no_aa_in_anions,
    no_nu_in_anions,
    no_nu_in_h_donors,
    no_aa_in_h_donors,
)


class HBondsGUI(SupramolecularGUI):
    def __init__(self, page, parallel_selection_function, name):
        SupramolecularGUI.__init__(self, page, parallel_selection_function, name)
        self.set_numerical_parameters(
            {
                "R": {"header": "Distance Don Acc"},
                "RH": {"header": "Distance H Acc"},
                "Angle": {"header": "Angle"},
            }
        )

        self.set_list_parameters(
            {
                "Acceptor": {"header": "Anion code"},
                "Acceptor el.": {"header": "Acceptor atom"},
                "Donor": {"header": "Donor code"},
                "Donor el.": {"header": "Donor atom"},
                "PDB": {"header": "PDB Code"},
                "H orig.": {"header": "H from Experm"},
            }
        )

        self.set_sorting_parameters(
            {
                "R": "Distance Don Acc",
                "Donor": "Donor code",
                "Don. el.": "Donor atom",
                "Acceptor": "Anion code",
                "Acc. el.": "Acceptor atom",
                "RH": "Distance H Acc",
                "Angle": "Angle",
            },
            ["R", "Donor", "Don. el.", "Acceptor", "Acc. el.", "RH", "Angle"],
        )

        self.set_tree_data(
            [
                "ID",
                "PDB",
                "Acceptor",
                "Anion id",
                "Acc. el.",
                "Anion gr. id",
                "Donor",
                "Donor id",
                "Donor el.",
                "R",
                "RH",
                "Angle",
                "H orig",
            ]
        )

        self.set_additional_checkboxes(
            [
                {"label": "No AA in acceptors", "func": no_aa_in_anions},
                {"label": "No AA in donors", "func": no_aa_in_h_donors},
                {"label": "No NU in acceptors", "func": no_nu_in_anions},
                {"label": "No NU in donors", "func": no_nu_in_h_donors},
            ]
        )

        self.set_unique_parameters(
            {
                "PDB": ["PDB Code"],
                "Acceptor": ["Anion code", "Anion chain", "Anion id"],
                "Acceptor id": ["Anion group id"],
                "Donor": ["Donor code", "Donor chain", "Donor id"],
                "Model": ["Model No"],
            },
            ["PDB", "Acceptor", "Acceptor id", "Donor", "Model"],
        )

        self.interaction_buttons = [
            {
                "name": "All acc int.",
                "headers": [
                    "PDB Code",
                    "Anion code",
                    "Anion chain",
                    "Anion id",
                    "Model No",
                ],
            },
            {
                "name": "All don int.",
                "headers": [
                    "PDB Code",
                    "Donor code",
                    "Donor chain",
                    "Donor id",
                    "Model No",
                ],
            },
        ]

        self.arrow_name = "HBondArrow"
        self.arrow_color = "red violet"

    def get_values(self, row_id, row):
        return (
            row_id,
            row["PDB Code"],
            row["Anion code"],
            row["Anion chain"] + str(row["Anion id"]),
            row["Acceptor atom"],
            row["Anion group id"],
            row["Donor code"],
            row["Donor chain"] + str(row["Donor id"]),
            row["Donor atom"],
            str(row["Distance Don Acc"])[:3],
            str(row["Distance H Acc"])[:3],
            str(row["Angle"]),
            str(row["H from Experm"]),
        )

    def get_selection(self, data):
        res1_id = data["Anion id"].values[0]
        res1_chain = data["Anion chain"].values[0]

        res2_id = data["Donor id"].values[0]
        res2_chain = data["Donor chain"].values[0]
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

        res2_id = data["Donor id"]
        res2_chain = data["Donor chain"]
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
            data["Acceptor x coord"].values[0],
            data["Acceptor y coord"].values[0],
            data["Acceptor z coord"].values[0],
        ]
        point2_coords = [
            data["Hydrogen x coord"].values[0],
            data["Hydrogen y coord"].values[0],
            data["Hydrogen z coord"].values[0],
        ]

        return point1_coords, point2_coords

    def get_arrow_from_row(self, data):
        point1_coords = [
            data["Acceptor x coord"],
            data["Acceptor y coord"],
            data["Acceptor z coord"],
        ]
        point2_coords = [
            data["Hydrogen x coord"],
            data["Hydrogen y coord"],
            data["Hydrogen z coord"],
        ]

        return point1_coords, point2_coords
