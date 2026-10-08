"""
Created on Fri Aug  3 10:33:21 2018

@author: michal
"""

try:
    from pymol import cmd
except ImportError:
    pass
from supramolecular_gui import SupramolecularGUI
from simple_filters import only_anions, only_complexes


class MetalLigandGUI(SupramolecularGUI):
    def __init__(self, page, parallel_selection_function, name):
        SupramolecularGUI.__init__(self, page, parallel_selection_function, name)
        self.set_numerical_parameters({"R": {"header": "Distance"}})

        self.set_list_parameters(
            {
                "Cation": {"header": "Cation code"},
                "Cat. el.": {"header": "Cation element"},
                "Ligand": {"header": "Anion code"},
                "Lig. el..": {"header": "Ligand element"},
                "Anion type": {"header": "anionType"},
                "coordNo": {"header": "CoordNo"},
                "summary": {"header": "Summary"},
                "PDB": {"header": "PDB Code"},
            }
        )

        self.set_sorting_parameters(
            {
                "R": "Distance",
                "Cation": "Cation code",
                "Cat. el.": "Cation element",
                "Ligand": "Anion code",
                "Lig. el.": "Ligand element",
                "Anion type": "anionType",
            },
            ["R", "Cation", "Cat. el.", "Ligand", "Lig. el.", "Anion type"],
        )

        self.set_tree_data(
            [
                "ID",
                "PDB",
                "Cation",
                "Cation id",
                "Cat. el.",
                "Ligand",
                "Ligand id",
                "Ligand el.",
                "Is Anion",
                "Anion type",
                "Anion gr. id",
                "R",
                "Complex",
                "Summary",
                "Coord No",
            ]
        )

        self.set_additional_checkboxes(
            [
                {"label": "Only anions as ligands", "func": only_anions},
                {"label": "Only complexes", "func": only_complexes},
            ]
        )

        self.set_unique_parameters(
            {
                "PDB": ["PDB Code"],
                "Metal": ["Cation Code", "Cation chain", "Cation id"],
                "Ligand": ["Anion code", "Anion chain", "Anion id"],
                "Model": ["Model No"],
            },
            ["PDB", "Metal", "Ligand", "Model"],
        )

        self.interaction_buttons = [
            {
                "name": "All metal int.",
                "headers": [
                    "PDB Code",
                    "Cation code",
                    "Cation chain",
                    "Cation id",
                    "Model No",
                ],
            }
        ]

        self.arrow_name = "MetalLigandArrow"
        self.arrow_color = "red orange"
        self.num_of_sorting_menu = 1

    def get_values(self, row_id, row):
        return (
            row_id,
            row["PDB Code"],
            row["Cation code"],
            row["Cation chain"] + str(row["Cation id"]),
            row["Cation element"],
            row["Anion code"],
            row["Anion chain"] + str(row["Anion id"]),
            row["Ligand element"],
            str(row["isAnion"]),
            row["anionType"],
            row["Anion group id"],
            str(row["Distance"])[:3],
            str(row["Complex"]),
            str(row["Summary"]),
            str(row["CoordNo"]),
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
