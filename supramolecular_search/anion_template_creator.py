"""
Created on Wed May 23 12:34:08 2018

@author: michal
"""

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.readwrite.json_graph import node_link_data
from os.path import isdir, join, isfile, dirname, abspath
from os import mkdir
import json
from glob import glob
import shutil

ANION_TEMPLATES_DIR = join(dirname(abspath(__file__)), "anion_templates")


class AnionMatcher(GraphMatcher):
    def semantic_feasibility(self, g1_node, g2_node):
        if "charged" in self.G1.nodes[g1_node]:
            if self.G1.nodes[g1_node]["charged"] != self.G2.nodes[g2_node]["charged"]:
                return False
        elif self.G2.nodes[g2_node]["charged"]:
            return False

        if self.G2.nodes[g2_node]["terminating"]:
            if len(list(self.G2.neighbors(g2_node))) != len(
                list(self.G1.neighbors(g1_node))
            ):
                return False

        if (
            self.G2.nodes[g2_node]["element"] == "X"
            or "X" in self.G2.nodes[g2_node]["aliases"]
        ) and self.G1.nodes[g1_node]["element"] not in self.G2.nodes[g2_node][
            "notAliases"
        ]:
            return True

        return (
            self.G1.nodes[g1_node]["element"] == self.G2.nodes[g2_node]["element"]
            or self.G1.nodes[g1_node]["element"] in self.G2.nodes[g2_node]["aliases"]
        )


def add_atribute(graph, nodes, key):
    if isinstance(nodes, list):
        for node_id in nodes:
            graph.nodes[node_id][key] = True

    else:
        graph.nodes[nodes][key] = True


def save_anion(
    atoms,
    bonds,
    charged,
    name,
    priority,
    terminating=[],
    aliases={},
    not_aliases={},
    geometry={},
    full_isomorphism=False,
    name_mapping={},
    non_unique_charge=[],
    properties2measure=[],
):
    graph = nx.Graph()
    non_unique_charge = set(non_unique_charge)

    for i, el in enumerate(atoms):
        graph.add_node(
            i, element=el, terminating=False, bonded=False, aliases=[], charged=False
        )

    graph.add_edges_from(bonds)

    add_atribute(graph, terminating, "terminating")

    for node_id in aliases:
        graph.nodes[node_id]["aliases"] = aliases[node_id]

    for node_id in not_aliases:
        graph.nodes[node_id]["notAliases"] = not_aliases[node_id]

    if not geometry:
        graph.graph["geometry"] = "no restrictions"
    else:
        graph.graph["geometry"] = geometry

    graph.graph["fullIsomorphism"] = full_isomorphism
    graph.graph["name"] = name
    graph.graph["nameMapping"] = name_mapping
    graph.graph["priority"] = priority
    graph.graph["properties2measure"] = properties2measure

    file_name = str(priority) + "_" + name

    if isinstance(charged, list):
        unique_charges = set(charged)

        for node_id in charged:
            nuc = unique_charges | non_unique_charge
            nuc.remove(node_id)
            save_anion_json(graph, file_name, node_id, nuc)
    else:
        save_anion_json(graph, file_name, charged, non_unique_charge)


def save_anion_json(graph, file_name, charged, non_unique_charges=[]):
    main_element = graph.nodes[charged]["element"]
    elements = [main_element]

    if "aliases" in graph.nodes[charged]:
        elements += graph.nodes[charged]["aliases"]
        graph.nodes[charged]["aliases"] = []

    graph.nodes[charged]["charged"] = True
    graph.graph["charged"] = charged
    graph.graph["otherCharges"] = list(non_unique_charges)

    old_name = ""
    name_mapping = False
    if "X" in graph.graph["name"] and charged in graph.graph["nameMapping"]:
        old_name = graph.graph["name"]
        name_mapping = graph.graph["nameMapping"][charged]
        graph.graph["nameMapping"].pop(charged)

    for element in elements:
        graph.nodes[charged]["element"] = element

        if name_mapping:
            graph.graph["name"] = old_name.replace(name_mapping, element)

        dir_path = join(ANION_TEMPLATES_DIR, element)
        if not isdir(dir_path):
            mkdir(dir_path)

        path2save = get_unique_path(dir_path, file_name)
        output = open(path2save, "w")

        json.dump(node_link_data(graph, edges="links"), output)
        output.close()

    graph.nodes[charged]["charged"] = False


def get_unique_path(dir_path, file_name):
    path2save = join(dir_path, file_name + ".json")
    if not isfile(path2save):
        return path2save

    similar_files = glob(join(dir_path, file_name) + "_*.json")
    if not similar_files:
        return join(dir_path, file_name + "_0.json")

    max_number = -1
    for s in similar_files:
        new_number = int(s[:-5].split("_")[-1])
        max_number = max(max_number, new_number)

    return join(dir_path, file_name + "_" + str(max_number + 1) + ".json")


def clear_anion_templates():
    if isdir(ANION_TEMPLATES_DIR):
        shutil.rmtree(ANION_TEMPLATES_DIR)
    mkdir(ANION_TEMPLATES_DIR)


if __name__ == "__main__":
    clear_anion_templates()

    # OXYGEN

    #    #RCOOH
    save_anion(
        ["C", "C", "O", "O"],
        [(0, 1), (1, 2), (1, 3)],
        2,
        "RCOO",
        0,
        terminating=[1, 2, 3],
        geometry="planar",
        non_unique_charge=[3],
        properties2measure=[
            {
                "kind": "plane",
                "atoms": [1, 2, 3],
                "directionalVector": [{"atom": 1}, {"center": [2, 3]}],
            }
        ],
    )

    # ClO, BrO, IO,
    save_anion(
        ["CL", "O"],
        [(0, 1)],
        1,
        "XO",
        5,
        full_isomorphism=True,
        aliases={0: ["BR", "I"]},
        name_mapping={0: "X"},
        properties2measure=[{"kind": "line", "atoms": [0, 1]}],
    )

    # NO2, ClO2, BRO2,
    save_anion(
        ["N", "O", "O"],
        [(0, 1), (0, 2)],
        1,
        "XO2",
        10,
        full_isomorphism=True,
        aliases={0: ["CL", "BR"]},
        name_mapping={0: "X"},
        non_unique_charge=[2],
        properties2measure=[
            {
                "kind": "plane",
                "atoms": [0, 1, 2],
                "directionalVector": [{"atom": 0}, {"center": [1, 2]}],
            }
        ],
    )

    # NO3, CO3, PO3, SO3, AsO3, BO3, ClO3, BRO3
    save_anion(
        ["N", "O", "O", "O"],
        [(0, 1), (0, 2), (0, 3)],
        1,
        "XO3",
        15,
        full_isomorphism=True,
        aliases={0: ["C", "P", "B", "S", "AS", "CL", "BR", "I"]},
        name_mapping={0: "X"},
        non_unique_charge=[2, 3],
        properties2measure=[
            {
                "kind": "plane",
                "atoms": [1, 2, 3],
                "directionalVector": [{"closest": [1, 2, 3]}, {"center": [1, 2, 3]}],
            }
        ],
    )

    # PO4, SO4, AsO4, ClO4, BRO4
    save_anion(
        ["P", "O", "O", "O", "O"],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
        1,
        "XO4",
        20,
        full_isomorphism=True,
        aliases={0: ["S", "AS", "CL", "BR", "I"]},
        name_mapping={0: "X"},
        non_unique_charge=[2, 3, 4],
    )

    #   Ph-OH
    # save_anion(
    #     ["C", "C", "C", "C", "C", "C", "O"],
    #     [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (5, 6)],
    #     6, "PhOH", 25, terminating=[6], geometry="planarWithSubstituents",
    # )

    #    #RBOOH
    save_anion(
        ["X", "B", "O", "O"],
        [(0, 1), (1, 2), (1, 3)],
        2,
        "RBOO",
        30,
        terminating=[2, 3],
        not_aliases={0: ["O"]},
        non_unique_charge=[3],
        properties2measure=[
            {
                "kind": "plane",
                "atoms": [1, 2, 3],
                "directionalVector": [{"atom": 1}, {"center": [2, 3]}],
            }
        ],
    )

    # COO
    save_anion(
        ["C", "O", "O"],
        [(0, 1), (0, 2)],
        1,
        "COO",
        35,
        terminating=[1, 2],
        non_unique_charge=[2],
        properties2measure=[
            {
                "kind": "plane",
                "atoms": [0, 1, 2],
                "directionalVector": [{"atom": 0}, {"center": [1, 2]}],
            }
        ],
    )

    # R-PO4, R-SO4, R-AsO4
    save_anion(
        ["P", "O", "O", "O", "O"],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
        1,
        "R-XO4",
        45,
        terminating=[1, 2, 3],
        aliases={0: ["S", "AS"]},
        name_mapping={0: "X"},
        non_unique_charge=[2, 3],
    )

    # R2-PO4, R2-SO4, R2-AsO4
    save_anion(
        ["P", "O", "O", "O", "O"],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
        1,
        "R2-XO4",
        47,
        terminating=[1, 2],
        aliases={0: ["S", "AS"]},
        name_mapping={0: "X"},
        non_unique_charge=[2],
    )

    # R3-PO4, R3-SO4, R3-AsO4
    #    saveAnion( ["P", "O", "O", "O", "O"], [(0,1), (0,2), (0,3), (0, 4)],
    #              1, "R2-XO4", 48, terminating = [ 1 ] ,
    #              aliases = { 0 : [ "S", "AS" ] }, nameMapping = { 0 : "X" } )

    # RAsO3, RPO3, RSO3
    save_anion(
        ["P", "O", "O", "O", "C"],
        [(0, 1), (0, 2), (0, 3), (0, 4)],
        1,
        "RXO3",
        50,
        terminating=[1, 2, 3],
        aliases={0: ["S", "AS"]},
        name_mapping={0: "X"},
        non_unique_charge=[2, 3],
    )

    # R2AsO2, R2PO2, RRSO2
    #    saveAnion( ["P", "O", "O", "C", "C"], [(0,1), (0,2), (0,3), (0, 4)],
    #              1, "R2XO2", 55, terminating = [1, 2],
    #              aliases = { 0 : [ "S", "AS" ] }, nameMapping = { 0 : "X" } )

    # F, CL, BR, I, S
    save_anion(
        ["F"],
        [],
        0,
        "X",
        55,
        aliases={0: ["CL", "BR", "I", "S"]},
        full_isomorphism=True,
        name_mapping={0: "X"},
    )

    # SCN
    save_anion(
        ["S", "C", "N"],
        [(0, 1), (0, 2)],
        [0, 1, 2],
        "SCN",
        62,
        full_isomorphism=True,
        properties2measure=[{"kind": "line", "atoms": [0, 2]}],
    )

    #    #RSH
    #    saveAnion( [ "X" , "S" ], [ (0,1)],
    #              1, "RSH", 60, terminating = [1],
    #              notAliases = {0 : [ "O" ] } )
    #

    # N3
    save_anion(
        ["N", "N", "N"],
        [(0, 1), (0, 2)],
        [0, 1],
        "N3",
        70,
        full_isomorphism=True,
        non_unique_charge=[2],
        properties2measure=[{"kind": "lineSymmetric", "atoms": [0, 2]}],
    )

    # CN
    save_anion(
        ["C", "N"],
        [(0, 1)],
        [0, 1],
        "CN",
        75,
        full_isomorphism=True,
        properties2measure=[{"kind": "line", "atoms": [0, 1]}],
    )


#    #RSSR
#    saveAnion( [ "X" , "S", "S" ], [ (0,1), (1,2)],
#              1, "RSS", 80 ,
#              notAliases = {0 : [ "O" ] } )
