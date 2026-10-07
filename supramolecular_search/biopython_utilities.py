"""
Created on Fri Jun 29 12:53:24 2018

@author: michal

In older versions of Biopython method get_fullid() returns empty tuple
"""


def create_res_id(residue):
    id_no = residue.get_id()[1]
    chain = residue.get_parent().get_id()
    name = residue.get_resname()
    first_coord = list(residue.get_atoms())[0].get_coord()
    coord_str = str(first_coord[0]) + str(first_coord[1]) + str(first_coord[2])

    return chain + str(id_no) + name + coord_str


def create_res_id_from_atom(atom):
    return create_res_id(atom.get_parent())
