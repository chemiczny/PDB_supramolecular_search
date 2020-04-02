from configure import configure
config = configure()

from Bio.PDB import PDBList
from os.path import dirname

pdb_dir = dirname(config["cif"])
pl = PDBList(pdb=pdb_dir)
pl.flat_tree = True
pl.update_pdb(file_format = "mmCif")
