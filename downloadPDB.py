import sys
sys.path.insert(0, "/net/archive/groups/plggsuprm/pythonPackages/" )

from Bio.PDB import PDBList

configurationFileName = "config.json"
configFile = open(configurationFileName)
config = json.load(configFile)
configFile.close()

pdb_dir = dirname(config["cif"])
pl = PDBList(pdb=pdb_dir)
pl.flat_tree = True
pl.download_entire_pdb(file_format = "mmCif")
